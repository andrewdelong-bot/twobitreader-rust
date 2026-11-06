// A specialized stand-in for `HashMap<String, usize>`, mapping sequence names to indices into
// `TwobitReader::seqs`.
//
// Human genomes are accessed almost entirely by the primary names chr1..chr22, chrX, chrY and chrM, so
// those names are mapped to an array slot by inspecting the characters, with no hashing or comparison.
// Every other name (alternate contigs such as chr1_KI270706v1_random, or non-human names) falls back to
// a HashMap.

use std::collections::HashMap;

// The slots are the chromosome numbers 1..=22 (slot 0 is unused), then chrX, chrY and chrM.
const SLOT_X: usize = 23;
const SLOT_Y: usize = 24;
const SLOT_M: usize = 25;
const NUM_SLOTS: usize = 26;

// Marks a slot that no sequence maps to.
const EMPTY: usize = usize::MAX;

// Returns the slot for a primary sequence name, or None for any other name.
// Each name maps to its own slot, and no two names share one.
fn direct_slot(name: &str) -> Option<usize> {
    match name.strip_prefix("chr")?.as_bytes() {
        [digit @ b'1'..=b'9'] => Some(usize::from(digit - b'0')), // chr1..chr9
        [tens @ b'1'..=b'2', ones @ b'0'..=b'9'] => {
            let number = usize::from(tens - b'0') * 10 + usize::from(ones - b'0');
            (number <= 22).then_some(number) // chr10..chr22, but not chr23..chr29
        }
        b"X" => Some(SLOT_X),
        b"Y" => Some(SLOT_Y),
        b"M" => Some(SLOT_M),
        _ => None,
    }
}

#[derive(Debug)]
pub(crate) struct SeqIndexMap {
    direct: Vec<usize>,            // direct[direct_slot(name)] is the sequence index, or EMPTY.
    other: HashMap<String, usize>, // The sequence index for each name that has no direct slot.
}

impl SeqIndexMap {
    // Returns the sequence index for the name, if there is one.
    pub(crate) fn get(&self, name: &str) -> Option<&usize> {
        match direct_slot(name) {
            Some(slot) => Some(&self.direct[slot]).filter(|&&index| index != EMPTY),
            None => self.other.get(name),
        }
    }

    pub(crate) fn contains_key(&self, name: &str) -> bool {
        self.get(name).is_some()
    }

    // Returns the number of distinct names.
    pub(crate) fn len(&self) -> usize {
        self.direct.iter().filter(|&&index| index != EMPTY).count() + self.other.len()
    }
}

// As with HashMap, when a name occurs more than once the last index wins and len() counts it once.
impl FromIterator<(String, usize)> for SeqIndexMap {
    fn from_iter<I: IntoIterator<Item = (String, usize)>>(iter: I) -> Self {
        let mut map = Self { direct: vec![EMPTY; NUM_SLOTS], other: HashMap::new() };
        for (name, index) in iter {
            match direct_slot(&name) {
                Some(slot) => map.direct[slot] = index,
                None => {
                    map.other.insert(name, index);
                }
            }
        }
        map
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn primary_names() -> Vec<String> {
        let numbered = (1..=22).map(|n| format!("chr{n}"));
        numbered.chain(["chrX", "chrY", "chrM"].map(String::from)).collect()
    }

    #[test]
    fn primary_names_have_distinct_slots() {
        let slots: Vec<_> = primary_names().iter().map(|name| direct_slot(name).unwrap()).collect();
        assert_eq!(slots, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]);
    }

    #[test]
    fn other_names_have_no_slot() {
        let names = [
            "",
            "chr",
            "chr0",
            "chr00",
            "chr01",
            "chr23",
            "chr29",
            "chr30",
            "chr100",
            "chrx",
            "chry",
            "chrm",
            "chrXY",
            "chrMT",
            "Chr1",
            "CHR1",
            "1",
            "X",
            "chr1_KI270706v1_random",
            "chrUn_GL000195v1",
            "chrEBV",
            "seq_abc",
        ];
        for name in names {
            assert_eq!(direct_slot(name), None, "{name:?}");
        }
    }

    // Behaves like a HashMap, whether the names have direct slots or not.
    #[test]
    fn matches_hashmap() {
        let mut names = primary_names();
        names.extend(
            ["chr1_KI270706v1_random", "chrUn_GL000195v1", "chrEBV", "seq_abc", "chr23", "chr01", ""].map(String::from),
        );
        let pairs: Vec<_> = names.iter().cloned().zip(0..).collect();
        let expected: HashMap<String, usize> = pairs.iter().cloned().collect();
        let map: SeqIndexMap = pairs.into_iter().collect();
        assert_eq!(map.len(), expected.len());
        for name in names.iter().map(String::as_str).chain(["chr24", "chr", "chrZ", "nope", "Chr1"]) {
            assert_eq!(map.get(name), expected.get(name), "{name:?}");
            assert_eq!(map.contains_key(name), expected.contains_key(name), "{name:?}");
        }
    }

    #[test]
    fn duplicates_count_once() {
        let map: SeqIndexMap =
            [("chr1", 0), ("chr1", 7), ("seq", 1), ("seq", 2)].map(|(n, i)| (n.to_string(), i)).into_iter().collect();
        assert_eq!(map.len(), 2);
        assert_eq!((map.get("chr1"), map.get("seq")), (Some(&7), Some(&2)));
    }
}
