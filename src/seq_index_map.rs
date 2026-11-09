// A specialized stand-in for `HashMap<String, usize>`, mapping sequence names to indices into
// `TwobitReader::seqs`.
//
// Genomes are accessed almost entirely by their primary sequences (the chromosomes and the mitochondrion),
// whose UCSC-style names are short: chr1, chrX, chr2L, chrM. Any name of the form "chr" plus one or two
// characters in the ranges below is found by array lookup, with no hashing or string comparison. Every other
// name (alternate contigs such as chr1_KI270706v1_random, scaffolds, accessions, and the longer chrXII) uses
// a HashMap, so every name is found whichever way it is stored.

use std::collections::HashMap;

// Marks an entry that no sequence maps to.
const EMPTY: u32 = u32::MAX;

// Names of length 4: "chr" and one more byte, from '1' to 'Z' (so digits, upper case letters, and a few punctuation marks).
const NUM_SLOTS_LEN4: usize = (b'Z' - b'1' + 1) as usize; // 42

// Names of length 5: "chr" and two more bytes: the first from '1' to 'X', the second from '0' to 'X'. That covers chr10..chr99,
// the arms such as chr2L and chr3R, and roman numerals up to chrXI, chrXV and chrXX.
const NUM_FIRST_LEN5: usize = (b'X' - b'1' + 1) as usize; // 40
const NUM_SECOND_LEN5: usize = (b'X' - b'0' + 1) as usize; // 41
const NUM_SLOTS_LEN5: usize = NUM_FIRST_LEN5 * NUM_SECOND_LEN5; // 1640

// Which table, and which slot of it, a name belongs to, if it is of length 4 or 5 and in the ranges.
enum Slot {
    Len4(usize),
    Len5(usize),
}

#[inline(always)]
fn slot(name: &[u8]) -> Option<Slot> {
    match *name {
        [b'c', b'h', b'r', a @ b'1'..=b'Z'] => Some(Slot::Len4(usize::from(a - b'1'))),
        [b'c', b'h', b'r', a @ b'1'..=b'X', b @ b'0'..=b'X'] => {
            Some(Slot::Len5(usize::from(a - b'1') * NUM_SECOND_LEN5 + usize::from(b - b'0')))
        }
        _ => None,
    }
}

#[derive(Debug)]
pub(crate) struct SeqIndexMap {
    slots_len4: Vec<u32>,          // The sequence index of each length 4 name, by slot, or EMPTY.
    slots_len5: Vec<u32>,          // The sequence index of each length 5 name, by slot, or EMPTY.
    other: HashMap<String, usize>, // The sequence index of every name that has no slot.
}

impl SeqIndexMap {
    // Returns the sequence index for the name, if there is one.
    // A name with a slot is in its slot or nowhere, so a miss there is final.
    #[inline(always)]
    pub(crate) fn get(&self, name: &str) -> Option<usize> {
        let index = match slot(name.as_bytes()) {
            Some(Slot::Len4(slot)) => self.slots_len4[slot],
            Some(Slot::Len5(slot)) => self.slots_len5[slot],
            None => return self.other.get(name).copied(),
        };
        (index != EMPTY).then_some(index as usize)
    }

    pub(crate) fn contains_key(&self, name: &str) -> bool {
        self.get(name).is_some()
    }

    // Returns the number of distinct names.
    pub(crate) fn len(&self) -> usize {
        let used = |slots: &Vec<u32>| slots.iter().filter(|&&index| index != EMPTY).count();
        used(&self.slots_len4) + used(&self.slots_len5) + self.other.len()
    }
}

// As with HashMap, when a name occurs more than once the last index wins and len() counts it once.
impl FromIterator<(String, usize)> for SeqIndexMap {
    fn from_iter<I: IntoIterator<Item = (String, usize)>>(iter: I) -> Self {
        let mut map = Self {
            slots_len4: vec![EMPTY; NUM_SLOTS_LEN4],
            slots_len5: vec![EMPTY; NUM_SLOTS_LEN5],
            other: HashMap::new(),
        };
        for (name, index) in iter {
            // The number of sequences is a u32 in the file, so an index always fits (and is never EMPTY).
            match slot(name.as_bytes()) {
                Some(Slot::Len4(slot)) => map.slots_len4[slot] = index as u32,
                Some(Slot::Len5(slot)) => map.slots_len5[slot] = index as u32,
                None => {
                    map.other.insert(name, index);
                }
            }
        }
        map
    }
}
