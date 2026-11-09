use twobitreader::{TwobitReader, reverse_complement};

// Standard library
use std::collections::HashMap;
use std::error::Error;
use std::io;
use std::iter::zip;
use std::path::Path;

// Modules
mod util;
use util::*;

fn open_tiny() -> TwobitReader {
    let r = TwobitReader::open("tests/assets/tiny-lilend.2bit");
    assert!(r.is_ok());
    r.unwrap()
}

#[test]
fn test_malformed_header_error() {
    // File has bad Twobit signature
    let r = TwobitReader::open("tests/assets/malformed-bad-signature.2bit");
    assert!(r.is_err() && r.err().unwrap().kind() == io::ErrorKind::InvalidData);

    // File has a duplicate sequence name
    let r = TwobitReader::open("tests/assets/malformed-duplicate-name.2bit");
    assert!(r.is_err() && r.err().unwrap().kind() == io::ErrorKind::InvalidData);
}

#[test]
fn test_single_record_offset() {
    // The only record's data starts immediately after the index, so its offset equals the end of the index.
    let tbr = TwobitReader::open("tests/assets/single-record.2bit").unwrap();
    assert_eq!(tbr.get("chr1", 0, 4), "ACGT");
}

#[test]
#[should_panic(expected = "num_blocks too large")]
fn test_truncated_blocks_error() {
    // File block indices were truncated by end of file.
    let tbr = TwobitReader::open_masked("tests/assets/malformed-truncated-blocks.2bit").unwrap();
    for name in tbr.iter_names() {
        let seq_len = tbr.seq_len(name);
        _ = tbr.get(name, 0, seq_len);
    }
}

#[test]
#[should_panic(expected = "Failed to read DNA")]
fn test_truncated_dna_error() {
    // File DNA sequence was truncated by end of file
    let tbr = TwobitReader::open_masked("tests/assets/malformed-truncated-dna.2bit").unwrap();
    for name in tbr.iter_names() {
        let seq_len = tbr.seq_len(name);
        _ = tbr.get(name, 0, seq_len);
    }
}

#[rustfmt::skip]
fn test_tiny<P: AsRef<Path>>(path: P, masked: bool) -> Result<(), Box<dyn Error>> {
    let tbr = if masked { TwobitReader::open_masked(path)} else { TwobitReader::open(path) }?;

    // Reproduce the expected contents of tiny.2bit file here, as strings.
    let expect_names = vec!["seq_xy", "seq_abc"];
    let expect_seqs = if masked {
        vec!["gTCCTGCTccaGAAGCAATAACTGATAACNNNNnnnngatcaGCAAGACAATTGAAGAAt",
             "NNGCTgtgcanNNNNNGACTCCTAcctcnn"]
    } else {
        vec!["GTCCTGCTCCAGAAGCAATAACTGATAACNNNNNNNNGATCAGCAAGACAATTGAAGAAT",
             "NNGCTGTGCANNNNNNGACTCCTACCTCNN"]
    };
    let expect_twobit = zip(expect_names.iter().copied(), expect_seqs.iter().copied()).collect::<HashMap<_, _>>();

    // Test names and iter_names
    assert_eq!(expect_names, tbr.names());
    assert_eq!(expect_names, tbr.iter_names().collect::<Vec<_>>());

    // Test name lookups (&str, String, &String) and sequence len
    for name in expect_names.iter().copied() {
        assert!(tbr.contains_name(name));
    }

    // Test sequence len
    for name in expect_names.iter().copied() {
        assert_eq!(tbr.seq_len(name), expect_twobit[name].len());
    }

    // Test contains_name negatives (&str, String, &String)
    assert!(!tbr.contains_name(""));
    assert!(!tbr.contains_name("no_such_name"));
    assert!(!tbr.contains_name("no_such_name".to_string()));
    assert!(!tbr.contains_name(&"no_such_name".to_string()));

    // Test iteration over &TwobitReader
    assert_eq!(tbr.num_seqs(), expect_twobit.len());

    // Iterator over all possible subsequences
    let windows = expect_twobit.iter().flat_map(|(name, seq)| {
        (0..=seq.len()).flat_map(move |window_size| {
            (0..seq.len() - window_size + 1).map(move |start| {
                (*name, start, start + window_size)
            })
        })
    });

    // Each window's expected sequence
    let expect_batch = windows.clone().map(|(name, start, end)| &expect_twobit[name][start..end]).collect::<Vec<_>>();

    // Test get for every sub-window
    for ((name, start, end), expect_seq) in zip(windows.clone(), expect_batch.iter().copied()) {
        // TwobitReader get methods
        assert_eq!(expect_seq, tbr.get(name, start, end));
        assert_eq!(expect_seq, tbr.get_inclusive(name, start + 1, end));
        assert_eq!(expect_seq, { let mut dst = String::new(); tbr.get_into(name, start, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::new(); tbr.get_inclusive_into(name, start + 1, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::from("abc"); tbr.get_into(name, start, end, &mut dst); dst });
        assert_eq!(expect_seq, { let mut dst = String::from("abc"); tbr.get_inclusive_into(name, start + 1, end, &mut dst); dst });
    }

    // Test prefetch with all valid ranges, both as reference and via move.
    let windows_vec: Vec<_> = windows.clone().collect();
    tbr.prefetch(&windows_vec);     // Not consumed
    tbr.prefetch(windows_vec);      // Consumed
    tbr.prefetch(windows.clone());  // Consumed
    
    // Test concat with zip (from iterator, from collection)
    let expect = if masked { "aACNNNcaGCATTGA" } else { "AACNNNCAGCATTGA" };
    let starts = [0, 10, 20, 30, 40, 50];
    let ends   = [0, 11, 22, 33, 44, 55];
    let ranges_iter = zip(starts, ends);
    let ranges = ranges_iter.clone().collect::<Vec<_>>();
    let ranges_iter_inclusive = ranges_iter.clone().map(|(start, end)| (start + 1, end));
    let ranges_inclusive = ranges_iter_inclusive.clone().collect::<Vec<_>>();
    assert_eq!(expect, tbr.concat("seq_xy", &ranges));  // Not consumed
    assert_eq!(expect, tbr.concat("seq_xy", ranges));   // Consumed
    assert_eq!(expect, tbr.concat_iter("seq_xy", ranges_iter));
    assert_eq!(expect, tbr.concat_inclusive("seq_xy", &ranges_inclusive));
    assert_eq!(expect, tbr.concat_inclusive("seq_xy", ranges_inclusive));
    assert_eq!(expect, tbr.concat_iter_inclusive("seq_xy", ranges_iter_inclusive));

    // Test concat with array of (start, end) pairs (from array, from iterator, from collection)
    let ranges = [(0, 0), (10, 11), (20, 22), (30, 33), (40, 44), (50, 55)];
    assert_eq!(expect, tbr.concat("seq_xy", ranges));
    assert_eq!(expect, tbr.concat("seq_xy", &ranges));
    assert_eq!(expect, tbr.concat_iter("seq_xy", ranges.into_iter()));
    assert_eq!(expect, tbr.concat("seq_xy", Vec::from(ranges)));

    // Test concat with total_len=0, and different name types (&str, String, &String)
    assert_eq!(tbr.concat("seq_xy", zip(starts, starts).collect::<Vec<_>>()), "");
    assert_eq!(tbr.concat_iter("seq_xy", zip(starts, starts)), "");
    assert_eq!(tbr.concat_iter("seq_xy".to_string(), zip(starts, starts)), "");
    assert_eq!(tbr.concat_iter(&"seq_xy".to_string(), zip(starts, starts)), "");

    // Test get_batch with vecs and iterators
    let windows_vec: Vec<_> = windows.clone().collect();
    assert_eq!(expect_batch, tbr.get_batch(&windows_vec).collect::<Vec<_>>());     // Not consumed
    assert_eq!(expect_batch, tbr.get_batch(windows_vec).collect::<Vec<_>>());      // Consumed
    assert_eq!(expect_batch, tbr.get_batch(windows.clone()).collect::<Vec<_>>());  // Consumed

    // Test get_batch_inclusive with vecs and iterators
    let windows_inclusive = windows.clone().map(|(name, start, end)| (name, start + 1, end));
    let windows_inclusive_vec: Vec<_> = windows_inclusive.clone().collect();
    assert_eq!(expect_batch, tbr.get_batch_inclusive(&windows_inclusive_vec).collect::<Vec<_>>());     // Not consumed
    assert_eq!(expect_batch, tbr.get_batch_inclusive(windows_inclusive_vec).collect::<Vec<_>>());      // Consumed
    assert_eq!(expect_batch, tbr.get_batch_inclusive(windows_inclusive.clone()).collect::<Vec<_>>());  // Consumed
    Ok(())
}

// Run test_tiny on the little-endian version of test.2bit
#[test]
fn test_tiny_lilend() -> Result<(), Box<dyn Error>> {
    test_tiny("tests/assets/tiny-lilend.2bit", false)
}

#[test]
fn test_tiny_lilend_masked() -> Result<(), Box<dyn Error>> {
    test_tiny("tests/assets/tiny-lilend.2bit", true)
}

// Run test_tiny on the big-endian version of test.2bit
#[test]
fn test_tiny_bigend() -> Result<(), Box<dyn Error>> {
    test_tiny("tests/assets/tiny-bigend.2bit", false)
}

#[test]
fn test_tiny_bigend_masked() -> Result<(), Box<dyn Error>> {
    test_tiny("tests/assets/tiny-bigend.2bit", true)
}

#[test]
#[should_panic(expected = "invalid start")]
fn test_get_inclusive_invalid_start_panic() {
    open_tiny().get_inclusive("seq_xy", 0, 10);
}

#[test]
#[should_panic(expected = "invalid range")]
fn test_get_invalid_range_panic() {
    open_tiny().get("seq_xy", 11, 10);
}

#[test]
#[should_panic(expected = "invalid range")]
fn test_prefetch_invalid_range_panic() {
    open_tiny().prefetch([("seq_xy", 11, 10)]);
}

#[test]
#[should_panic(expected = "invalid range")]
fn test_prefetch_inclusive_invalid_range_panic() {
    open_tiny().prefetch_inclusive([("seq_xy", 12, 10)]);
}

#[test]
#[should_panic(expected = "invalid start")]
fn test_prefetch_inclusive_invalid_start_panic() {
    open_tiny().prefetch_inclusive([("seq_xy", 0, 10)]);
}

// An empty 1-based range is written (start, start - 1), and must be accepted by
// prefetch_inclusive just as it is by get_inclusive and concat_inclusive.
#[test]
fn test_prefetch_inclusive_empty_range() {
    let tbr = open_tiny();
    tbr.prefetch_inclusive([("seq_xy", 1, 0)]);
    tbr.prefetch_inclusive([("seq_xy", 11, 10)]);
    assert_eq!(tbr.get_inclusive("seq_xy", 11, 10), "");
    assert_eq!(tbr.concat_inclusive("seq_xy", [(11, 10)]), "");
}

#[test]
#[should_panic(expected = "invalid end")]
fn test_get_invalid_end_panic() {
    open_tiny().get("seq_xy", 0, 61);
}

#[test]
#[should_panic(expected = "sequence name not found")]
fn test_seq_len_invalid_name_panic() {
    open_tiny().seq_len("seq_XY");
}

#[test]
#[should_panic(expected = "sequence name not found")]
fn test_get_invalid_name_panic() {
    open_tiny().get("seq_XY", 0, 10);
}

#[test]
#[should_panic(expected = "sequence name not found")]
fn test_concat_invalid_name_panic() {
    open_tiny().concat("seq_XY", [(0, 1), (2, 3)]);
}

#[test]
#[should_panic(expected = "invalid end")]
fn test_concat_invalid_end_panic() {
    open_tiny().concat("seq_xy", [(0, 1), (2, 61)]);
}

#[test]
#[should_panic(expected = "sequence name not found")]
fn test_concat_inclusive_invalid_name_panic() {
    open_tiny().concat_inclusive("seq_XY", [(1, 1), (3, 3)]);
}

#[test]
#[should_panic(expected = "invalid end")]
fn test_concat_inclusive_invalid_end_panic() {
    open_tiny().concat("seq_xy", [(1, 1), (3, 61)]);
}

#[test]
#[ignore] // hg38.2bit is not committed to git; run only if requested by user.
fn test_hg38() -> Result<(), Box<dyn Error>> {
    // Check that all sample intervals match reference sequence
    let r = TwobitReader::open_masked(download_hg38()?.as_path())?;
    let fasta = read_test_fasta("tests/assets/hg38-expected-seqs.fasta")?;
    for (chrom, start, end, seq) in fasta {
        assert_eq!(r.get(&chrom, start, end), seq);
    }
    Ok(())
}

#[test]
fn test_reverse_complement() {
    let check_pair = |x: &str, y: &str| {
        assert_eq!(*x, reverse_complement(y.to_string()));
        assert_eq!(*y, reverse_complement(x.to_string()));
    };

    check_pair("", "");
    check_pair("A", "T");
    check_pair("C", "G");
    check_pair("N", "N");
    check_pair("a", "t");
    check_pair("c", "g");
    check_pair("n", "n");
    check_pair("AC", "GT");
    check_pair("ac", "gt");
    check_pair("ACG", "CGT");
    check_pair("acg", "cgt");
    check_pair("ACGT", "ACGT");
    check_pair("acgt", "acgt");
    check_pair("ACGTN", "NACGT");
    check_pair("acgtn", "nacgt");
    check_pair(
        "ACGTNaaccggttnnAAACCCGGGTTTNNNaaaaccccggggttttnnnn",
        "nnnnaaaaccccggggttttNNNAAACCCGGGTTTnnaaccggttNACGT",
    );
}

#[cfg(debug_assertions)]
#[test]
#[should_panic(expected = "invalid character")]
fn test_reverse_complement_invalid() {
    let _ = reverse_complement("ACGTCXACGT".to_string());
}

#[test]
fn test_try_get_into() {
    let tbr = open_tiny();
    let mut dst = String::new();
    tbr.try_get_into("seq_abc", 2, 6, &mut dst).unwrap();
    assert_eq!(dst, tbr.get("seq_abc", 2, 6));
    tbr.try_get_into("seq_abc", 5, 5, &mut dst).unwrap(); // Empty ranges are valid.
    assert_eq!(dst, "");
    let mut message = |name, start, end| tbr.try_get_into(name, start, end, &mut dst).unwrap_err().to_string();
    assert_eq!(message("nope", 0, 1), "sequence nope not found");
    assert!(message("seq_abc", 6, 2).contains("seq_abc: invalid range (start 6 > end 2)"));
    assert!(message("seq_abc", 0, 31).contains("seq_abc: end (31) is past the end of the sequence (length 30)"));
}

// ---- The lookup of sequence names, which has a fast path for short UCSC-style names ----
//
// tests/assets/make-name-assets.py writes the files. In each, sequence i has i + 1 bases, so a name that
// finds the wrong sequence is detected by the length of what it finds.

// Opens a file, and checks that every name in it is found, and finds its own sequence.
fn open_names(file_name: &str) -> TwobitReader {
    let tbr = TwobitReader::open(format!("tests/assets/{file_name}")).unwrap();
    for (i, name) in tbr.names().into_iter().enumerate() {
        assert!(tbr.contains_name(name), "{name:?}");
        assert_eq!(tbr.seq_len(name), i + 1, "{name:?}");
    }
    tbr
}

#[test]
fn test_names_with_and_without_fast_slots() {
    let tbr = open_names("names-mixed.2bit");
    assert_eq!(tbr.num_seqs(), 61);
    // Names that the fast path finds (including the first and last slot of each table, and the names either
    // side of a change of row in the table of five character names), and names that use the hash map.
    let present = ["chr1", "chr22", "chr33", "chrX", "chrM", "chr2L", "chr2A", "chrIV", "chr:", "chr@", "chrA"];
    let present_at_edges = ["chr10", "chr99", "chr1X", "chr1@", "chrX0", "chrXX"];
    let present_in_hash_map = ["chrIII", "chrXVI", "chr100", "chr01", "chrx", "Chr1", "chr1_KI270706v1_random"];
    for name in present.into_iter().chain(present_at_edges).chain(present_in_hash_map) {
        assert!(tbr.contains_name(name), "{name:?}");
    }
    // Names that have a slot that no sequence uses, and others that are near misses.
    let absent =
        ["chr23", "chr34", "chrY0", "chr3A", "chrB", "chr1Z", "chr0", "chr[", "chr`", "chr1_", "chr2R_", "chrIIII"];
    let absent_in_hash_map = ["", "chr", "Chr2", "chr100x", "chrEBV", "CM010927.1", "NC_003070.8"];
    for name in absent.into_iter().chain(absent_in_hash_map) {
        assert!(!tbr.contains_name(name), "{name:?}");
    }
}

#[test]
fn test_names_that_differ_only_in_a_nul() {
    // "chr1" and "chr1\0" have different lengths, so they must not share a slot.
    let tbr = open_names("names-nul.2bit");
    assert_eq!(tbr.names(), ["chr1", "chr1\0", "chr12", "chr12\0", "chr1\0\0", "chr\0"]);
    assert!(!tbr.contains_name("chr") && !tbr.contains_name("chr\0\0") && !tbr.contains_name("chr13\0"));
}

#[test]
fn test_duplicate_names_are_rejected() {
    // The duplicate is in the table of five character names, and of four character names. (The hash map is
    // tested above, in test_malformed_header_error.)
    for file_name in ["malformed-duplicate-chr12.2bit", "malformed-duplicate-chrX.2bit"] {
        let r = TwobitReader::open(format!("tests/assets/{file_name}"));
        assert!(r.is_err() && r.err().unwrap().kind() == io::ErrorKind::InvalidData, "{file_name}");
    }
}

#[test]
fn test_no_names() {
    let tbr = open_names("no-sequences.2bit");
    assert_eq!(tbr.num_seqs(), 0);
    assert!(!tbr.contains_name("chr1") && !tbr.contains_name("chr12") && !tbr.contains_name("seq"));
}
