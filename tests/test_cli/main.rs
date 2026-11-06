#![cfg(feature = "cli")]

use std::fs;
use std::process::{Command, Output};

const TWOBIT: &str = "tests/assets/tiny-lilend.2bit";
// seq_xy  = gTCCTGCTccaGAAGCAATAACTGATAACNNNNnnnngatcaGCAAGACAATTGAAGAAt (60 nt)
// seq_abc = NNGCTgtgcanNNNNNGACTCCTAcctcnn (30 nt)
const SEQ_ABC: &str = "NNGCTgtgcanNNNNNGACTCCTAcctcnn";

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_twobitreader")).args(args).output().unwrap()
}

// Runs the command on TWOBIT with the given options and returns its output, which must succeed.
fn fasta(options: &[&str]) -> String {
    let out = run(&[options, &[TWOBIT, "stdout"]].concat());
    assert!(out.status.success(), "{options:?}: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

// Whether the command fails, with a message, on TWOBIT with the given options.
fn fails(options: &[&str]) -> bool {
    let out = run(&[options, &[TWOBIT, "stdout"]].concat());
    !out.status.success() && !out.stderr.is_empty()
}

// Writes a file in a scratch folder and returns its path.
fn tmp(name: &str, contents: &str) -> String {
    let dir = std::env::temp_dir().join("twobitreader-cli-tests");
    fs::create_dir_all(&dir).unwrap();
    let path = dir.join(name).to_str().unwrap().to_string();
    fs::write(&path, contents).unwrap();
    path
}

#[test]
fn all_sequences() {
    let masked =
        format!(">seq_xy\ngTCCTGCTccaGAAGCAATAACTGATAACNNNNnnnngatcaGCAAGACA\nATTGAAGAAt\n>seq_abc\n{SEQ_ABC}\n");
    assert_eq!(fasta(&[]), masked);
    let unmasked = fasta(&["-noMask"]);
    assert!(unmasked.lines().filter(|l| !l.starts_with('>')).all(|l| l == l.to_uppercase()));
    assert_eq!(unmasked.to_uppercase(), masked.to_uppercase());
    // Writing to a file gives the same output as stdout.
    let path = tmp("all.fa", "");
    assert!(run(&[TWOBIT, &path]).status.success());
    assert_eq!(fs::read_to_string(path).unwrap(), masked);
}

#[test]
fn seq_ranges() {
    assert_eq!(fasta(&["-seq=seq_xy", "-start=1", "-end=9"]), ">seq_xy:1-9\nTCCTGCTc\n");
    assert_eq!(fasta(&["-seq=seq_xy:1-9"]), ">seq_xy:1-9\nTCCTGCTc\n");
    assert_eq!(fasta(&["-seq=seq_abc", "-end=3"]), ">seq_abc:0-3\nNNG\n"); // -start defaults to 0
    assert_eq!(fasta(&["-seq=seq_abc"]), format!(">seq_abc\n{SEQ_ABC}\n"));
    // Like twoBitToFa, -start/-end give a range only when end > start; otherwise the whole sequence.
    for range in [&["-start=3"][..], &["-start=9", "-end=5"], &["-start=5", "-end=5"]] {
        assert_eq!(fasta(&[&["-seq=seq_abc"], range].concat()), format!(">seq_abc\n{SEQ_ABC}\n"));
    }
    // ...and without -seq they are ignored, with a warning.
    let out = run(&["-start=1", "-end=4", TWOBIT, "stdout"]);
    assert_eq!(String::from_utf8(out.stdout).unwrap(), fasta(&[]));
    assert!(!out.stderr.is_empty());
    // The range can also be part of the input name.
    let out = run(&[&format!("{TWOBIT}:seq_abc:2-6"), "stdout"]);
    assert_eq!(String::from_utf8(out.stdout).unwrap(), ">seq_abc:2-6\nGCTg\n");
}

#[test]
fn seq_list() {
    let list = tmp("list.txt", "# comment\nseq_abc:0-3\n\nseq_xy:1-4\n");
    assert_eq!(fasta(&[&format!("-seqList={list}")]), ">seq_abc:0-3\nNNG\n>seq_xy:1-4\nTCC\n");
}

#[test]
fn bed_formats() {
    // BED3, BED4, and (below) BED12, where blocks 1..4 and 5..8 of seq_xy are joined without the intron.
    let bed = "tests/assets/tiny-bed3-bed4.bed"; // Starts with a track line.
    assert_eq!(fasta(&[&format!("-bed={bed}")]), ">seq_xy:1-4\nTCC\n>foo\nNNG\n");
    assert_eq!(fasta(&[&format!("-bed={bed}"), "-bedPos"]), ">seq_xy:1-4\nTCC\n>seq_abc:0-3\nNNG\n");
    let bed12 = "tests/assets/tiny-bed12.bed";
    assert_eq!(fasta(&[&format!("-bed={bed12}")]), ">tx\nTCCGCT\n");
}

#[test]
fn bed_minus_strand_is_reverse_complemented() {
    // Expected values verified against twoBitToFa.
    let bed = "tests/assets/tiny-strand.bed";
    assert_eq!(fasta(&[&format!("-bed={bed}")]), ">m\nTCtgGGA\n>b\ngcacAGC\n>b\nGCTgtgc\n");
}

#[test]
fn wraps_at_50() {
    let lens: Vec<_> = fasta(&["-seq=seq_xy"]).lines().skip(1).map(str::len).collect();
    assert_eq!(lens, [50, 10]);
}

#[test]
fn failures() {
    assert!(fails(&["-seq=nope"]));
    assert!(fails(&["-seq=seq_abc", "-end=31"])); // Past the end
    assert!(fails(&["-seq=seq_xy:6-2"]));
    assert!(fails(&["-seq=seq_xy:2-99"]));
    assert!(fails(&["-seq=seq_xy:2-6", "-start=1", "-end=4"]));
    assert!(fails(&["-bed=tests/assets/tiny.fasta", "-seq=seq_xy"])); // Conflicting modes
    assert!(fails(&["-bedPos"]));
    assert!(fails(&["-bpt=x.bpt"]));
    assert!(fails(&["-bogus"]));
    assert!(!run(&["tests/assets/malformed-bad-signature.2bit", "stdout"]).status.success());
    assert!(!run(&["missing.2bit", "stdout"]).status.success());
    assert!(!run(&[TWOBIT]).status.success()); // Missing output argument
}

#[test]
fn corrupt_file_exits_cleanly_without_panic_output() {
    let out = run(&["tests/assets/malformed-truncated-dna.2bit", "stdout"]);
    assert_eq!(out.status.code(), Some(1)); // Not 101, which is a panic.
    let err = String::from_utf8(out.stderr).unwrap();
    assert!(err.starts_with("twobitreader: ") && !err.contains("panicked at"), "{err}");
}

#[test]
fn help_succeeds() {
    assert!(run(&["-help"]).stdout.starts_with(b"twobitreader"));
}

// A BED file large enough to be divided into many blocks, mixing BED3, BED4, BED6 and BED12 records.
// Every i gives exactly one record.
fn big_bed() -> String {
    let mut bed = String::new();
    for i in 0..30_000 {
        let (s, e) = (i % 20, 21 + i % 30);
        let strand = if i % 3 == 0 { '-' } else { '+' };
        bed += &match i % 4 {
            0 => format!("seq_xy {s} {e} r{i}\n"),
            1 => format!("seq_abc {} {} r{i} 0 {strand}\n", s % 10, 11 + i % 19),
            2 => format!("seq_xy 0 40 r{i} 0 {strand} 0 40 0 2 5,6, 1,20,\n"),
            _ => format!("#comment{i}\n\nseq_xy {s} {e}\n"),
        };
    }
    bed.replace(' ', "\t")
}

#[test]
fn threads_give_identical_output() {
    let bed = format!("-bed={}", tmp("big.bed", &big_bed()));
    let expected = fasta(&[&bed]);
    assert_eq!(expected.matches('>').count(), 30_000);
    for threads in ["2", "3", "8", "0"] {
        assert_eq!(fasta(&[&bed, &format!("-threads={threads}")]), expected, "threads={threads}");
    }
    let path = tmp("threads.fa", "");
    assert!(run(&[&bed, "-threads=4", TWOBIT, &path]).status.success());
    assert_eq!(fs::read_to_string(path).unwrap(), expected);

    let list: String = (0..20_000)
        .map(|i| format!("seq_xy:{}-{}\nseq_abc:{}-{}\n", i % 20, 21 + i % 30, i % 10, 11 + i % 19))
        .collect();
    let list = format!("-seqList={}", tmp("big-list.txt", &list));
    assert_eq!(fasta(&[&list, "-threads=5"]), fasta(&[&list]));
}

#[test]
fn threads_report_errors() {
    // A bad line in the middle of a large file fails the run, whichever thread meets it.
    let mut text = big_bed();
    let middle = text.len() / 2;
    text.insert_str(middle + text[middle..].find('\n').unwrap() + 1, &"seq_xy 5 999 bad\n".replace(' ', "\t"));
    let out = run(&[&format!("-bed={}", tmp("big-bad.bed", &text)), "-threads=4", TWOBIT, "stdout"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("past the end"));
    // -threads needs a list of regions, and a number.
    assert!(fails(&["-threads=2"]));
    assert!(fails(&["-threads=2", "-seq=seq_xy"]));
    assert!(fails(&["-threads=x", "-bed=tests/assets/tiny.fasta"]));
}

// Differential test against UCSC twoBitToFa, skipped unless TWOBITTOFA names its executable, e.g.
// `TWOBITTOFA=$(pixi run which twoBitToFa) cargo test --features cli` after `pixi add ucsc-twobittofa`.
#[test]
fn matches_twobittofa() {
    let Some(reference) = std::env::var_os("TWOBITTOFA") else {
        eprintln!("TWOBITTOFA not set; skipping differential test");
        return;
    };
    let bed = "-bed=tests/assets/tiny-mixed.bed";
    let list = format!("-seqList={}", tmp("diff-list.txt", "seq_abc:0-3\nseq_xy\nseq_xy:1-40\n"));
    let whole_file = [vec![], vec!["-noMask"]];
    let regions = [
        vec!["-seq=seq_xy"],
        vec!["-seq=seq_xy", "-start=3", "-end=20"],
        vec!["-seq=seq_xy", "-start=3"],
        vec!["-seq=seq_xy", "-end=7"],
        vec!["-seq=seq_abc", "-start=5", "-end=5"],
        vec!["-seq=seq_abc", "-start=9", "-end=5"],
        vec!["-seq=seq_xy:2-6"],
        vec!["-noMask", "-seq=seq_abc", "-start=1", "-end=5"],
        vec![&list],
        vec![&bed],
        vec![&bed, "-bedPos"],
        vec![&bed, "-noMask"],
        // Error cases: both tools must fail.
        vec!["-seq=nope"],
        vec!["-seq=seq_abc", "-end=99"],
        vec!["-seq=seq_xy:6-2"],
    ];
    let run_to_file = |exe: &std::ffi::OsStr, args: &[&str]| {
        let path = tmp("out.fa", "");
        let ok = Command::new(exe).args(args).arg(&path).output().unwrap().status.success();
        (ok, fs::read_to_string(path).unwrap())
    };
    let ours = env!("CARGO_BIN_EXE_twobitreader").as_ref();
    let check = |args: Vec<&str>| {
        let (ref_ok, ref_out) = run_to_file(&reference, &args);
        let (our_ok, our_out) = run_to_file(ours, &args);
        assert_eq!(ref_ok, our_ok, "exit status differs for {args:?}");
        assert!(!ref_ok || ref_out == our_out, "output differs for {args:?}");
    };
    for file in ["tests/assets/tiny-lilend.2bit", "tests/assets/tiny-bigend.2bit", "tests/assets/single-record.2bit"] {
        for options in whole_file.iter().chain(if file.contains("single") { &[][..] } else { &regions[..] }) {
            check([options.as_slice(), &[file]].concat());
        }
    }
    check(vec!["tests/assets/tiny-lilend.2bit:seq_abc:2-6"]); // Range as part of the input name.
}
