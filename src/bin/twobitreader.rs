//! `twobitreader`: a command-line replacement for UCSC `twoBitToFa`.
//!
//! Build with `cargo build --release --features cli`.

use std::fmt::Arguments;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::ops::Range;
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
use std::sync::{Condvar, Mutex};

use twobitreader::{TwobitReader, reverse_complement};

const LINE_WIDTH: usize = 50;

const USAGE: &str = "\
twobitreader - convert all or part of a .2bit file to FASTA (compatible with twoBitToFa)
usage:
   twobitreader [options] input.2bit output.fa
Use 'stdout' as output.fa to write to standard output.
The input can be a 32-bit or 64-bit .2bit file.
options:
   -seq=name       Restrict output to this sequence. May be name:start-end.
   -start=N        Start at this position in the sequence (zero-based). Requires -seq.
   -end=N          End at this position in the sequence (non-inclusive). Requires -seq.
                   Default is the end of the sequence.
   -seqList=file   File of sequence specs, one per line: name or name:start-end,
                   with half-open zero-based coordinates.
   -noMask         Convert sequence to all upper case. Default keeps lowercase mask.
   -bed=file       Extract the regions in a BED file (BED3 to BED12). BED12 blocks are
                   concatenated, excluding introns, and minus-strand
                   records are reverse-complemented.
   -bedPos         With -bed, name each record chrom:start-end instead of the BED name.
   -threads=N      With -bed or -seqList, extract with N threads (0 = one per core). The
                   whole list is read into memory first. Output order is unchanged.
   -prefetch       With -bed or -seqList, first ask the operating system to page in all the
                   sequence data the list needs. Speeds up a 2bit file that is not in the
                   page cache; slightly slows one that is. The whole list is read first.
   -h, -help       Show this message.
Sequences are written 50 bases per line. A record extracted with a range
is named name:start-end.
";

#[derive(Default)]
struct Options {
    seq: Option<String>,
    start: Option<usize>,
    end: Option<usize>,
    seq_list: Option<String>,
    bed: Option<String>,
    bed_pos: bool,
    no_mask: bool,
    threads: Option<usize>, // None = sequential and streaming; Some(0) = one per core.
    prefetch: bool,
}

fn main() -> ExitCode {
    // The library panics on corrupt files that pass validation at open time; report cleanly.
    std::panic::set_hook(Box::new(|info| {
        let msg = info
            .payload()
            .downcast_ref::<String>()
            .map(String::as_str)
            .or_else(|| info.payload().downcast_ref::<&str>().copied())
            .unwrap_or("unexpected failure");
        eprintln!("twobitreader: {msg}");
        std::process::exit(1);
    }));
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) if e.kind() == io::ErrorKind::BrokenPipe => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("twobitreader: {e}");
            ExitCode::FAILURE
        }
    }
}

fn err<T>(msg: impl Into<String>) -> io::Result<T> {
    Err(io::Error::other(msg.into()))
}

fn usage_error<T>(msg: &str) -> io::Result<T> {
    err(format!("{msg}\n\n{USAGE}"))
}

fn parse_args() -> io::Result<Option<(Options, String, String)>> {
    let mut opts = Options::default();
    let mut positional = Vec::new();
    for arg in std::env::args().skip(1) {
        let Some(body) = arg.strip_prefix('-').filter(|b| !b.is_empty()) else {
            positional.push(arg);
            continue;
        };
        let body = body.strip_prefix('-').unwrap_or(body);
        let (key, value) = match body.split_once('=') {
            Some((k, v)) => (k, Some(v.to_string())),
            None => (body, None),
        };
        let need = |v: Option<String>| v.map_or_else(|| usage_error(&format!("option -{key} requires =value")), Ok);
        let num = |v: Option<String>| -> io::Result<usize> {
            let v = need(v)?;
            v.parse().or_else(|_| usage_error(&format!("bad number for -{key}: {v}")))
        };
        match key {
            "h" | "help" => return Ok(None),
            "seq" => opts.seq = Some(need(value)?),
            "start" => opts.start = Some(num(value)?),
            "end" => opts.end = Some(num(value)?),
            "seqList" => opts.seq_list = Some(need(value)?),
            "bed" => opts.bed = Some(need(value)?),
            "bedPos" => opts.bed_pos = true,
            "noMask" => opts.no_mask = true,
            "threads" => opts.threads = Some(num(value)?),
            "prefetch" => opts.prefetch = true,
            "bpt" | "udcDir" => return usage_error(&format!("option -{key} is not supported")),
            _ => return usage_error(&format!("unrecognized option -{key}")),
        }
    }
    let [input, output] =
        <[String; 2]>::try_from(positional).or_else(|_| usage_error("expected two arguments: input.2bit output.fa"))?;
    let modes = [opts.seq.is_some(), opts.seq_list.is_some(), opts.bed.is_some()];
    if modes.iter().filter(|&&m| m).count() > 1 {
        return usage_error("-seq, -seqList and -bed are mutually exclusive");
    }
    if (opts.start.is_some() || opts.end.is_some()) && opts.seq.is_none() {
        eprintln!("twobitreader: warning: -start and -end are ignored without -seq");
        opts.start = None;
        opts.end = None;
    }
    if opts.seq.as_deref().is_some_and(|s| parse_spec(s).1.is_some()) && (opts.start.is_some() || opts.end.is_some()) {
        return usage_error("-seq=name:start-end cannot be combined with -start or -end");
    }
    if opts.bed_pos && opts.bed.is_none() {
        return usage_error("-bedPos requires -bed");
    }
    if opts.threads.is_some() && opts.bed.is_none() && opts.seq_list.is_none() {
        return usage_error("-threads requires -bed or -seqList");
    }
    if opts.prefetch && opts.bed.is_none() && opts.seq_list.is_none() {
        return usage_error("-prefetch requires -bed or -seqList");
    }
    Ok(Some((opts, input, output)))
}

fn run() -> io::Result<()> {
    let Some((opts, input, output)) = parse_args()? else {
        print!("{USAGE}");
        return Ok(());
    };
    let (input, inline_spec) = match input.find(".2bit:") {
        Some(i) if !std::path::Path::new(&input).exists() => {
            (input[..i + 5].to_string(), Some(input[i + 6..].to_string()))
        }
        _ => (input, None),
    };
    if inline_spec.is_some() && (opts.seq.is_some() || opts.seq_list.is_some() || opts.bed.is_some()) {
        return usage_error("input.2bit:name cannot be combined with -seq, -seqList or -bed");
    }
    let tbr = if opts.no_mask { TwobitReader::open(&input) } else { TwobitReader::open_masked(&input) }
        .map_err(|e| io::Error::new(e.kind(), format!("{input}: {e}")))?;

    let dst: Box<dyn Write> = if output == "stdout" {
        Box::new(io::stdout().lock())
    } else {
        Box::new(File::create(&output).map_err(|e| io::Error::new(e.kind(), format!("{output}: {e}")))?)
    };
    let mut out = BufWriter::with_capacity(1 << 20, dst);

    let mut seq = String::new(); // Reused for every record, to avoid reallocating.
    if let Some(path) = &opts.bed {
        if opts.prefetch {
            prefetch_list(&tbr, path, true)?;
        }
        let process = |line: &str, seq: &mut String, out: &mut Vec<u8>| {
            if is_skippable(line, true) {
                return Ok(());
            }
            write_bed_record(&tbr, &input, line, opts.bed_pos, seq, out)
        };
        run_lines(path, opts.threads, &process, &mut out)?;
    } else if let Some(path) = &opts.seq_list {
        if opts.prefetch {
            prefetch_list(&tbr, path, false)?;
        }
        let process = |line: &str, seq: &mut String, out: &mut Vec<u8>| {
            if is_skippable(line, false) {
                return Ok(());
            }
            write_spec(&tbr, &input, line.trim(), None, None, seq, out)
        };
        run_lines(path, opts.threads, &process, &mut out)?;
    } else if let Some(spec) = inline_spec.as_ref().or(opts.seq.as_ref()) {
        write_spec(&tbr, &input, spec, opts.start, opts.end, &mut seq, &mut out)?;
    } else {
        for name in tbr.iter_names() {
            tbr.get_into(name, 0, tbr.seq_len(name), &mut seq);
            write_record(&mut out, format_args!("{name}"), &seq)?;
        }
    }
    out.flush()
}

// Asks the operating system to page in the sequence data that every line of the BED or seqList file will read.
// Lines that are malformed or out of range are left for the extraction to report.
fn prefetch_list(tbr: &TwobitReader, path: &str, bed: bool) -> io::Result<()> {
    let text = std::fs::read_to_string(path).map_err(|e| io::Error::new(e.kind(), format!("{path}: {e}")))?;
    let mut ranges = Vec::new();
    for line in text.lines().filter(|line| !is_skippable(line, bed)) {
        let (name, range) = if bed {
            let mut cols = line.split_ascii_whitespace();
            let (Some(chrom), Some(start), Some(end)) = (cols.next(), cols.next(), cols.next()) else { continue };
            let (Ok(start), Ok(end)) = (start.parse(), end.parse()) else { continue };
            (chrom, Some((start, end)))
        } else {
            parse_spec(line.trim())
        };
        if tbr.contains_name(name) {
            let len = tbr.seq_len(name);
            let (start, end) = range.unwrap_or((0, len));
            if start <= end && end <= len {
                ranges.push((name, start, end));
            }
        }
    }
    tbr.prefetch(&ranges);
    Ok(())
}

// Processes one input line, writing any output. The String is a buffer to reuse between lines.
type LineFn<'a> = dyn Fn(&str, &mut String, &mut Vec<u8>) -> io::Result<()> + Sync + 'a;

// Runs `process` over every line of the file, writing the output to `out` in input order.
// Sequentially, the file is streamed. With more than one thread, it is read into memory and
// divided into blocks that threads claim in turn; the blocks' outputs are then written in order.
fn run_lines(path: &str, threads: Option<usize>, process: &LineFn, out: &mut dyn Write) -> io::Result<()> {
    let threads = match threads {
        Some(0) => std::thread::available_parallelism().map_or(1, |n| n.get()),
        Some(n) => n,
        None => 1,
    };
    if threads == 1 {
        let (mut seq, mut buf) = (String::new(), Vec::new());
        for_each_line(path, |line| {
            process(line, &mut seq, &mut buf)?;
            if buf.len() >= 1 << 20 {
                out.write_all(&buf)?;
                buf.clear();
            }
            Ok(())
        })?;
        return out.write_all(&buf);
    }

    let text = std::fs::read_to_string(path).map_err(|e| io::Error::new(e.kind(), format!("{path}: {e}")))?;
    let blocks = split_blocks(&text, threads);
    let next_block = AtomicUsize::new(0);
    // One slot per block for its output (or error), filled by a worker and drained in order by this thread.
    let slots = Mutex::new((0..blocks.len()).map(|_| None).collect::<Vec<Option<io::Result<Vec<u8>>>>>());
    let slot_filled = Condvar::new();
    std::thread::scope(|scope| {
        for _ in 0..threads {
            scope.spawn(|| {
                let mut seq = String::new();
                loop {
                    let i = next_block.fetch_add(1, Relaxed);
                    let Some(block) = blocks.get(i) else { break };
                    let mut buf = Vec::with_capacity(block.len() * 5); // Output is several times the size of the input.
                    let result = text[block.clone()].lines().try_for_each(|line| process(line, &mut seq, &mut buf));
                    slots.lock().unwrap()[i] = Some(result.map(|()| buf));
                    slot_filled.notify_one();
                }
            });
        }
        for i in 0..blocks.len() {
            let mut guard = slots.lock().unwrap();
            while guard[i].is_none() {
                guard = slot_filled.wait(guard).unwrap();
            }
            let result = guard[i].take().unwrap();
            drop(guard);
            out.write_all(&result?)?;
        }
        Ok(())
    })
}

// Divides the text into blocks of whole lines. Blocks are small enough that every thread gets several,
// which balances uneven work, and large enough that handing them out costs nothing.
fn split_blocks(text: &str, threads: usize) -> Vec<Range<usize>> {
    let target = (text.len() / (4 * threads)).clamp(16 << 10, 4 << 20);
    let mut blocks = Vec::new();
    let mut start = 0;
    while start < text.len() {
        let mut end = (start + target).min(text.len());
        if end < text.len() {
            end = text[end..].find('\n').map_or(text.len(), |i| end + i + 1);
        }
        blocks.push(start..end);
        start = end;
    }
    blocks
}

// Calls `f` with each line of the file, stripped of its line ending. One String is reused for all lines.
fn for_each_line(path: &str, mut f: impl FnMut(&str) -> io::Result<()>) -> io::Result<()> {
    let file = File::open(path).map_err(|e| io::Error::new(e.kind(), format!("{path}: {e}")))?;
    let mut reader = BufReader::with_capacity(1 << 20, file);
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line)? == 0 {
            return Ok(());
        }
        f(line.trim_end_matches(['\n', '\r']))?;
    }
}

fn is_skippable(line: &str, bed: bool) -> bool {
    let t = line.trim();
    t.is_empty() || t.starts_with('#') || (bed && (t.starts_with("track") || t.starts_with("browser")))
}

// Splits "name" or "name:start-end" (zero-based, half-open). A colon that is not followed by
// a valid range is treated as part of the name.
fn parse_spec(spec: &str) -> (&str, Option<(usize, usize)>) {
    if let Some((name, range)) = spec.rsplit_once(':')
        && let Some((s, e)) = range.split_once('-')
        && let (Ok(s), Ok(e)) = (s.replace(',', "").parse(), e.replace(',', "").parse())
    {
        return (name, Some((s, e)));
    }
    (spec, None)
}

// Adds the input file to an error from the library.
fn in_file(e: io::Error, file: &str) -> io::Error {
    io::Error::new(e.kind(), format!("{file}: {e}"))
}

fn check_range(name: &str, start: usize, end: usize, len: usize) -> io::Result<()> {
    if start >= end {
        return err(format!("{name}: start ({start}) must be less than end ({end})"));
    }
    if end > len {
        return err(format!("{name}: end ({end}) is past the end of the sequence (length {len})"));
    }
    Ok(())
}

// Writes a record: ">header", then the sequence wrapped at LINE_WIDTH columns.
fn write_record(out: &mut (impl Write + ?Sized), header: Arguments, seq: &str) -> io::Result<()> {
    writeln!(out, ">{header}")?;
    for line in seq.as_bytes().chunks(LINE_WIDTH) {
        out.write_all(line)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

// Writes the record for `name` or `name:start-end`. `start`/`end` (from -start/-end) apply only if the
// spec has no range. Like twoBitToFa, they select a range only when end > start; otherwise the whole
// sequence is written.
fn write_spec(
    tbr: &TwobitReader,
    file: &str,
    spec: &str,
    start: Option<usize>,
    end: Option<usize>,
    seq: &mut String,
    out: &mut impl Write,
) -> io::Result<()> {
    let (name, range) = parse_spec(spec);
    if !tbr.contains_name(name) {
        return err(format!("sequence {name} not found in {file}"));
    }
    let len = tbr.seq_len(name);
    let (start, end) = (start.unwrap_or(0), end.unwrap_or(0));
    let range = range.or((end > start).then_some((start, end)));
    let (s, e) = range.unwrap_or((0, len));
    check_range(name, s, e, len)?;
    tbr.get_into(name, s, e, seq);
    match range {
        Some(_) => write_record(out, format_args!("{name}:{s}-{e}"), seq),
        None => write_record(out, format_args!("{name}"), seq),
    }
}

// Writes the record for one BED line (BED3 to BED12). BED12 blocks are joined, and minus-strand
// records are reverse-complemented.
fn write_bed_record(
    tbr: &TwobitReader,
    file: &str,
    line: &str,
    bed_pos: bool,
    seq: &mut String,
    out: &mut impl Write,
) -> io::Result<()> {
    let bad = |what: &str| err(format!("bad BED line ({what}): {line}"));
    let mut cols = line.split_ascii_whitespace();
    let mut next = || cols.next();
    let (Some(chrom), Some(start), Some(end)) = (next(), next(), next()) else {
        return bad("need at least 3 columns");
    };
    let (Ok(start), Ok(end)) = (start.parse::<usize>(), end.parse::<usize>()) else {
        return bad("non-numeric start or end");
    };
    let (name, _score, strand) = (next(), next(), next());
    let (_thick_start, _thick_end, _rgb, count, sizes, starts) = (next(), next(), next(), next(), next(), next());

    if start >= end {
        return err(format!("{chrom}: start ({start}) must be less than end ({end})"));
    }

    if let (Some(count), Some(sizes), Some(starts)) = (count, sizes, starts) {
        let list = |s: &str| s.split(',').filter(|x| !x.is_empty()).map(|x| x.parse::<usize>().ok()).collect();
        let (Ok(count), Some::<Vec<_>>(sizes), Some::<Vec<_>>(starts)) =
            (count.parse::<usize>(), list(sizes), list(starts))
        else {
            return bad("blockCount, blockSizes or blockStarts");
        };
        if sizes.len() != count || starts.len() != count {
            return bad("blockCount does not match blockSizes/blockStarts");
        }
        let ranges: Vec<_> = starts.iter().zip(&sizes).map(|(b, size)| (start + b, start + b + size)).collect();
        if ranges.iter().any(|&(s, e)| s == e || e > end) {
            return bad("block outside the feature");
        }
        tbr.try_get_into(chrom, end, end, seq).map_err(|e| in_file(e, file))?; // Checks the name and end; extracts nothing.
        *seq = tbr.concat(chrom, &ranges);
    } else {
        tbr.try_get_into(chrom, start, end, seq).map_err(|e| in_file(e, file))?;
    }
    if strand == Some("-") {
        *seq = reverse_complement(std::mem::take(seq)); // In place, so no reallocation.
    }
    match name {
        Some(name) if !bed_pos => write_record(out, format_args!("{name}"), seq),
        _ => write_record(out, format_args!("{chrom}:{start}-{end}"), seq),
    }
}
