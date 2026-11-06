"""Benchmarks for the `twobitreader` command-line tool versus UCSC `twoBitToFa`.

This mirrors the Rust and Python benchmarks (same exons and transcripts workloads, same hot and cold
cache scenarios, same timing harness and output format), but each timing is the wall-clock time of a
whole command-line run: process start-up, opening the 2bit file, parsing the BED file, extraction and
FASTA formatting. Output goes to /dev/null.

Decompression is not timed: the gzipped BED assets are expanded once, beforehand, into plain BED files
under target/bench_cli/, which both tools read.

  exons        BED4 (no strand column), so every exon is extracted on the + strand.
  transcripts  BED12, one line per transcript, whose blocks are the exons. Both tools join the blocks
               and reverse-complement transcripts on the - strand.

Variants, each starting from the one whole BED file and producing one output, in input order:
  basic     one process.
  parallel  twobitreader: one process with -threads=N, which splits the BED file itself.
            twoBitToFa: has no parallel option, so the BED file is split into N equal chunks of lines
            (`split`), `xargs -P` runs one process per chunk, each writing a temporary file, and the
            files are concatenated in order (into /dev/null, standing in for the final .fa) and
            deleted. All of that is timed.

Prerequisites, from the repository root:

    cargo build --release --features cli
    pixi add ucsc-twobittofa     # from bioconda; or point $TWOBITTOFA at the executable
    uv run --no-project benches/bench_cli/bench_cli.py [name filter ...]
"""

import gzip
import itertools
import mmap
import os
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from enum import Enum
from pathlib import Path
from time import perf_counter, sleep

HG38_PATH = Path("tests/assets/hg38.p13.2bit")
HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/p13/hg38.p13.2bit"
COLD_PATH = Path("tests/assets/hg38.p13.cold.2bit")  # Scratch copy used by cold benchmarks.
WORK_DIR = Path("target/bench_cli")  # Decompressed and pre-split BED files (git-ignored).

# Threads (or processes) per tool: the count that was fastest for each on the benchmark machine.
NUM_PROCS = {"cli": 16, "twobittofa": 8}


class Cache(Enum):
    HOT = "hot"  # The 2bit file is already in the page cache.
    COLD = "cold"  # The 2bit file is not in the page cache at all.


class Parallel(Enum):
    BASIC = "basic"
    PARALLEL = "parallel"


def find_tools():
    """Return {label: argv-prefix} for the tools that can be found."""
    cli = Path("target/release/twobitreader")
    if not cli.exists():
        sys.exit("target/release/twobitreader not found; run: cargo build --release --features cli")
    tools = {"cli": [str(cli.resolve())]}
    ref = os.environ.get("TWOBITTOFA") or shutil.which("twoBitToFa") or shutil.which(".pixi/envs/default/bin/twoBitToFa")
    if ref:
        tools["twobittofa"] = [str(Path(ref).resolve())]
    else:
        print("     twoBitToFa not found (set $TWOBITTOFA): benchmarking only the cli")
    return tools


def bench(name, f, num_iter):
    """Invoke f() num_iter times and print the mean/median/min/max of the reported seconds, in ms."""
    # Any command line arguments act as filters, so that a subset can be re-run quickly.
    filters = [arg for arg in sys.argv[1:] if not arg.startswith("-")]
    if filters and not any(filter in name for filter in filters):
        return
    times = []
    for _ in range(num_iter):
        sleep(0.1)  # Let the machine quiet down between runs.
        times.append(int(f() * 1000))
    stats = f"(med={sorted(times)[len(times) // 2]} min={min(times)} max={max(times)}) {times[:5]}"
    print(f"     \x1b[32m{name}\x1b[0m: avg={sum(times) // len(times)} ms {stats}")


def prepare_beds():
    """Expand the BED assets into plain files. Not timed. Returns {workload: bed_path}."""
    def read_lines(path):
        with gzip.open(path, "rt", encoding="ascii") as f:
            return [line.rstrip("\n").split("\t") for line in f]

    exons = [f"{c[0]}\t{c[1]}\t{c[2]}\t{i}" for i, c in enumerate(read_lines("benches/assets/gencode-exons.bed.gz"))]

    # The transcripts asset has one exon per line. Group them by transcript ID into one BED12 line each.
    # BED12 blocks ascend, so exons are sorted; the amount of DNA read and joined is unchanged.
    groups = defaultdict(list)
    for chrom, start, end, name, _score, strand in (r[:6] for r in read_lines("benches/assets/gencode-transcripts.bed.gz")):
        groups[(name, chrom, strand)].append((int(start), int(end)))
    transcripts = []
    for (name, chrom, strand), blocks in groups.items():
        blocks.sort()
        lo, hi = blocks[0][0], max(end for _, end in blocks)
        sizes = "".join(f"{e - s}," for s, e in blocks)
        starts = "".join(f"{s - lo}," for s, _ in blocks)
        transcripts.append(f"{chrom}\t{lo}\t{hi}\t{name}\t0\t{strand}\t{lo}\t{hi}\t0\t{len(blocks)}\t{sizes}\t{starts}")

    shutil.rmtree(WORK_DIR, ignore_errors=True)
    WORK_DIR.mkdir(parents=True)
    beds = {}
    for workload, lines in (("exons", exons), ("transcripts", transcripts)):
        beds[workload] = WORK_DIR / f"{workload}.bed"
        beds[workload].write_text("\n".join(lines) + "\n")
    return beds


def split_command(label, tool, twobit, bed, sink="/dev/null"):
    """Shell command for a tool with no parallel option: split the BED file into one equal chunk of
    lines per process, extract each chunk to its own temporary file in parallel (the processes never
    wait for one another), concatenate the outputs in order into `sink`, and remove them all."""
    q = shlex.quote
    n, folder = NUM_PROCS[label], q(str(WORK_DIR / "split"))
    job = " ".join(map(q, tool + ["-noMask"])) + f' -bed="$1" {q(str(twobit))} "$1.fa"'
    return (
        f"rm -rf {folder} && mkdir {folder}"
        f" && split -l $(( ($(wc -l < {q(str(bed))}) + {n} - 1) / {n} )) -a 3 {q(str(bed))} {folder}/c_"
        f" && ls {folder}/c_??? | xargs -P {n} -I{{}} sh -c {q(job)} _ {{}}"
        f" && cat {folder}/c_???.fa > {q(sink)} && rm -rf {folder}"
    )


def parallel_command(label, tool, twobit, bed, sink="/dev/null"):
    """Argv of one whole-file parallel run: -threads for the cli, split and xargs for twoBitToFa."""
    if label == "cli":
        return tool + ["-noMask", f"-bed={bed}", f"-threads={NUM_PROCS[label]}", str(twobit), sink]
    return ["sh", "-c", split_command(label, tool, twobit, bed, sink)]


def verify(tools, beds, twobit):
    """Check, untimed, that the tools agree, and that each parallel run's output equals a plain run's."""
    for workload, bed in beds.items():
        sample = WORK_DIR / "sample.bed"
        sample.write_text("".join(itertools.islice(open(bed), 2000)))
        outputs = {
            label: subprocess.run(argv + ["-noMask", f"-bed={sample}", str(twobit), "stdout"], capture_output=True, check=True).stdout
            for label, argv in tools.items()
        }
        if len(set(outputs.values())) != 1:
            sys.exit(f"{workload}: the tools produce different output")
        for label, tool in tools.items():
            plain, joined = WORK_DIR / "plain.fa", WORK_DIR / "joined.fa"
            subprocess.run(tool + ["-noMask", f"-bed={bed}", str(twobit), str(plain)], check=True)
            subprocess.run(parallel_command(label, tool, twobit, bed, str(joined)), check=True)
            if plain.read_bytes() != joined.read_bytes():
                sys.exit(f"{label}/{workload}: parallel output differs from a plain run")
    print("     verified: the tools agree, and parallel output equals a plain run's")


def setup(cache):
    """Return the 2bit path to benchmark against, with the page cache prepared. Not timed."""
    if not HG38_PATH.exists():
        subprocess.run(["curl", "--fail", "--location", "-o", str(HG38_PATH), HG38_URL], check=True)
    if cache is Cache.HOT:
        with open(HG38_PATH, "rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as mm:
            mm[:: mmap.PAGESIZE]  # Touch every page.
        return HG38_PATH
    return make_cold_copy(HG38_PATH)


def make_cold_copy(src):
    """Write a copy of the 2bit file that is not in the page cache (macOS only; None elsewhere)."""
    if sys.platform != "darwin":
        return None
    import fcntl

    with open(src, "rb") as fsrc, open(COLD_PATH, "wb") as fdst:
        fcntl.fcntl(fdst.fileno(), 48, 1)  # F_NOCACHE: the bytes go straight to disk.
        while chunk := fsrc.read(8 << 20):
            fdst.write(chunk)
        fdst.flush()
        os.fsync(fdst.fileno())
    return COLD_PATH


def bench_extract(label, tool, workload, cache, parallel, beds):
    """Time one run of the tool over the workload."""
    twobit, bed = setup(cache), beds[workload]
    if parallel is Parallel.BASIC:
        argv = tool + ["-noMask", f"-bed={bed}", str(twobit), "/dev/null"]  # Lowercase masks were not used in the other benchmarks either.
    else:
        argv = parallel_command(label, tool, twobit, bed)
    tic = perf_counter()
    subprocess.run(argv, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    return perf_counter() - tic


def main():
    tools = find_tools()
    beds = prepare_beds()
    twobit = setup(Cache.HOT)
    verify(tools, beds, twobit)

    caches = [(Cache.HOT, 10)]
    if make_cold_copy(twobit):  # Each cold run rewrites an 850MB copy of the genome, so fewer iterations.
        caches.append((Cache.COLD, 5))
    for cache, iterations in caches:
        for workload in beds:
            for parallel in Parallel:
                for label, tool in tools.items():
                    name = f"hg38_{workload}_{cache.value}_{parallel.value}_{label}"
                    bench(name, lambda: bench_extract(label, tool, workload, cache, parallel, beds), iterations)
    COLD_PATH.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
