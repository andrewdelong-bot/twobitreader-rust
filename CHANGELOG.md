# Changelog

All notable changes to this project are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). The Rust
crate `twobitreader` and the Python package `twobitreader-rs` share a version number.

## [Unreleased]

## [0.3.0] - 2026-11-06

### Added

- `twobitreader` command-line tool, behind the `cli` feature, a drop-in replacement for the core
  features of UCSC `twoBitToFa` that is several times faster on BED extraction (see `BENCH.md`). Its
  `-threads=N` option extracts a `-bed` or `-seqList` file in parallel, with output in input order,
  and its `-prefetch` option asks the operating system to page in the needed data first, which
  greatly speeds up a 2bit file that is not already in the page cache.
- `TwobitReader::try_get_into`, a version of `get_into` that returns an error instead of panicking.
- Python 3.15 wheels, including free-threaded 3.15t.

### Changed

- Dropped support for Python 3.9, which is end-of-life.
- Sequence names chr1-22/X/Y/M are now looked up without hashing, which is ~2.5x faster per lookup.
  Other names fall back to a hash map as before.
- Narrowed the version requirements of dependencies from `"0"` to the ranges that were tested, and
  committed `Cargo.lock`.

## [0.2.3] - 2026-11-03

### Fixed

- `TwobitReader::open` no longer rejects a valid 2bit file that contains exactly one sequence
  record. The record's data starts immediately after the index, so its offset equals the end
  of the index, which the offset check wrongly treated as invalid.
- Missing word and a misplaced sentence in the Benchmarks section of the documentation.

### Changed

- README.md is now generated from the crate documentation with `cargo rdme`, and CI checks
  that it is up to date. The badges moved from the crate documentation into README.md, so
  they no longer appear on docs.rs.

## [0.2.2] - 2026-10-01

### Added

- README-PY.md to document Python usage and to provide a description to PyPI.

## [0.2.1] - 2026-10-01

First release of `twobitreader-rs` on PyPI. There are no library changes; 0.2.0 reached
crates.io before the wheel-building workflow was in place, so it was never published for
Python.

### Added

- Pre-built wheels on PyPI for CPython 3.9 through 3.14, including the free-threaded
  3.14t build, across Linux (glibc and musl), macOS and Windows, each on x86_64 and
  aarch64. A source distribution is published alongside them for any platform or
  interpreter not covered, PyPy included.
- A release workflow that builds and publishes those wheels when a version tag is pushed.
- A statement on Use of AI.

## [0.2.0] - 2026-10-01

`twobitreader` is now available from Python as **`twobitreader-rs`**, built with PyO3 and
published as a wheel. The Python API mirrors the Rust one, with speed close to pure Rust.

### Added

- **Python bindings** (`import twobitreader_rs`), exposing `TwobitReader` and
  `reverse_complement`. Available for CPython 3.9+ on Linux, macOS and Windows.
- **Free-threaded CPython support.** The extension is marked `gil_used = false` and the
  reader is immutable and shareable. A thread pool over batches scales across cores on a
  free-threaded build.
- **Type stubs** (`py.typed` and `__init__.pyi`), so editors and type checkers see full
  signatures for the compiled module.
- **`prefetch_inclusive`**, the 1-based counterpart of `prefetch`, which was missing.

### Changed

- CI moved from CircleCI to GitHub Actions, now also building and testing the Python wheel.

### Security

Hardened against malformed 2bit files; the issues were unbounded allocations triggered by
a crafted file, not memory-safety defects.

### Known issues

- Reading a 2bit file that another process modifies concurrently yields undefined behaviour.

## [0.1.1] - 2026-09-25

### Changed

- Documentation fixes and updated crate keywords.

## [0.1.0] - 2026-09-24

Initial release of the `twobitreader` crate: sequence extraction with `get`, batched
extraction with `get_batch`, exon concatenation with `concat`, reverse complement, and
`prefetch` for cold files. Sequence records are parsed lazily so that opening a large file
stays cheap. No git tag was made for this version.

[Unreleased]: https://github.com/andrewdelong/twobitreader-rust/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/andrewdelong/twobitreader-rust/compare/v0.2.3...v0.3.0
[0.2.3]: https://github.com/andrewdelong/twobitreader-rust/compare/v0.2.2...v0.2.3
[0.2.2]: https://github.com/andrewdelong/twobitreader-rust/compare/v0.2.1...v0.2.2
[0.2.1]: https://github.com/andrewdelong/twobitreader-rust/compare/v0.2.0...v0.2.1
[0.2.0]: https://github.com/andrewdelong/twobitreader-rust/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/andrewdelong/twobitreader-rust/releases/tag/v0.1.1
[0.1.0]: https://crates.io/crates/twobitreader/0.1.0
