"""Writes the small 2bit files that test the lookup of sequence names (see "names" in tests/test_twobitreader).

Run from the repository root, with any Python 3 and no packages:

    python3 tests/assets/make-name-assets.py

In every file, sequence i (counting from 0, in file order) is i + 1 bases of "A". The tests find a sequence by its
name and check its length, so a name that finds the wrong sequence is detected.
"""

import struct
from pathlib import Path

FOLDER = Path(__file__).parent


def write_2bit(file_name, names):
    """Writes a little-endian, version 0 2bit file with a sequence for each name (bytes), in this order."""
    records = []
    for i in range(len(names)):
        n_bases = i + 1
        packed = b"\xaa" * -(-n_bases // 4)  # "A" is 0b10; the unused bases at the end are ignored.
        # Length, no N blocks, no mask blocks, reserved.
        records.append(struct.pack("<4I", n_bases, 0, 0, 0) + packed)

    header = struct.pack("<4I", 0x1A412743, 0, len(names), 0)  # Signature, version, sequence count, reserved.
    index_size = sum(1 + len(name) + 4 for name in names)
    offset = len(header) + index_size
    index = b""
    for name, record in zip(names, records):
        index += bytes([len(name)]) + name + struct.pack("<I", offset)
        offset += len(record)
    (FOLDER / file_name).write_bytes(header + index + b"".join(records))


def names_mixed():
    """Names that have a fast slot (chr plus one or two characters), names at the edges of the slot ranges, and
    names that look like them but use the hash map. In alphabetical order, as UCSC writes them."""
    in_slots = (
        [f"chr{n}" for n in range(1, 23)]  # Human, mouse, and so on.
        + ["chr33", "chrX", "chrY", "chrM", "chrW", "chrZ"]  # Chicken, sex chromosomes, mitochondrion.
        + ["chr2L", "chr2R", "chr3L", "chr3R", "chr2A", "chr2B"]  # Fly, chimp.
        + ["chrI", "chrII", "chrIV", "chrV", "chrXI"]  # Worm, yeast.
    )
    # The first and last slot of the table of four character names are chr1 and chrZ, and of five character names
    # chr10 and chrXX. The edges are also on either side of a change of row in the table of five character names:
    # chr1X and chr20 would share a slot if the rows were one slot too short.
    edges = ["chr:", "chr@", "chrA", "chr99", "chr1X", "chr1@", "chrX0", "chrXX"]
    in_hash_map = [
        "chrIII", "chrXVI", "chr100",  # Too long.
        "chr01", "chrx", "chr2l", "Chr1",  # Not in the ranges.
        "chr1_KI270706v1_random", "chr12_KI270904v1_alt", "chrUn_GL000195v1",  # Alternate sequences.
        "CM010927.2", "NC_003070.9", "PGSH02012206.1", "seq_abc",  # Accessions, and something else.
    ]
    names = in_slots + edges + in_hash_map
    assert len(names) == len(set(names))
    return [name.encode() for name in sorted(names)]


def main():
    write_2bit("names-mixed.2bit", names_mixed())
    # Names that differ only in NULs have different lengths, so they must have different sequences.
    write_2bit("names-nul.2bit", [b"chr1", b"chr1\0", b"chr12", b"chr12\0", b"chr1\0\0", b"chr\0"])
    # A name repeated in a file is an error, whichever way the name is stored.
    write_2bit("malformed-duplicate-chr12.2bit", [b"chr3", b"chr12", b"chr4", b"chr12"])
    write_2bit("malformed-duplicate-chrX.2bit", [b"chr3", b"chrX", b"chr4", b"chrX"])
    write_2bit("no-sequences.2bit", [])


if __name__ == "__main__":
    main()
