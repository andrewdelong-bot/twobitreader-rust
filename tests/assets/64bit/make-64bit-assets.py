"""Writes the 2bit files that test the 64-bit ("long") format, in this folder, from the tiny files in the folder above.

Run from the repository root, with any Python 3 and no packages:

    python3 tests/assets/64bit/make-64bit-assets.py

A 64-bit file, which UCSC writes with `faToTwoBit -long`, differs from a 32-bit file in only one way: each
sequence's offset in the index is 64 bits, not 32. (The format's header calls the 32-bit format version 0 and the
64-bit format version 1.) The files here have the names of their 32-bit counterparts, such as tiny-lilend.2bit.
"""

import struct
from pathlib import Path

FOLDER = Path(__file__).parent  # tests/assets/64bit
MAIN = FOLDER.parent  # tests/assets, which has the 32-bit files


def read_index(data):
    """Returns (byte order, [(name, offset)], the position after the index) of a 32-bit file."""
    order = "<" if struct.unpack("<I", data[:4])[0] == 0x1A412743 else ">"
    version, count, reserved = struct.unpack(order + "3I", data[4:16])
    assert version == 0 and reserved == 0
    position, entries = 16, []
    for _ in range(count):
        name = data[position + 1 : position + 1 + data[position]]
        position += 1 + len(name)
        entries.append((name, struct.unpack(order + "I", data[position : position + 4])[0]))
        position += 4
    return order, entries, position


def to_64bit(file_name, patch=None):
    """The 64-bit form of a 32-bit file in tests/assets. The index gets 4 bytes longer per sequence, so the offsets
    are shifted. `patch(offsets, file_length)` may change the offsets, to make a malformed file."""
    data = (MAIN / file_name).read_bytes()
    order, entries, index_end = read_index(data)
    shift = 4 * len(entries)
    offsets = [offset + shift for _, offset in entries]
    file_length = len(data) + shift
    if patch:
        patch(offsets, file_length)
    header = data[:4] + struct.pack(order + "3I", 1, len(entries), 0)
    index = b"".join(bytes([len(name)]) + name + struct.pack(order + "Q", o) for (name, _), o in zip(entries, offsets))
    return header + index + data[index_end:]


def main():
    for name in ["tiny-lilend.2bit", "tiny-bigend.2bit"]:
        (FOLDER / name).write_bytes(to_64bit(name))

    # The first sequence's offset is somewhere that no record can be: far too large, past the end of the file,
    # and so near the end that the record's four 32-bit header fields do not fit.
    def set_first(value):
        return lambda offsets, file_length: offsets.__setitem__(0, value(file_length))

    malformed = {
        "malformed-offset-huge.2bit": set_first(lambda length: 2**64 - 16),
        "malformed-offset-past-end.2bit": set_first(lambda length: length + 1000),
        "malformed-offset-near-end.2bit": set_first(lambda length: length - 8),
    }
    for name, patch in malformed.items():
        (FOLDER / name).write_bytes(to_64bit("tiny-lilend.2bit", patch))

    # A count of 4 sequences, which 31 bytes would hold if offsets were 4 bytes, but not with 8 byte offsets:
    # header, one 13 byte index entry for "chr1" (offset 29), and one 18 byte record.
    entry = bytes([4]) + b"chr1" + struct.pack("<Q", 29)
    record = struct.pack("<4I", 8, 0, 0, 0) + bytes([0xAA, 0xAA])
    header = struct.pack("<4I", 0x1A412743, 1, 4, 0)
    (FOLDER / "malformed-count.2bit").write_bytes(header + entry + record)

    # A format that does not exist: neither 32-bit (0) nor 64-bit (1). This one goes in the folder above, since it
    # is a 32-bit file otherwise.
    data = bytearray((MAIN / "tiny-lilend.2bit").read_bytes())
    data[4:8] = struct.pack("<I", 2)
    (MAIN / "malformed-unknown-format.2bit").write_bytes(data)


if __name__ == "__main__":
    main()
