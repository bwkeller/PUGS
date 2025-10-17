#!/usr/bin/env python3
import struct
from sys import argv

if __name__ == "__main__":
    fname = argv[1]
    n = 2048**3
    pad = ((n>>32)&0xFF) + ((n>>16)&0xFF0000) 
    n &= 0xFFFFFFFF
    with open(fname, 'r+b') as f:
        t, n_, ndim, ng_, nd_, ns_, pad_ = struct.unpack("<dIIIIII", f.read(32))
        f.seek(0)
        f.write(struct.pack("<dIIIIII", t, n, ndim, 0, n, 0, pad))
