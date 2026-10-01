#!/usr/bin/env python3
"""Count the scatter markers drawn in matplotlib PDFs.

A group of points of one marker and colour is drawn as uses of a marker
XObject ('/Mn Do' or '/Pn Do'). A group of ONE point is cheaper to draw inline,
so matplotlib's PDF backend writes it as a plain closed path ('h' then 'B' for
fill and stroke, or 'S' for stroke only) and a count of XObject uses misses it.
This counts both kinds inside the axes clip.

  count_markers.py PDF...
"""
import re, sys, zlib, collections
for f in sys.argv[1:]:
    data = open(f, "rb").read()
    xo, inline = collections.Counter(), 0
    for m in re.finditer(rb"stream\r?\n", data):
        try: s = zlib.decompressobj().decompress(data[m.end():])   # stops at the end of the zlib data
        except Exception: continue
        for k in re.findall(rb"/([MP]\d+) Do", s): xo[k.decode()] += 1
        # inline closed paths drawn after an axes clip ('re W n'), i.e. one-point groups
        for seg in re.split(rb"re W n", s)[1:]:
            seg = seg.split(b"Q", 1)[0]
            if b" Do" in seg: continue
            inline += len(re.findall(rb"\bh\s+(?:B|S)\b", seg))
    print(f"{f.split('/')[-1]}: {sum(xo.values())} XObject markers + {inline} inline = {sum(xo.values()) + inline}  {dict(xo)}")
