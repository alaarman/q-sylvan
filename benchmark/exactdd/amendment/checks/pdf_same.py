#!/usr/bin/env python3
"""Same figure? Compare two matplotlib PDFs with the Info dictionary (dates)
removed: every object's bytes, streams decompressed, in object order.
  pdf_same.py A.pdf B.pdf       prints "same" or "DIFF (...)"; exit status 0 or 1"""
import re, sys, zlib
def objs(path):
    b = open(path, "rb").read()
    out = []
    for m in re.finditer(rb"(\d+) 0 obj\s*(.*?)endobj", b, re.S):
        body = m.group(2)
        if b"/CreationDate" in body or b"/Producer" in body: continue      # the Info dict
        s = re.search(rb"stream\r?\n(.*)\r?\nendstream", body, re.S)
        if s and b"/FlateDecode" in body:
            try: body = body[:s.start()] + zlib.decompress(s.group(1))
            except zlib.error: pass
        out.append(re.sub(rb"/Length \d+", b"", body))
    return out
a, b = objs(sys.argv[1]), objs(sys.argv[2])
print("same" if a == b else f"DIFF ({len(a)} vs {len(b)} objects, {sum(x != y for x, y in zip(a, b))} differ)")
sys.exit(0 if a == b else 1)
