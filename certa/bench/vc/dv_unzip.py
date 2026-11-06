#!/usr/bin/env python3
"""dv_unzip.py ZIP DEST: unpack a Bazel python zip as `unzip` would: symlinks
restored (zipfile writes them as text files), mode bits kept; members may not
leave DEST. Run with python3 -I."""
import os, stat, sys, zipfile

z = zipfile.ZipFile(sys.argv[1])
dest = os.path.realpath(sys.argv[2])
links = []
for info in z.infolist():
    p = os.path.realpath(os.path.join(dest, info.filename))
    if not (p == dest or p.startswith(dest + os.sep)):
        sys.exit(f"{info.filename} leaves {dest}")
    mode = info.external_attr >> 16
    if stat.S_ISLNK(mode):
        links.append((info, p))
        continue
    z.extract(info, dest)
    if mode and not info.is_dir():
        os.chmod(p, mode & 0o7777)
for info, p in links:
    os.makedirs(os.path.dirname(p), exist_ok=True)
    if os.path.lexists(p):
        os.remove(p)
    os.symlink(z.read(info).decode(), p)
print(f"{len(z.infolist())} members, {len(links)} symlinks")
