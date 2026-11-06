#!/usr/bin/env python3
"""dv_fetch.py REPO TAG OUT: download a Docker Hub image (registry API, no
daemon) and unpack its layers in order into OUT/rootfs, applying whiteouts.

Every layer blob is checked against its sha256 digest. Paths are resolved
inside the rootfs by hand (an absolute symlink target counts from the rootfs),
so no member, and no symlink a layer creates, can write outside it. Device
nodes are skipped (no root). Run with python3 -I."""
import hashlib, json, os, shutil, stat, sys, tarfile, urllib.request

repo, tag, out = sys.argv[1:4]
root = os.path.realpath(os.path.join(out, "rootfs"))
blobs = os.path.join(out, "blobs")
os.makedirs(root, exist_ok=True)
os.makedirs(blobs, exist_ok=True)


def get(url, headers=None):
    req = urllib.request.Request(url, headers=headers or {})
    return urllib.request.urlopen(req, timeout=120)


tok = json.load(get(f"https://auth.docker.io/token?service=registry.docker.io&scope=repository:{repo}:pull"))["token"]
H = {"Authorization": f"Bearer {tok}"}
man = json.load(get(f"https://registry-1.docker.io/v2/{repo}/manifests/{tag}",
                    dict(H, Accept="application/vnd.docker.distribution.manifest.v2+json")))
json.dump(man, open(os.path.join(out, "manifest.json"), "w"), indent=1)


def fetch(digest):
    path = os.path.join(blobs, digest.replace(":", "_"))
    if not os.path.exists(path):
        h = hashlib.sha256()
        with get(f"https://registry-1.docker.io/v2/{repo}/blobs/{digest}", H) as r, open(path + ".part", "wb") as f:
            while True:
                b = r.read(1 << 22)
                if not b:
                    break
                h.update(b)
                f.write(b)
        if "sha256:" + h.hexdigest() != digest:
            os.remove(path + ".part")
            sys.exit(f"digest mismatch for {digest}")
        os.rename(path + ".part", path)
    return path


cfg = fetch(man["config"]["digest"])
shutil.copy(cfg, os.path.join(out, "config.json"))


def resolve(rel, follow_last=False, depth=0):
    """Host path for rootfs path `rel`, following symlinks inside the rootfs."""
    if depth > 40:
        raise ValueError(f"symlink loop at {rel}")
    parts = [p for p in rel.split("/") if p not in ("", ".")]
    cur = root
    for i, p in enumerate(parts):
        if p == "..":
            cur = os.path.dirname(cur) if cur != root else root
            continue
        nxt = os.path.join(cur, p)
        last = i == len(parts) - 1
        if os.path.islink(nxt) and (not last or follow_last):
            t = os.readlink(nxt)
            base = "" if t.startswith("/") else os.path.relpath(cur, root)
            cur = resolve(os.path.join(base, t) if base != "." else t, True, depth + 1)
        else:
            cur = nxt
        if not (cur == root or cur.startswith(root + os.sep)):
            raise ValueError(f"{rel} escapes the rootfs")
    return cur


def remove(path):
    if os.path.islink(path) or os.path.isfile(path):
        os.remove(path)
    elif os.path.isdir(path):
        shutil.rmtree(path)


skipped = 0
for n, layer in enumerate(man["layers"]):
    path = fetch(layer["digest"])
    with tarfile.open(path, "r:*") as tf:
        for m in tf:
            name = m.name.lstrip("./")
            if not name or name.startswith("/") or ".." in name.split("/"):
                skipped += 1
                continue
            d, b = os.path.split(name)
            if b.startswith(".wh."):
                parent = resolve(d) if d else root
                if b == ".wh..wh..opq":
                    if os.path.isdir(parent):
                        for e in os.listdir(parent):
                            remove(os.path.join(parent, e))
                else:
                    remove(os.path.join(parent, b[4:]))
                continue
            parent = resolve(d) if d else root
            os.makedirs(parent, exist_ok=True)
            dst = os.path.join(parent, b)
            if m.isdir():
                if os.path.islink(dst) or (os.path.exists(dst) and not os.path.isdir(dst)):
                    remove(dst)
                os.makedirs(dst, exist_ok=True)
                os.chmod(dst, (m.mode & 0o777) | 0o700)
            elif m.issym():
                remove(dst) if os.path.lexists(dst) else None
                os.symlink(m.linkname, dst)
            elif m.islnk():
                src = resolve(m.linkname.lstrip("./"))
                remove(dst) if os.path.lexists(dst) else None
                try:
                    os.link(src, dst)
                except OSError:
                    shutil.copy2(src, dst)
            elif m.isreg():
                remove(dst) if os.path.lexists(dst) else None
                with tf.extractfile(m) as s, open(dst, "wb") as f:
                    shutil.copyfileobj(s, f, 1 << 22)
                os.chmod(dst, (m.mode & 0o777) | 0o600)
            else:
                skipped += 1  # devices, FIFOs
    print(f"layer {n + 1}/{len(man['layers'])} {layer['digest'][:19]} {layer['size'] / 1e6:.0f} MB", flush=True)
print(f"rootfs {root}; {skipped} members skipped", flush=True)
