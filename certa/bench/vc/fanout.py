#!/usr/bin/env python3
"""fanout.py URL1 URL2 DIR NAME...

Streams two gzipped FASTQ URLs (mate 1, mate 2) to every consumer NAME through
the FIFOs DIR/NAME.1.fq and DIR/NAME.2.fq, without storing them on disk.

The mates are read in byte lockstep (their records have the same lengths, so
equal byte offsets are nearly equal record counts), and reading stops while
any consumer lags more than HIGH bytes on a mate. A consumer blocked on one
mate has therefore read that mate up to the input frontier and lags at most
one batch on the other, which it has in full: no deadlock while HIGH exceeds
a consumer's batch. Writes are non-blocking, so a slow consumer never holds
up the others by more than HIGH. A consumer that exits early is dropped (and
reported); the exit status is non-zero if any was dropped or a download
failed."""
import bisect, errno, fcntl, os, selectors, subprocess, sys, time

HIGH = 6 << 30          # max bytes a consumer may lag on a mate
LOCK = 256 << 20        # max bytes one mate's input may run ahead of the other
CHUNK = 4 << 20

urls, d, names = sys.argv[1:3], sys.argv[3], sys.argv[4:]
for n in names:
    for m in (1, 2):
        p = f"{d}/{n}.{m}.fq"
        if not os.path.exists(p):
            os.mkfifo(p)

src = [subprocess.Popen(["bash", "-c", f"set -o pipefail; curl -fsS '{u}' | pigz -dc"], stdout=subprocess.PIPE)
       for u in urls]
infd = [p.stdout.fileno() for p in src]
for fd in infd:
    os.set_blocking(fd, False)

# Each FIFO is opened for writing as soon as its reader has opened it (a
# reader blocked in open() counts), inside the main loop, so a consumer may
# open its mates in any order and at any time; -1 = not open yet (or closed).
OPEN_WAIT = 1800
fds = {(n, m): -1 for n in names for m in (1, 2)}
unopened = set(fds)
t0 = time.time()

starts = [[], []]                                 # committed chunks per mate: start offsets
chunks = [[], []]                                 # and their bytes
head = [0, 0]                                     # first chunk still needed
pend = [bytearray(), bytearray()]                 # read but not yet committed
top = [0, 0]                                      # bytes committed per mate
eof = [False, False]
pos = {k: 0 for k in fds}                         # bytes written per (consumer, mate)
live = set(names)
dropped = []
sel = selectors.DefaultSelector()
reading = [False, False]
writing = set()
closed = set()                                    # written in full and closed


def drop(n, why):
    live.discard(n)
    dropped.append((n, why))
    for m in (1, 2):
        k = (n, m)
        if k in writing:
            sel.unregister(fds[k])
            writing.discard(k)
        if fds[k] >= 0:
            os.close(fds[k])
            fds[k] = -1
        unopened.discard(k)
    print(f"fanout: dropped {n}: {why}", file=sys.stderr, flush=True)


def try_open():
    for k in sorted(unopened):
        n, m = k
        try:
            fds[k] = os.open(f"{d}/{n}.{m}.fq", os.O_WRONLY | os.O_NONBLOCK)
        except OSError as e:
            if e.errno != errno.ENXIO:
                raise
            if time.time() - t0 > OPEN_WAIT:
                drop(n, f"did not open mate {m} within {OPEN_WAIT} s")
            continue
        unopened.discard(k)
        try:
            fcntl.fcntl(fds[k], 1031, 1 << 20)  # F_SETPIPE_SZ: 1 MiB where allowed
        except OSError:
            pass


def lag(m):
    return top[m - 1] - min((pos[(n, m)] for n in live), default=top[m - 1])


def update_reads():
    for i in (0, 1):
        want = (not eof[i] and lag(i + 1) < HIGH and live and
                (eof[1 - i] or top[i] + len(pend[i]) <= top[1 - i] + len(pend[1 - i]) + LOCK))
        if want and not reading[i]:
            sel.register(infd[i], selectors.EVENT_READ, ("in", i))
        elif not want and reading[i]:
            sel.unregister(infd[i])
        reading[i] = want


def update_writes():
    for k in fds:
        n, m = k
        if n not in live or k in unopened or k in closed:
            continue
        want = pos[k] < top[m - 1]
        if want and k not in writing:
            sel.register(fds[k], selectors.EVENT_WRITE, ("out", k))
            writing.add(k)
        elif not want and k in writing:
            sel.unregister(fds[k])
            writing.discard(k)
        if not want and eof[m - 1] and pos[k] == top[m - 1]:
            os.close(fds[k])
            fds[k] = -1
            closed.add(k)


def commit(i):
    if pend[i]:
        starts[i].append(top[i])
        chunks[i].append(bytes(pend[i]))
        top[i] += len(pend[i])
        pend[i] = bytearray()


def trim():
    for i in (0, 1):
        lo = min((pos[(n, i + 1)] for n in live), default=top[i])
        h = head[i]
        while h < len(starts[i]) and starts[i][h] + len(chunks[i][h]) <= lo:
            chunks[i][h] = b""
            h += 1
        head[i] = h
        if h > 4096:
            del starts[i][:h], chunks[i][:h]
            head[i] = 0


def chunk_at(i, off):
    j = bisect.bisect_right(starts[i], off, head[i]) - 1
    if j < head[i] or off >= starts[i][j] + len(chunks[i][j]):
        raise RuntimeError(f"fanout: offset {off} of mate {i + 1} not buffered")
    return memoryview(chunks[i][j])[off - starts[i][j]:]


last = time.time()
update_reads()
while True:
    if unopened:
        try_open()
    if not reading[0] and not reading[1] and not writing and not unopened:
        if all(eof) or not live:
            break
    got = [False, False]
    for key, _ in sel.select(timeout=0.2 if unopened else 5):
        kind, x = key.data
        if kind == "in":
            try:
                b = os.read(infd[x], CHUNK)
            except BlockingIOError:
                continue
            got[x] = True
            if not b:
                commit(x)
                eof[x] = True
            else:
                pend[x] += b
        else:
            n, m = x
            if n not in live or fds[x] < 0:
                continue
            try:
                w = os.write(fds[x], chunk_at(m - 1, pos[x]))
                pos[x] += w
            except BlockingIOError:
                pass
            except BrokenPipeError:
                drop(n, f"closed mate {m} at byte {pos[x]:,}")
    for i in (0, 1):
        if pend[i] and (len(pend[i]) >= CHUNK or not got[i]):
            commit(i)
    trim()
    update_writes()
    update_reads()
    if time.time() - last > 300:
        last = time.time()
        print(f"fanout {time.strftime('%H:%M:%S')}: read {top[0] / 2**30:.1f} + {top[1] / 2**30:.1f} GiB; "
              f"lag {lag(1) / 2**30:.2f} / {lag(2) / 2**30:.2f} GiB; live {sorted(live)}", file=sys.stderr, flush=True)

status = [p.wait() for p in src]
print(f"fanout: read {top[0]:,} + {top[1]:,} bytes; inputs exit {status}; dropped {dropped}", file=sys.stderr)
sys.exit(1 if any(status) or dropped else 0)
