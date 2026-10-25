"""Scores a fault-injection run of certa check (Aim 4).

  aim4_score.py faults.tsv failures.txt

faults.tsv from aim4_tamper.py; failures.txt from certa check --failures
(verdict, record number, record). Prints, per fault, how many injected
records were rejected and with which verdicts, and the rejections of
untouched records (false alarms; there should be none).
"""
import collections
import sys


def main():
    faults = {}
    sure = {}
    with open(sys.argv[1]) as f:
        next(f)
        for line in f:
            rec, fault, s = line.rstrip("\n").split("\t")
            faults[int(rec)] = fault
            sure[fault] = s == "1"
    rejected = {}
    with open(sys.argv[2]) as f:
        for line in f:
            verdict, rec = line.split("\t", 2)[:2]
            rejected[int(rec)] = verdict
    injected = collections.Counter(faults.values())
    caught = collections.Counter()
    how = collections.defaultdict(collections.Counter)
    for rec, fault in faults.items():
        if rec in rejected:
            caught[fault] += 1
            how[fault][rejected[rec]] += 1
    print("fault\tfalse_by_construction\tinjected\trejected\trate\tverdicts")
    for fault in sorted(injected, key=lambda x: (not sure[x], x)):
        n, c = injected[fault], caught[fault]
        verdicts = ", ".join("%s %d" % kv for kv in how[fault].most_common())
        print("%s\t%s\t%d\t%d\t%.4f\t%s" % (fault, "yes" if sure[fault] else "no", n, c, c / n, verdicts))
    false_alarms = sum(1 for rec in rejected if rec not in faults)
    print("untouched records rejected (false alarms): %d" % false_alarms)


if __name__ == "__main__":
    main()
