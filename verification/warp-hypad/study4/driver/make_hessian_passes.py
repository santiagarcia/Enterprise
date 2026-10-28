"""Study-4 OTI Hessian pass schedule (pre-registration artifact; no objective evaluation).

Builds a cyclic Steiner triple system STS(109) from a Heffter difference partition of {1..54}
(found by seeded random search), then drops point 108. The result is a list of 1962 seed sets over
the 108 design coefficients (1908 triples, 54 pairs) that contains every unordered coefficient pair
exactly once. One OTI m3n2 pass per seed set yields all second derivatives among its seeds.
Usage: python -I make_hessian_passes.py OUT.txt
"""
import hashlib
import itertools
import random
import sys

V, M, N = 109, 18, 54


def heffter() -> list[tuple[int, int, int]]:
    """Partition {1..54} into 18 triples with b + c = a or a + b + c = 109 (deterministic seeds)."""
    for seed in range(200000):
        rnd = random.Random(seed)
        used: set[int] = set()
        out = []
        for _ in range(M):
            free = [i for i in range(1, N + 1) if i not in used]
            a = max(free)
            opts = [(b, c) for b in free if b != a for c in (a - b, V - a - b)
                    if 1 <= c <= N and c not in (a, b) and c in free and b < c]
            if not opts:
                break
            b, c = rnd.choice(opts)
            used |= {a, b, c}
            out.append((b, c, a))
        else:
            return out
    raise SystemExit("no partition found")


def main() -> None:
    H = heffter()
    T = {tuple(sorted((x + s) % V for x in (0, b, (b + c) % V))) for b, c, _ in H for s in range(V)}
    assert len(T) == M * V
    cnt: dict[tuple[int, int], int] = {}
    for t in T:
        for p in itertools.combinations(t, 2):
            cnt[p] = cnt.get(p, 0) + 1
    assert len(cnt) == V * (V - 1) // 2 and set(cnt.values()) == {1}
    passes = sorted(tuple(x for x in t if x != V - 1) for t in T)
    c2: dict[tuple[int, int], int] = {}
    for t in passes:
        for p in itertools.combinations(t, 2):
            c2[p] = c2.get(p, 0) + 1
    assert len(c2) == 108 * 107 // 2 and set(c2.values()) == {1}
    assert set(itertools.chain.from_iterable(passes)) == set(range(108))
    s = "\n".join(",".join(map(str, t)) for t in passes) + "\n"
    with open(sys.argv[1], "w", newline="\n") as fh:
        fh.write(s)
    print("difference triples:", H)
    print("passes", len(passes), "two-seed", sum(len(t) == 2 for t in passes))
    print("sha256", hashlib.sha256(s.encode()).hexdigest())


if __name__ == "__main__":
    main()
