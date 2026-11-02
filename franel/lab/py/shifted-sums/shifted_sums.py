import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

TOP = 26
LOW = 16
WORKERS = 8
SHOW = 40
CONTROL_TOP = 16
CHECK = 18
PAGE_T = [0, -1, -3, -4, 0, -6, -5, -27, -21, -25, -51, -118, -168, -178, -257, -169]
PAGE_ENDPOINT = [
    0,
    -1,
    -1,
    -1,
    3,
    -6,
    0,
    -22,
    6,
    -4,
    -27,
    -67,
    -49,
    -10,
    -78,
    88,
    82,
    209,
    543,
    858,
    335,
    -407,
    3030,
]
LOWS = f"L=vector(2^{LOW},r,6*fromdigits(binary(r-1),3));"

# PARI


def gp(script):
    out = subprocess.run(
        ["gp", "-q"], input=script + "\n", capture_output=True, text=True, check=True
    ).stdout
    return [[int(x) for x in line.split()] for line in out.split("\n") if line.strip()]


def fan(scripts):
    with ThreadPoolExecutor(len(scripts)) as pool:
        return [row for part in pool.map(gp, scripts) for row in part]


def head(j):
    return f"h=6*fromdigits(binary({j}),3)*3^{LOW}+1;"


def sums_script(js):
    body = [LOWS]
    for j in js:
        body.append(
            f"{head(j)} p=0; mn=0; mx=0; "
            f"for(r=1,2^{LOW}, p+=moebius(h+L[r]); if(p<mn,mn=p); if(p>mx,mx=p)); "
            f'print({j}," ",p," ",mn," ",mx);'
        )
    return "\n".join(body)


def cross_script(jobs):
    body = [LOWS]
    for j, start, ls in jobs:
        body.append(
            f"{head(j)} p={start}; ls={ls}; "
            f'for(r=1,2^{LOW}, p+=moebius(h+L[r]); if(p!=0&&sign(p)!=ls, if(ls!=0, print({j}," ",r," ",sign(p)," ",h+L[r])); ls=sign(p)));'
            f'print({j}," ",0," ",ls," ",p);'
        )
    return "\n".join(body)


def low_levels():
    return gp(
        LOWS
        + f' p=0; s=0; for(r=1,2^{LOW}, p+=moebius(1+L[r]); if(r==2^s, print(s," ",p); s++));'
    )


def sequential(top):
    return gp(
        f'p=0; ls=0; k=0; at=0; for(i=0,2^{top}-1, p+=moebius(6*fromdigits(binary(i),3)+1); if(p!=0&&sign(p)!=ls, if(ls!=0, k++; at=i+1); ls=sign(p))); print(p," ",k," ",at);'
    )[0]


def spread(items):
    return [items[k::WORKERS] for k in range(WORKERS) if items[k::WORKERS]]


# THE WALK


def chunk_sums():
    got = fan([sums_script(js) for js in spread(list(range(2 ** (TOP - LOW))))])
    return {j: (p, mn, mx) for j, p, mn, mx in got}


def candidates(sums):
    todo, starts, ls, start = [], {}, 0, 0
    for j in range(len(sums)):
        p, mn, mx = sums[j]
        starts[j] = (start, ls)
        if ls == 0 or (ls > 0 and start + mn < 0) or (ls < 0 and start + mx > 0):
            todo.append(j)
            end = start + p
            if end != 0:
                ls = 1 if end > 0 else -1
            else:
                ls = gp(cross_script([(j, start, ls)]))[-1][2]
        start += p
    return todo, starts


def crossings(sums):
    todo, starts = candidates(sums)
    got = fan([cross_script([(j, *starts[j]) for j in js]) for js in spread(todo)])
    flips = sorted((j * 2**LOW + r, sg, n) for j, r, sg, n in got if r > 0)
    return flips, len(todo)


def positive_share(flips, n):
    pos, sign, at = 0, 1, 1
    for m, sg, _ in flips:
        if m > n:
            break
        if sign > 0:
            pos += m - at
        sign, at = sg, m
    if sign > 0:
        pos += n + 1 - at
    return pos / n


def verb_walk():
    t0 = time.time()
    low = dict(low_levels())
    sums = chunk_sums()
    t1 = time.time()
    flips, todo = crossings(sums)
    t2 = time.time()
    level = dict(low)
    run = 0
    for j in range(len(sums)):
        run += sums[j][0]
        if (j + 1) & j == 0:
            level[LOW + (j + 1).bit_length() - 1] = run
    assert [level[s] for s in range(1, 17)] == PAGE_T
    mu = dict(
        gp("\n".join(f'print({t}," ",moebius(3^{t}+1));' for t in range(2, TOP + 2)))
    )
    ends = [level[t - 1] - level[t - 2] + mu[t] for t in range(2, 25)]
    assert ends == PAGE_ENDPOINT
    early = [m for m, _, _ in flips if m <= 2**CHECK]
    assert sequential(CHECK) == [level[CHECK], len(early), early[-1]]
    print(
        f"T_s = sum_(w in B_s) mu(6w + 1), B_s the 2^s base 3 {{0,1}}-strings of length s, walked in ascending w to s = {TOP}"
    )
    print(
        f"checks: T_1..T_16 equal the page's list; T_(t-1) - T_(t-2) + mu(3^t + 1) equals the page's M_F(x_t; R_t) at t = 2..24; a sequential walk to 2^{CHECK} repeats T_{CHECK}, the flip count and the last flip"
    )
    print(
        "  s |        2^s |       top element |       T_s | T_s/2^(s/2) | flips n <= 2^s | share of n <= 2^s with last nonzero sign +"
    )
    for s in range(TOP + 1):
        n = 2**s
        k = sum(1 for m, _, _ in flips if m <= n)
        print(
            f"{s:3d} | {n:10d} | {3 ** (s + 1) - 2:17d} | {level[s]:9d} | {level[s] / 2 ** (s / 2):11.4f} | {k:14d} | {positive_share(flips, n):.4f}"
        )
    neg = [s for s in range(1, TOP + 1) if level[s] < 0]
    pos = [s for s in range(1, TOP + 1) if level[s] > 0]
    print(f"T_s < 0 at s = {neg}")
    print(f"T_s > 0 at s = {pos}")
    print(f"T_s = 0 at s = {[s for s in range(1, TOP + 1) if level[s] == 0]}")
    print(
        f"cut: n = 2^{TOP} = {2**TOP} terms, largest element 3^{TOP + 1} - 2 = {3 ** (TOP + 1) - 2}, P = {level[TOP]}"
    )
    print(
        f"sign changes of the prefix walk P_n (last nonzero sign flips): {len(flips)} in all, {todo} chunks of 2^{LOW} rewalked"
    )
    print(f"the last {min(SHOW, len(flips))}: n, new sign, element 6w + 1")
    for m, sg, n in flips[-SHOW:]:
        print(f"  n = {m:10d}  {'+' if sg > 0 else '-'}  6w + 1 = {n}")
    print(f"walk {t1 - t0:.1f} s, rewalk {t2 - t1:.1f} s, {WORKERS} gp processes")


# THE CONTROL


def smooth23(x):
    count, a = 0, 1
    while a <= x:
        b = a
        while b <= x:
            count, b = count + 1, 3 * b
        a *= 2
    return count


def verb_control():
    t0 = time.time()
    got = gp(
        "p=0; ls=0; ch=0; ps=0; c=0; s=0; nx=1; "
        f"forstep(n=1,3^{CONTROL_TOP}-2,6, p+=moebius(n); c++; if(p!=0&&sign(p)!=ls, if(ls!=0,ch++); ls=sign(p)); if(ls>0,ps++); "
        'if(n==nx, print(s," ",p," ",ch," ",ps," ",c); s++; nx=3^(s+1)-2));'
    )
    print(
        "control: the whole class, M(x; 6, 1) = sum_(n <= x, n = 1 mod 6) mu(n), at x = 3^(s+1) - 2"
    )
    print(
        "drift = 3/4 - Psi_(2,3)(x), the s = 0 residue of M(x; 6, 1) with M's constant -2 and 1/(L(0, chi_-3)(1 + 2^0)) = 3/2"
    )
    print(
        "  s |                x | M(x; 6, 1) |   drift | M/sqrt(x/6) | flips to x | share with last nonzero sign +"
    )
    for s, p, ch, ps, c in got:
        x = 3 ** (s + 1) - 2
        print(
            f"{s:3d} | {x:16d} | {p:10d} | {0.75 - smooth23(x):7.2f} | {p / (x / 6) ** 0.5:11.4f} | {ch:10d} | {ps / c:.4f}"
        )
    print(f"control {time.time() - t0:.1f} s")


def main():
    verbs = {"walk": verb_walk, "control": verb_control}
    for v in sys.argv[1:] or list(verbs):
        verbs[v]()


if __name__ == "__main__":
    main()
