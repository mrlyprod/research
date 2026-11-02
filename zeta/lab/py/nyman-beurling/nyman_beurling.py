import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", "..", "..", "..", ".."))
DATA_DIR = os.path.join(REPO, "data", os.path.relpath(HERE, REPO))
PAIRS = [
    (1, 1),
    (1, 2),
    (2, 3),
    (3, 4),
    (5, 7),
    (1, 9),
    (4, 9),
    (13, 27),
    (10, 81),
    (40, 81),
    (6, 10),
]
DESIGNS = [
    (3, (0, 1)),
    (3, (0, 2)),
    (3, (1, 2)),
    (4, (0, 1)),
    (4, (1, 2, 3)),
    (5, (0, 1)),
    (10, tuple(range(9))),
    (10, (0, 1, 2, 3, 4, 5, 6, 7, 9)),
]

ENGINE = r"""
default(parisizemax, 6000000000);
default(realprecision, PREC);
LC = log(2*Pi) - Euler;
CT = Map();
cots(k) = my(c); if(mapisdefined(CT, k, &c), c, c = vector(k - 1, m, cotan(Pi*m/k)); mapput(CT, k, c); c);
VT = Map();
vv(h, k) = {
  my(r = h % k, key, val);
  if(k < 2, return(0));
  key = [r, k];
  if(mapisdefined(VT, key, &val), return(val));
  val = vector(k - 1, m, (m*r) % k) * cots(k)~ / k;
  mapput(VT, key, val);
  val;
}
nu(h, k) = LC/2 * (1/h + 1/k) + (k - h)/(2*h*k) * log(h/k) - Pi/(2*h*k) * (vv(h, k) + vv(k, h));
gram(m, n) = my(g = gcd(m, n)); nu(m/g, n/g)/g;
quad(h, k) = {
  my(pts = vecsort(concat(vector(h+1, j, (j-1)/h), vector(k+1, j, (j-1)/k)), , 8), tot = 0);
  for(i = 1, #pts - 1, tot += intnum(w = pts[i], pts[i+1], frac(k*w)*frac(h*w)*zetahurwitz(2, w)));
  tot/(h*k);
}
grammat(S) = my(n = #S, G = matrix(n, n)); for(i = 1, n, for(j = 1, i, G[i, j] = gram(S[i], S[j]); G[j, i] = G[i, j])); G;
target(S) = vector(#S, i, (log(S[i]) + 1 - Euler)/S[i])~;
cond1(G) = my(Gi = G^(-1), n = #G); vecmax(vector(n, j, vecsum(abs(G[, j])))) * vecmax(vector(n, j, vecsum(abs(Gi[, j]))));
dist(S) = {
  my(G = grammat(S), b = target(S), v = vector(#S, i, 1/S[i])~, x = matsolve(G, b), y = matsolve(G - v*v~, b));
  [1 - b~*x, 1 - b~*y, cond1(G), x, y];
}
design(base, F, N) = select(n -> setminus(Set(digits(n, base)), Set(F)) == [], [1..N]);
resid(S, c) = {
  my(f(u) = sum(i = 1, #S, c[i]*frac(u/S[i])), P = lcm(S), tot = intnum(u = 0, 1, f(u)^2/u^2));
  for(k = 1, P - 1, tot += intnum(u = k, k+1, (1 - f(u))^2/u^2));
  for(k = 0, P - 1, tot += intnum(v = k, k+1, (1 - f(v))^2*zetahurwitz(2, 1 + v/P))/P^2);
  tot;
}
floorq(S, K, tail) = {
  my(E = select(n -> n < K, S), d = #E + 1, L = vector(K - 1, k, log((k+1)/k)));
  my(Q(x) = my(s = x[1], t = if(tail, s^2, 0)); for(k = 1, K - 1, my(B = 1 + sum(i = 1, #E, x[i+1]*floor(k/E[i]))); t += s^2 - 2*s*B*L[k] + B^2/(k*(k+1))); t);
  my(e(i) = vector(d, j, j == i), c0 = Q(vector(d)), g = vector(d, i, (Q(e(i)) - Q(-e(i)))/4)~, H = matrix(d, d));
  for(i = 1, d, H[i, i] = (Q(e(i)) + Q(-e(i)))/2 - c0; for(j = 1, i - 1, H[i, j] = (Q(e(i) + e(j)) - Q(e(i)) - Q(e(j)) + c0)/2; H[j, i] = H[i, j]));
  c0 - g~*matsolve(H, g);
}
qzero(S, N) = my(r = 0); for(n = 1, N, if(issquarefree(n) && !setsearch(Set(S), n), r = n; break)); r;
"""


def gp(lines, prec=60):
    if not shutil.which("gp"):
        raise SystemExit("PARI is needed: gp is not on PATH")
    os.makedirs(DATA_DIR, exist_ok=True)
    path = os.path.join(DATA_DIR, "run.gp")
    with open(path, "w") as fh:
        fh.write(f"PREC = {prec};" + ENGINE + "\n".join(lines) + "\nquit;\n")
    done = subprocess.run(
        ["gp", "-q", path],
        capture_output=True,
        text=True,
        stdin=subprocess.DEVNULL,
        check=False,
    )
    if done.returncode != 0:
        raise RuntimeError(done.stderr)
    return [line for line in done.stdout.splitlines() if not line.startswith("  ***")]


def verb_gram():
    lines = []
    for h, k in PAIRS:
        lines.append(
            f'my(a = gram({h}, {k}), b = quad({h}, {k})); print("{h} {k} ", a, " ", b, " ", abs(a - b))'
        )
    lines.append('print("norm ", gram(1, 1), " ", LC)')
    for S in ("[1, 2, 3]", "[1, 3]", "[1, 2, 3, 4]", "[1, 3, 4]"):
        lines.append(
            f'my(S = {S}, r = dist(S)); print("d2 {S} ", r[1], " ", resid(S, r[4]))'
        )
    lines.append('print("d2 1 ", dist([1])[1], " ", 1 - (1 - Euler)^2/LC)')
    for row in gp(lines, 40):
        print(row)


def verb_table(levels="2-6", extra=7):
    lo, hi = (int(v) for v in levels.split("-"))
    for prec in (60, 90):
        lines = ['C = 2 + Euler - log(4*Pi); print("burnol C ", C)']
        for level in range(lo, hi + 1):
            lines.append(
                f'my(N = 3^{level}, r = dist([1..N])); print("full {level} ", N, " ", N, " ", r[1], " ", r[2], " ", r[3], " ", C/log(N), " ", r[1]*log(N))'
            )
        for level in range(1, max(hi, extra) + 1):
            lines.append(
                f'my(N = 3^{level}, D = design(3, [0, 1], N), r = dist(D)); print("design {level} ", N, " ", #D, " ", r[1], " ", r[2], " ", r[3], " ", C/log(N), " ", r[1]*log(N))'
            )
        start = time.time()
        rows = gp(lines, prec)
        print("precision", prec, f"digits, {time.time() - start:.1f} s")
        for row in rows:
            print(row)


def verb_zeros(height=1000):
    lines = [
        f'my(z = lfunzeros(1, {height}), s = 2*sum(i = 1, #z, 1/(1/4 + z[i]^2)), tail = 2*(log({height}/(2*Pi)) + 1)/(2*Pi*{height})); print("zeros ", #z, " partial ", s, " tail ", tail, " sum ", s + tail, " C ", 2 + Euler - log(4*Pi))'
    ]
    for row in gp(lines, 30):
        print(row)


def verb_floor(kmax=40):
    lines = []
    for base, F in DESIGNS:
        lines.append(
            'my(S = design({0}, {1}, 100000)); print("q0 {0} {1} ", qzero(S, 100000), " ", S[1..12])'.format(
                base, list(F)
            )
        )
    lines.append(
        f'my(S = design(3, [0, 1], 1000)); for(K = 2, {kmax}, print("floor 3 [0,1] ", K, " ", floorq(S, K, 1), " ", floorq(S, K, 0)))'
    )
    lines.append(
        'my(S = [1..1000]); for(K = 2, 12, print("floor full ", K, " ", floorq(S, K, 1), " ", floorq(S, K, 0)))'
    )
    for base, F in DESIGNS[1:]:
        lines.append(
            'my(S = design({0}, {1}, 1000)); foreach([{2}], K, print("floor {0} {1} ", K, " ", floorq(S, K, 1), " ", floorq(S, K, 0)))'.format(
                base, list(F), kmax
            )
        )
    for row in gp(lines, 60):
        print(row)


if __name__ == "__main__":
    verb = sys.argv[1] if len(sys.argv) > 1 else "table"
    start = time.time()
    if verb == "gram":
        verb_gram()
    elif verb == "table":
        verb_table(*(sys.argv[2:3]))
    elif verb == "floor":
        verb_floor(*(int(v) for v in sys.argv[2:3]))
    elif verb == "zeros":
        verb_zeros(*(int(v) for v in sys.argv[2:3]))
    else:
        raise SystemExit("verbs: gram, table, floor, zeros")
    print(f"runtime {time.time() - start:.1f} s")
