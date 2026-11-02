import numpy as np

KNOWN = [
    14.1347,
    21.0220,
    25.0109,
    30.4249,
    32.9351,
    37.5862,
    40.9187,
    43.3271,
]

N_MERTENS = 50000
N_STACK = 200
N_READOUT = 20000
N_SAMPLES = 8192
BAND = (8.0, 55.0)
PEAK_FACTOR = 3.0


def mobius(n):
    mu = np.ones(n + 1, dtype=np.int64)
    composite = np.zeros(n + 1, dtype=bool)
    for p in range(2, n + 1):
        if composite[p]:
            continue
        composite[p::p] = True
        mu[p::p] *= -1
        square = p * p
        if square <= n:
            mu[square::square] = 0
    mu[0] = 0
    return mu


def mertens(mu):
    running = np.cumsum(mu[1:].astype(np.int64))
    return np.concatenate((np.zeros(1, dtype=np.int64), running))


def stack(mu, n):
    nodes = 0
    brightness = np.empty(n, dtype=np.int64)
    for b in range(1, n + 1):
        brightness[b - 1] = mu[b::b][: n // b].sum()
        nodes += int((np.gcd(np.arange(b + 1), b) == 1).sum())
    return nodes, brightness


def readout(m, n):
    breaches = 0
    for x in range(1, n + 1):
        if int(m[x // np.arange(1, x + 1)].sum()) != 1:
            breaches += 1
    return breaches


def normalised(m, n):
    x = np.arange(1, n + 1, dtype=np.float64)
    return x, m[1 : n + 1] / np.sqrt(x)


def spectrum(x, values, samples):
    log_x = np.log(x)
    grid = np.linspace(log_x[0], log_x[-1], samples)
    resampled = np.interp(grid, log_x, values)
    step = (grid[-1] - grid[0]) / (samples - 1)
    amplitude = np.fft.rfft(resampled * np.hanning(samples))
    gamma = np.fft.rfftfreq(samples, d=step) * 2.0 * np.pi
    return gamma, np.abs(amplitude) ** 2, step


def peaks(gamma, power, band, factor):
    inside = (gamma > band[0]) & (gamma < band[1])
    threshold = np.median(power[inside]) * factor
    found = []
    for i in range(1, power.size - 1):
        if not inside[i]:
            continue
        if power[i] > power[i - 1] and power[i] > power[i + 1] and power[i] > threshold:
            found.append(i)
    return found, threshold


def main():
    mu = mobius(N_MERTENS)
    m = mertens(mu)

    nodes, brightness = stack(mu[: N_STACK + 1], N_STACK)
    depth = m[N_STACK // np.arange(1, N_STACK + 1)]
    agreements = int((brightness == depth).sum())
    breaches = readout(m, N_READOUT)

    x, values = normalised(m, N_MERTENS)
    gamma, power, step = spectrum(x, values, N_SAMPLES)
    found, threshold = peaks(gamma, power, BAND, PEAK_FACTOR)

    print(f"domain  N_mertens = {N_MERTENS}  N_stack = {N_STACK}")
    print(f"samples {N_SAMPLES}  log step {step:.8f}  bin width {gamma[1]:.6f}")
    print(
        f"stack   nodes {nodes}  brightness min {brightness.min()} max {brightness.max()}"
    )
    print(
        f"stack   sum mu(kb) equals M(floor(N/b)) at {agreements} of {N_STACK} denominators"
    )
    print(f"readout sum M(floor(x/n)) = 1 through x = {N_READOUT}, breaches {breaches}")
    print(
        f"band    {BAND[0]:g} to {BAND[1]:g}  peaks {len(found)}  threshold {threshold:.4g}"
    )
    print(f"mertens M({N_MERTENS}) = {int(m[-1])}")
    print()
    print("known      detected  error")
    for target in KNOWN:
        best = min(found, key=lambda i: abs(gamma[i] - target))
        print(f"{target:8.4f}   {gamma[best]:7.2f}   {abs(gamma[best] - target):5.2f}")


if __name__ == "__main__":
    main()
