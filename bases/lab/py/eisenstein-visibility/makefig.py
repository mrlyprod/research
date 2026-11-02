import math
from pathlib import Path

import numpy as np
from PIL import Image

WIDTH = 880
HEIGHT = 760
SCALE = 36.0
N = 30

SURF = (255, 255, 255)
BLUE = (0, 140, 255)
GRAY = (198, 198, 201)
INK = (0, 0, 0)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
DATA_DIR = ROOT / "data" / HERE.relative_to(ROOT)
OUT = DATA_DIR / "bases-fig.png"


def points():
    cx = WIDTH / 2.0
    cy = HEIGHT / 2.0
    root3 = math.sqrt(3.0)
    out = []
    for a in range(-N, N + 1):
        for b in range(-N, N + 1):
            x = cx + SCALE * (a - b / 2.0)
            y = cy - SCALE * (b * root3 / 2.0)
            if -8.0 <= x <= WIDTH + 8.0 and -8.0 <= y <= HEIGHT + 8.0:
                out.append((a, b, x, y))
    return out


def disc(img, cx, cy, r, colour):
    x0 = max(int(cx - r - 2.0), 0)
    x1 = min(int(cx + r + 2.0) + 1, WIDTH)
    y0 = max(int(cy - r - 2.0), 0)
    y1 = min(int(cy + r + 2.0) + 1, HEIGHT)
    if x0 >= x1 or y0 >= y1:
        return
    xs = np.arange(x0, x1, dtype=np.float64) + 0.5 - cx
    ys = np.arange(y0, y1, dtype=np.float64) + 0.5 - cy
    d = np.hypot(xs[None, :], ys[:, None])
    alpha = np.clip(r + 0.5 - d, 0.0, 1.0)
    mask = alpha > 0.0
    if not mask.any():
        return
    patch = img[y0:y1, x0:x1]
    tint = np.array(colour, dtype=np.float64)
    blended = np.rint(patch + (tint - patch) * alpha[:, :, None])
    img[y0:y1, x0:x1] = np.where(mask[:, :, None], blended, patch)


def render():
    pts = points()
    img = np.empty((HEIGHT, WIDTH, 3), dtype=np.float64)
    img[:, :] = SURF
    for a, b, x, y in pts:
        if (a, b) != (0, 0) and math.gcd(a, b) != 1:
            disc(img, x, y, 2.6, GRAY)
    for a, b, x, y in pts:
        if math.gcd(a, b) == 1:
            disc(img, x, y, 5.0, BLUE)
    cx = WIDTH / 2.0
    cy = HEIGHT / 2.0
    disc(img, cx, cy, 6.0, INK)
    disc(img, cx, cy, 3.4, SURF)
    seen = sum(1 for a, b, _, _ in pts if math.gcd(a, b) == 1)
    total = sum(1 for a, b, _, _ in pts if (a, b) != (0, 0))
    return img.astype(np.uint8), seen, total


def main():
    img, seen, total = render()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(img, mode="RGB").save(OUT, optimize=True)
    print(f"wrote {OUT.relative_to(ROOT)}: {WIDTH}x{HEIGHT} RGB")
    print(f"in-frame points visible = {seen}/{total} = {seen / total:.4f}")


main()
