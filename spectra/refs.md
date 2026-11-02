# REFS

## OEIS SEQUENCES

| ref | title | url |
|---|---|---|
| A299916 | Name is `a(n) = A299914(2n+1)` only, offset 0, from Sahin and Tan's divisibility paper; the six-pointed-star hole count is a comment, not the definition. Holes of the `n`-th size, not tiles: this tree's slice census is `A299916(n+1)` | https://oeis.org/A299916 |

## PAPERS AND BOOKS

| ref | title | url |
|---|---|---|
| Marstrand 1954 | Some Fundamental Geometrical Properties of Plane Sets of Fractional Dimensions, Proc. LMS s3-4(1), 257-302 - almost-every-line slice dimension in the plane, cited on [spectra](notes/spectra.md) | https://doi.org/10.1112/plms/s3-4.1.257 |
| Mattila 1975 | Hausdorff dimension, orthogonal projections and intersections with planes, Ann. Acad. Sci. Fenn. Ser. A I Math. 1, 227-244 - the hyperplane version of the same statement, which is the one [spectra](notes/spectra.md) needs in `R^3` | https://doi.org/10.5186/aasfm.1975.0110 |
| Cox, Sederberg and Chen 1998 | The moving line ideal basis of planar rational curves, Computer Aided Geometric Design 15(8), 803-827, doi:10.1016/S0167-8396(98)00014-4 - the mu-basis degree identity `mu_1 + mu_2 = n - deg(gcd)` that [CLAIMS](README.md) records as the prior art behind `delta_1 + delta_2 = 12R + 5`, so neither that identity nor the rank-2 freeness under it is claimable. Title, journal, volume, issue, pages and year are resolved at the Crossref record; ScienceDirect refuses fetchers, so the text is not read here and the Index Sum Theorem naming on that ledger line is carried unchecked | https://doi.org/10.1016/S0167-8396(98)00014-4 |
| slice-recurrence-order | Menger Diagonal Slices: A Recurrence of Order ceil(D/2); the shelf paper of the order law, the product formula over 3-adic angle towers and the pinning of rho_dim | ../cuts/papers/slice-sign-law.md |
| slice-sign-even-half | The Even Half of the Slice Sign Law; the shelf paper of the even half at bases 3 and 5, the certificate machines, the transient constant ln(R)/4, the tent rank law and the layer-2 window law | ../cuts/papers/slice-sign-law.md |

## PRIOR ART ON THE BASE-3 SLICE

| ref | title | url |
|---|---|---|
| Perez-Duarte, "Slice of Menger" | Flickr: the base-3 centroid diagonal cut, "a very interesting pattern of stars and hexagons"; Abel credits it as first. The companion animated cross-section at `1438621219` now returns 404 | https://www.flickr.com/photos/sbprzd/1432723128/ |
| Hart, "Mathematical Impressions" | Simons Foundation: The Surprising Menger Sponge Slice, the video that popularised the cut | https://www.simonsfoundation.org/2012/12/10/mathematical-impressions-the-surprising-menger-sponge-slice/ |
| Hart, mirrored | the same video at Scientific American; cite the Simons original | https://www.scientificamerican.com/article/mathematical-impressions-the-surprising-menger-sponge-slice/ |
| Cook 2011 | Code to slice open a Menger sponge - working Python from the base-3 digit predicate, point-sampled into a raster; its prose says the normal runs to `(1, 1, 1)` while the listing sets `normal = (1, 1, 0.5)` | https://www.johndcook.com/blog/2011/08/30/slice-a-menger-sponge/ |
| Abel, "Seeing Stars" | "replace each hexagon with 6 hexagons and 6 triangles, and replace each triangle with 1 hexagon and 3 triangles", and `d = log_3((9+sqrt(33))/2) = 1.8184`, with the author's own hedge that it is a computation and not yet a full proof. The substitution is reproved by exhaustion in [README](../README.md) | http://blog.zacharyabel.com/2012/02/seeing-stars/ |
| Hocking, Bridges 2023 | Three-Dimensional Diagonal Cross-Sections of Four-Dimensional Menger Sponges - generalises the cut along dimension `n` and along a hole-iness parameter `k`, with "the base three expansion" fixed throughout | https://archive.bridgesmathart.org/2023/bridges2023-291.pdf |
| Hocking, Bridges 2024 | Menger-Slice Inspired Fractals based on the Pentagon, Dodecahedron, and 120-Cell - frames the Menger slice as a two-tile closed fractal family and writes "The literature uses the term 'directed-graph iterated function system'", which is [spectra](notes/spectra.md)'s own structural claim at base 3. Different solid, same method: it owns the grammar move without touching the base generalisation | https://archive.bridgesmathart.org/2024/bridges2024-297.pdf |
