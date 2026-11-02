# slice-ladder-rate

- Settles the decay rate `r_inf` of the slice excess `rho_dim - fill/3` of the base-3 central diagonal slice: the 3-adic cosine product `1/prod_(k>=2) cos(2 pi/3^k)` against the retired Richardson fit, from exact Perron roots to `dim = 120`.
- `P(t) = (1 + t^2)^(dim-1)(1 + dim t + t^2)`, `fill = P(1)`, the carry matrix `M[c', c] = P[c + dim - 3c']` on `abs(c) <= floor((dim-1)/2)` and its reflection-even block `M_even` of size `ceil(dim/2)`, exact integers; `rho_dim` is the Perron root of `M_even`.
- `rho_dim` comes from inverse iteration at the shift `fill/3` in 100-digit mpmath, closed by the exact rational Collatz-Wielandt bounds `min_i (M_even v)_i/v_i <= rho_dim <= max_i (M_even v)_i/v_i` on the positive iterate `v`; every check raises on failure.

## VERBS

- `product`: `prod_(k>=2) cos(2 pi/3^k)` and `r_inf`, 139 cosines at 60 and 120 digits, both truncated to 40 digits and asserted equal; the tail bound `pi^2/(4 9^140)`, the geometric step ratio against `1/9`; `2 r_inf`, `A = 2/(3 prod)`, the residuals from `3/4`, `4/3`, `8/3`; the two candidates against `r_inf` and the digits the fit shares.
- `rate`: the determinants `2, 12, -66, 3120, -14310, 3894912` at `dim = 2..7`; `rho_dim` at every `dim = 2..120` enclosed below `1e-70` relative, numpy's dominant eigenvalue within `1e-8`, `rho_3` against `(9 + sqrt 33)/2`, `rho_120` to 60 digits; per `dim` the excess `delta_dim`, the one-step ratio `|delta_dim/delta_(dim-1)|`, the scaled excess `|delta_dim| r^dim/(dim-1)` under both candidates and the residual of the form `(2(dim-1)/3) prod cos^(dim-1) (dim + 2 cos)/(dim+2)`; Richardson in `1/dim` of the ratio by Neville at orders 1..4 at `dim = 60, 80, 100, 120`, on one parity at orders 1..5, and the least-squares fit of orders 1..3 over `dim = 60..100` and `60..120` with its parity split; the ratio deflated by the form's `dim`-dependence; Richardson of the scaled excess under both candidates; `|delta_dim| r^dim` at `dim = 60, 100` under both.

## RESULTS

- `prod_(k>=2) cos(2 pi/3^k) = 0.7428747134175312701067637570301888840694`, `r_inf = 1.346122006764217277640827069917397630838`, `2 r_inf = 2.69224401352843455`, `A = 0.8974146711761448`.
- The fit `1.34612251727283689` is `5.10509e-7` above `r_inf` and shares seven digits; the order-3 least-squares fit over `dim = 60..120` reads `0.7428745558107063`, `1.58e-7` below the product, the same size and sign.
- The deflated ratio is within `1.06e-17` of the product at `dim = 120`, the form residual `4.42e-18` there, falling by `0.71519` per dimension; the same-parity Neville at order 5 reads `0.7428747134180743`.
- Richardson of the scaled excess at `dim = 120`, order 3, lands `1.95e-8` from `A` under the product and `0.000161` from it under the fit.

## RUN

```
uv run python -I research/cuts/lab/py/slice-ladder-rate/rate.py
uv run python -I research/cuts/lab/py/slice-ladder-rate/rate.py product
```

- Verbs run in the order given, both without arguments; `product` is instant, `rate` about 45 s on one core.

## WITNESSES

- `claims/diagonal-slice-ladder.md`: the rate row, `r_inf`, `2 r_inf`, `A`, the retired fit and its digit count.
- `claims/slice-ladder-rate-and-dead-routes.md`: the rate row, the scaled-excess row and the settling row, verbs `product` and `rate`.
- `shared/constants.md`: the row `slice-ladder-rate`, verb `product`.
