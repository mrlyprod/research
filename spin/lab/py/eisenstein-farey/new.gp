default(parisizemax, 1000000000);
L(x) = #select(c->c, Vec(qfrep([2,1;1,2], x, 1)));
a(n) = sumdiv(n, d, moebius(d)*L(3*n^2/d^2));
n = 200;
print(vector(n, k, a(k)));
quit
