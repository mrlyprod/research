default(parisizemax, 1000000000);
ok(i, j) = gcd(binomial(i+j, i), 6) == 1;
K(m) = my(n=3^m, c=0); for(i=0, n-1, for(j=0, n-1, if(ok(i,j), c++))); c;
row(n) = sum(k=0, n, gcd(binomial(n,k), 6) == 1);
print("K(3^m) m=0..7");
print(vector(8, m, K(m-1)));
print("row(n) n=0..299");
print(vector(300, n, row(n-1)));
quit
