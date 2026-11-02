dpoly(D) = (1+t^2)^(D-1)*(1+D*t+t^2);
coef(P,i) = if(i<0 || i>poldegree(P), 0, polcoeff(P,i));
eblock(D) = {
  my(P=dpoly(D), w=(D-1)\2+1, s);
  matrix(w,w,r,c, s=D-3*(r-1); coef(P,s+(c-1)) + if(c==1,0,coef(P,s-(c-1))))
};
sig(D) = my(p=charpoly(eblock(D),x)); vector(poldegree(p),i,-polcoeff(p,poldegree(p)-i));
ladder(D,n) = {
  my(B=eblock(D), v=vector(#B,i,i==1), out=[1]);
  for(l=1,n-1, v=B*v~; v=v~; out=concat(out,v[1]));
  out
};
prodladder(D,n) = {
  my(P=subst(dpoly(D),t,'y), H=1, out=[1]);
  for(l=1,n-1, H=H*subst(P,'y,'y^(3^(l-1))); out=concat(out,polcoeff(H,D*(3^l-1)/2)));
  out
};
for(D=4,6, print("dim ",D); print("  even block ",eblock(D)); print("  charpoly ",charpoly(eblock(D),x)); print("  signature ",sig(D)); print("  terms ",ladder(D,8)); print("  product agrees ",prodladder(D,8)==ladder(D,8)));
quit
