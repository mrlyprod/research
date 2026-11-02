f(a,b,c,d,T)=my(s=2*(a+b+c+d),n=0,k);k=s-3*a;if(k<=T,n+=1+f(b,c,d,k,T));k=s-3*b;if(k<=T,n+=1+f(a,c,d,k,T));k=s-3*c;if(k<=T,n+=1+f(a,b,d,k,T));n;
a(n)=my(T=10^n,r=[-1,2,2,3],s=2*vecsum(r),t=0,k,o);for(i=1,4,if(r[i]<=T,t++);k=s-3*r[i];if(k<=T,o=vecextract(r,Str("^",i));t+=1+f(o[1],o[2],o[3],k,T)));t;
print("root excluded ",concat(vector(4,n,Str(a(n)-4,if(n<4,",","")))));
print("every circle ",concat(vector(5,n,Str(a(n-1),if(n<5,",","")))));
