syms x
u=(x-0.5)^3;
e=1;
a=1;
f=-e*diff(diff(u))+a*diff(u);
n=100;
L=-1;
R=1;
X_=L+1/n:(R-L)/n:R;
Y_=double(subs(u,X_));
U0=[Y_(1),Y_(2)];
U=double(Solve(e,a,f,n,U0,L,R));
hold on;
plot(X_,Y_);
plot(X_,U);


function u=Solve(e,a,f,n,u0,L,R)
u(1)=u0(1);
u(2)=u0(2);
h=(R-L)/n;
r=a/e*h;
for i=3:n
    Y_sum=0;
    df=f;
    for m=0:100
        fi_m=subs(df,L+(i-1)/n);
        if(fi_m==0)
            break;
        end
        df=diff(df);
        Y_m=double(Y(m,e,a,r,fi_m));
        Y_sum=Y_sum+Y_m;
    end
    u(i)=(Y_sum-u(i-2)+(exp(-r)+1)*u(i-1))*exp(r);
end
end

function Ym=Y(m,e,a,r,fi_m)
seg=0;
for l=1:m+1
    tmp=r^l/factorial(l)*(exp(-r)+(-1)^l);
    seg=seg+tmp;
end
Ym=e^(m+1)/a^(m+2)*seg*fi_m;
end