
syms x
u=exp(x)+2^(-1/1e-2)*(1+x)^(1+1/1e-2);
e=1e-2;
a=1/(1+x);
f=-e*diff(diff(u))+a*diff(u);
n=128;
i=1;
X_=0:1/n:1;
Y_=double(subs(u,X_));
%U0=[Y_(1),Y_(n+1)];
tic
%U=Solve2_(e,a,f,n,U0);
U0=[0,0];
e=1;
[U,p,e]=Solve1(n,e,U0);
toc
plot(X_,U);
norm(e)
% E=U-Y_;
% data(:,1)=X_;
% data(:,2)=U;
% data(:,3)=Y_;
% data(:,4)=E;
while n<=128
    X_=0:1/n:1;
    Y_=double(subs(u,X_));
    U0=[Y_(1),Y_(n+1)];
    U=double(Solve2_(e,a,f,n,U0));
    E(i)=norm(U-Y_,inf);
    n=n*2;
    i=i+1;
end
i=2;
while i<=5
    e(i-1)=log(E(i-1)/E(i))/log(2);
    i=i+1;
end
% plot(X_,E);
% hold on;
% plot(X_,Y_);
% plot(X_,U);

function [u,p,e]=Solve1(n,e,u0)
h=1/n;
x=0:h:1;
y=0:h:1;
a=zeros(1,n+1);
f=zeros(1,n+1);
p=zeros(1,n+1);
for i=1:n+1
    a(i)=aexact(x(i),y(2));
    f(i)=fexact(x(i),y(2))+e*uyyexact(x(i),y(2))-bexact(x(i),y(2))*uyexact(x(i),y(2));
    p(i)=uexact(x(i),y(2));
end
A1=zeros(1,n-1);
A2=zeros(1,n-1);
A3=zeros(1,n-1);
Y=zeros(1,n-1);
for i=2:n
    ri=h*a(i)/e;
    A1(i-1)=A1_1(ri,a(i),e);
    A2(i-1)=A2_1(ri,a(i),e);
    A3(i-1)=A3_1(ri,a(i),e);
    Y(i-1)=Y0_1(ri,a(i),e)*f(i);
end
%t=A1;
Y(1)=Y(1)-A1(1)*u0(1);
A1(1)=[];
Y(n-1)=Y(n-1)-A3(n-1)*u0(2);
A3(n-1)=[];
u=[u0(1),zhuiganfa(A1,A2,A3,Y),u0(2)];
e=u-p;
end

function u=Solve2(e,a,f,n,u0)
u(1)=u0(1);
u(2)=u0(2);
h=1/n;
df=diff(f);
a1=diff(a);
for i=3:n
    ai=subs(a,(i-1)/n);
    a1i=subs(a1,(i-1)/n);
    ri=h*ai/e;
    fi=subs(f,(i-1)/n);
    dfi=subs(df,(i-1)/n);
    u(i)=(Y0_2(ri,ai,a1i,e)*fi+Y1_2(ri,ai,a1i,e)*dfi-A1_2(ri,ai,a1i,e)*u(i-2)-A2_2(ri,ai,a1i,e)*u(i-1))/A3_2(ri,ai,a1i,e);
end
end

function u=Solve2_(e,a,f,n,u0)
h=1/n;
df=diff(f);
a1=diff(a);
for i=1:n-1
    ai=double(subs(a,i/n));
    a1i=double(subs(a1,i/n));
    ri=h*ai/e;
    fi=double(subs(f,i/n));
    dfi=double(subs(df,i/n));
    A1(i)=A1_2(ri,ai,a1i,e);
    A2(i)=A2_2(ri,ai,a1i,e);
    A3(i)=A3_2(ri,ai,a1i,e);
    Y(i)=Y0_2(ri,ai,a1i,e)*fi+Y1_2(ri,ai,a1i,e)*dfi;
end
%t=A1;
Y(1)=Y(1)-A1(1)*u0(1);
A1(1)=[];
Y(n-1)=Y(n-1)-A3(n-1)*u0(2);
A3(n-1)=[];
u=[u0(1),double(zhuiganfa(double(A1),double(A2),double(A3),double(Y))),u0(2)];
end

function u=Solve2__(n,e,u0)
h=1/n;
x=0:h:1;
a=zeros(1,n+1);
ax=zeros(1,n+1);
f=zeros(1,n+1);
fx=zeros(1,n+1);
for i=1:n+1
    a(i)=aexact(x(i));
    ax(i)=axexact(x(i));
    f(i)=fexact(x(i));
    fx(i)=fxexact(x(i));
end

A1=zeros(1,n-1);
A2=zeros(1,n-1);
A3=zeros(1,n-1);
Y=zeros(1,n-1);
for i=2:n
    ri=h*a(i)/e;
    A1(i-1)=A1_2(ri,a(i),ax(i),e);
    A2(i-1)=A2_2(ri,a(i),ax(i),e);
    A3(i-1)=A3_2(ri,a(i),ax(i),e);
    Y(i-1)=Y0_2(ri,a(i),ax(i),e)*f(i)+Y1_2(ri,a(i),ax(i),e)*fx(i);
end
%t=A1;
Y(1)=Y(1)-A1(1)*u0(1);
A1(1)=[];
Y(n-1)=Y(n-1)-A3(n-1)*u0(2);
A3(n-1)=[];
u=[u0(1),zhuiganfa(A1,A2,A3,Y),u0(2)];
end

function [r]=aexact(x,y)
r=exp(x+y);
end

function [r]=bexact(x,y)
r=exp(3*x+y);
end

function [r]=fexact(x,y)
r=2*exp(x)*sin(pi*y)*(exp(2*y) - 1) - exp(x + y)*(sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1) + exp(x)*sin(pi*y)*(exp(2*y) - 1)*(x - 1)) - exp(3*x + y)*(2*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) + pi*cos(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1)) + 4*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) + exp(x)*sin(pi*y)*(exp(2*y) - 1)*(x - 1) - pi^2*sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1) + 4*pi*exp(2*y)*cos(pi*y)*(exp(x) - 1)*(x - 1);
end

function [r]=uexact(x,y)
r=-sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1);
end

function [r]=uyexact(x,y)
r=- 2*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) - pi*cos(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1);
end

function [r]=uyyexact(x,y)
r=pi^2*sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1) - 4*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) - 4*pi*exp(2*y)*cos(pi*y)*(exp(x) - 1)*(x - 1);
end

function res=P_h_1_i(ri,ai,e)
res=e/ai*(1-exp(-ri));
end

function res=P_p_1_i(ri,ai,e)
res=e/ai*(exp(-2*ri)-exp(-ri));
end

function res=Q_h_1_0i(ri,ai,e)
res=-e/ai^2*(1-exp(-ri)-ri*exp(-ri));
end

function res=Q_p_1_0i(ri,ai,e)
res=-e/ai^2*(exp(-2*ri)-exp(-ri)+ri*exp(-ri));
end

function res=P_h_2_i(ri,ai,a1i,e)
res=a1i*e^2/2/ai^3*((ri^2-2*ri+2)-2*exp(-ri));
end

function res=P_p_2_i(ri,ai,a1i,e)
res=a1i*e^2/2/ai^3*(exp(-2*ri)*(ri^2+2*ri+2)-2*exp(-ri));
end

function res=Q_h_2_0i(ri,ai,a1i,e)
res=-a1i*e^2/2/ai^4*((1-exp(-ri)-ri*exp(-ri))*(ri^2-2*ri)+ri^3*exp(-ri));
end

function res=Q_p_2_0i(ri,ai,a1i,e)
res=-a1i*e^2/2/ai^4*((exp(-2*ri)-exp(-ri)+ri*exp(-ri))*(ri^2+2*ri)-ri^3*exp(-ri));
end

function res=Q_h_2_1i(ri,ai,e)
res=-e^2/ai^3*(1-exp(-ri)-ri*exp(-ri)-ri^2/2*exp(-ri));
end

function res=Q_p_2_1i(ri,ai,e)
res=-e^2/ai^3*(exp(-2*ri)-exp(-ri)+ri*exp(-ri)-ri^2/2*exp(-ri));
end

function res=A1_1(ri,ai,e)
res=-exp(-ri)*(P_h_1_i(ri,ai,e));
end

function res=A2_1(ri,ai,e)
res=exp(-ri)*(P_h_1_i(ri,ai,e)-P_p_1_i(ri,ai,e));
end

function res=A3_1(ri,ai,e)
res=exp(-ri)*(P_p_1_i(ri,ai,e));
end

function res=Y0_1(ri,ai,e)
res=P_p_1_i(ri,ai,e)*Q_h_1_0i(ri,ai,e)-P_h_1_i(ri,ai,e)*Q_p_1_0i(ri,ai,e);
end

function res=A1_2(ri,ai,a1i,e)
res=-exp(-ri)*(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e));
end

function res=A2_2(ri,ai,a1i,e)
res=exp(-ri)*(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))-exp(-ri)*(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e));
end

function res=A3_2(ri,ai,a1i,e)
res=exp(-ri)*(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e));
end

function res=Y0_2(ri,ai,a1i,e)
res=(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e))*(Q_h_1_0i(ri,ai,e)+Q_h_2_0i(ri,ai,a1i,e))-(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))*(Q_p_1_0i(ri,ai,e)+Q_p_2_0i(ri,ai,a1i,e));
end

function res=Y1_2(ri,ai,a1i,e)
res=(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e))*Q_h_2_1i(ri,ai,e)-(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))*Q_p_2_1i(ri,ai,e);
end



function [x]=zhuiganfa(a,b,c,d)
r=size(a);
m=r(2);
r=size(b);
n=r(2);
if size(a)~=size(c)|m~=n-1|size(b)~=size(d)
    error('');
end
u(1)=b(1);
for i=2:n
    l(i-1)=a(i-1)/u(i-1);
    u(i)=b(i)-l(i-1)*c(i-1);
end
y(1)=d(1);
for i=2:n
    y(i)=d(i)-l(i-1)*y(i-1);
end
x(n)=y(n)/u(n);
for i=n-1:-1:1
    x(i)=y(i)/u(i);
    x(i)=(y(i)-c(i)*x(i+1))/u(i);
end
end