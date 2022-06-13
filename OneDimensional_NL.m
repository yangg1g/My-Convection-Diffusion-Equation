clear
% syms x
% u=exp(x)+2^(-1/1e-2)*(1+x)^(1+1/1e-2);
% e=1e-2;
% a=1/(1+x);
% f=-e*diff(diff(u))+a*diff(u);
% n=8;
% i=1;
% X_=0:1/n:1;
% Y_=double(subs(u,X_));
% U0=[Y_(1),Y_(n+1)];
%U=Solve2_(e,a,f,n,U0);
%U0=[0,0];
%e=1;
% tic
n=128;
h=1/n;
x=0:h:1;
[U,pp,ee,Y1]=Solve1__(n,0.01);
% t1=norm(ee)
% [U,pp,ee,Y2]=Solve2__(n,1);
% t2=norm(ee)
% toc
% norm(Y1-Y2)
% plot(X_,U);
% norm(e)
% E=U-Y_;
% data(:,1)=X_;
% data(:,2)=U;
% data(:,3)=Y_;
% data(:,4)=E;
% plot(x,U);
% hold on
% plot(x,pp);
norm(ee,inf);
n=8;
i=1;
while n<=1024
%     X_=0:1/n:1;
%     Y_=double(subs(u,X_));
%     U0=[Y_(1),Y_(n+1)];
    %U=double(Solve1__(e,a,f,n,U0));
    [U,pp,E(i),tt(i)]=Solve1__(n,0.01);
    n=n*2;
    i=i+1;
end
i=2;
fprintf("%d & %d & & %d \\\\\n",i*4,E(i-1),tt(i-1));
while i<=8
    e(i-1)=log(E(i-1)/E(i))/log(2);
    fprintf("%d & %d & %d & %d \\\\\n",2^(i+2),E(i-1),e(i-1),tt(i-1));
    i=i+1;
end
% plot(X_,E);
% hold on;
% plot(X_,Y_);
% plot(X_,U);
function u=Solve1(e,a,f,n,u0)
u(1)=u0(1);
u(2)=u0(2);
h=1/n;
for i=3:n
    ai=subs(a,(i-1)/n);
    ri=h*ai/e;
    fi=subs(f,(i-1)/n);
    u(i)=(Y0_1(ri,ai,e)*fi-A1_1(ri,ai,e)*u(i-2)-A2_1(ri,ai,e)*u(i-1))/A3_1(ri,ai,e);
end
end

function u=Solve1_(e,a,f,n,u0)
h=1/n;
for i=1:n-1
    ai=subs(a,i/n);
    ri=h*ai/e;
    fi=subs(f,i/n);
    A1(i)=A1_1(ri,ai,e);
    A2(i)=A2_1(ri,ai,e);
    A3(i)=A3_1(ri,ai,e);
    Y(i)=Y0_1(ri,ai,e)*fi;
end
Y(1)=Y(1)-A1(1)*u0(1);
A1(1)=[];
Y(n-1)=Y(n-1)-A3(n-1)*u0(2);
A3(n-1)=[];
u=[u0(1),zhuiganfa(A1,A2,A3,Y),u0(2)];
end

function [u,p,e,j]=Solve1__(n,e)
h=1/n;
x=0:h:1;
a=zeros(1,n+1);
f=zeros(1,n+1);
p=zeros(1,n+1);
u_=ones(1,n+1)*1e-5;
u=zeros(1,n+1);
for i=1:n+1
    a(i)=aexact(x(i));
    f(i)=fexact(x(i));
    p(i)=uexact(x(i));
end
A1=zeros(1,n-1);
A2=zeros(1,n-1);
A3=zeros(1,n-1);
Y=zeros(1,n-1);
for j=1:100
    for i=2:n
        ai=u_(i);
        ri=h*ai/e;
        A1(i-1)=A1_1(ri,ai,e);
        A2(i-1)=A2_1(ri,ai,e);
        A3(i-1)=A3_1(ri,ai,e);
        Y(i-1)=Y0_1(ri,ai,e)*f(i);
    end
    %t=A1;
    Y(1)=Y(1)-A1(1)*p(1);
    A1(1)=[];
    Y(n-1)=Y(n-1)-A3(n-1)*p(n+1);
    A3(n-1)=[];
    u=[p(1),zhuiganfa(A1,A2,A3,Y),p(n+1)];
    if(max(abs(u-u_))<1e-8)
        break
    end
    u_=u;
end

e=max(abs(u-p));
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
    ai=subs(a,i/n);
    a1i=subs(a1,i/n);
    ri=h*ai/e;
    fi=subs(f,i/n);
    dfi=subs(df,i/n);
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

function [u,p,ee,A2]=Solve2__(n,e)
h=1/n;
x=0:h:1;
a=zeros(1,n+1);
ax=zeros(1,n+1);
f=zeros(1,n+1);
fx=zeros(1,n+1);
p=zeros(1,n+1);
px=zeros(1,n+1);
for i=1:n+1
    a(i)=aexact(x(i));
    ax(i)=axexact(x(i));
    f(i)=fexact(x(i));
    fx(i)=fxexact(x(i));
    p(i)=uexact(x(i));
    px(i)=uxexact(x(i));
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
Y(1)=Y(1)-A1(1)*p(1);
A1(1)=[];
Y(n-1)=Y(n-1)-A3(n-1)*p(n+1);
A3(n-1)=[];
u=[p(1),zhuiganfa(A1,A2,A3,Y),p(n+1)];
ee=p-u;

% i=6;
% ri=h*a(i)/e;
% exp(-ri)*(p(i+1)-p(i))
% exp(-ri)*(u(i+1)-u(i))
% P_h_1_i(ri,a(i),e)*px(i)+Q_h_1_0i(ri,a(i),e)*f(i)
% (P_h_1_i(ri,a(i),e)+P_h_2_i(ri,a(i),ax(i),e))*px(i)+(Q_h_1_0i(ri,a(i),e)+Q_h_2_0i(ri,a(i),ax(i),e))*f(i)+Q_h_2_1i(ri,a(i),e)*fx(i)
end

function [r]=fexact(x)
r=(exp(x) + (x + 1)^101/1267650600228229401496703205376)*(exp(x) + (101*(x + 1)^100)/1267650600228229401496703205376) - exp(x)/100 - (101*(x + 1)^99)/1267650600228229401496703205376;
%r=x*(2*x - 1)*(x - 1) + 2;
%r=((exp(x) + (101*(x + 1)^100)/1267650600228229401496703205376)*(exp(x) + (x + 1)^101/1267650600228229401496703205376 + 1))/(x + 1) - (101*(x + 1)^99)/1267650600228229401496703205376 - exp(x)/100;
%r=2*sin(x) + 2*cos(x)*(x - 1) + 2*cos(x)*(x - 3/10) - cos(x)*(sin(x)*(x - 1) + sin(x)*(x - 3/10) + cos(x)*(x - 1)*(x - 3/10)) - sin(x)*(x - 1)*(x - 3/10);
%r=(exp(x) + (101*(x + 1)^100)/1267650600228229401496703205376)/(x + 1) - (101*(x + 1)^99)/1267650600228229401496703205376 - exp(x)/100;
%r=((exp(x) + (x + 1)^101/1267650600228229401496703205376)*(exp(x) + (101*(x + 1)^100)/1267650600228229401496703205376))/(x + 1) - (63125*(x + 1)^99)/79228162514264337593543950336 - 100*exp(x);
%r=2 - x*(2*x - 1);
end

function [r]=fxexact(x)
%r=(exp(x) + (2525*(x + 1)^99)/316912650057057350374175801344)/(x + 1) - (9999*(x + 1)^98)/1267650600228229401496703205376 - (exp(x) + (101*(x + 1)^100)/1267650600228229401496703205376)/(x + 1)^2 - exp(x)/100;
%r=1 - 4*x;
end

function [r]=uexact(x)
%r=-sin(x)*(x - 1)*(x - 3/10);
r=exp(x)+2^(-100)*(1+x)^(1+100);
%r=x*(1-x);
end

function [r]=uxexact(x)
r=1 - 2*x;
%r=exp(x) + (101*(x + 1)^100)/1267650600228229401496703205376;
end

function [r]=aexact(x)
%r=cos(x);
%r=1/(1+x);
%r=x;
r=1;
end

function [r]=axexact(x)
r=0;
%r=-1/(x + 1)^2;
%r=1;
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
%res=-exp(-ri)*(P_h_1_i(ri,ai,e));
end

function res=A2_2(ri,ai,a1i,e)
res=exp(-ri)*(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))-exp(-ri)*(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e));
%res=exp(-ri)*(P_h_1_i(ri,ai,e))-exp(-ri)*(P_p_1_i(ri,ai,e));
end

function res=A3_2(ri,ai,a1i,e)
res=exp(-ri)*(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e));
%res=exp(-ri)*(P_p_1_i(ri,ai,e));
end

function res=Y0_2(ri,ai,a1i,e)
res=(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e))*(Q_h_1_0i(ri,ai,e)+Q_h_2_0i(ri,ai,a1i,e))-(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))*(Q_p_1_0i(ri,ai,e)+Q_p_2_0i(ri,ai,a1i,e));
%res=(P_p_1_i(ri,ai,e))*(Q_h_1_0i(ri,ai,e))-(P_h_1_i(ri,ai,e))*(Q_p_1_0i(ri,ai,e));
end

function res=Y1_2(ri,ai,a1i,e)
res=(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e))*Q_h_2_1i(ri,ai,e)-(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))*Q_p_2_1i(ri,ai,e);
%res=(P_p_1_i(ri,ai,e))*Q_h_2_1i(ri,ai,e)-(P_h_1_i(ri,ai,e))*Q_p_2_1i(ri,ai,e);
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

function [res]=mynorm(A,B,h)
T=(A-B).*(A-B);
res=sqrt(sum(sum(T)))*h;
end