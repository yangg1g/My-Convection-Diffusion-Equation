clc
clear

[p,e,x,y,u,k,et]=Solve(128,1e-6);
surf(x,y,u);
mynorm(e,0,1/128)
n=4;
i=1;
while n<=128
    [p,e,x,y,u,k,et]=Solve(n,1e-6);
    E(i)=mynorm(e,0,1/n);
    n=n*2;
    i=i+1;
end
i=2;
while i<=6
    mye(i-1)=log(E(i-1)/E(i))/log(2);
    i=i+1;
end
function [p,e,x,y,u0,i,et]=Solve(n,ep)
t=1/n;
h=1/n;
r=t/h^2;
x=0:h:1;
y=0:h:1;
u0=zeros(n+1,n+1);
for i=1:n+1
    for j=1:n+1
        u0(i,j)=uexact(x(i),y(j));
    end
end
f=zeros(n+1,n+1);
for i=1:n+1
    for j=1:n+1
        f(i,j)=fexact(x(i),y(j));
    end
end
u_=zeros(n+1,n+1);
u=zeros(n+1,n+1);
a1=-r/2*ones(1,n-2);
b1=(1+r)*ones(1,n-1);
a2=-1/2*r*ones(1,n);
c2=a2;
a2(n)=-1;c2(1)=-1;
b2=(1+r)*ones(1,n+1);
b2(1)=1;b2(n+1)=1;
% L2=zeros(n-1,n-1);
% L2(1,1)=-2;
% for i=2:n-1
%     L2(i,i)=(1+r);
%     L2(i,i-1)=-r/2;
%     L2(i-1,i)=-r/2;
% end
% L1=zeros(n-1,n-1);
% L1(1,1)=-2;
% for i=2:n-1
%     L1(i,i)=(1-r);
%     L1(i,i-1)=r/2;
%     L1(i-1,i)=r/2;
% end
% L2_=zeros(n+1,n+1);
% for i=2:n
%     L2_(i,i)=(1+r);
%     L2_(i,i-1)=-r/2;
%     L2_(i,i+1)=-r/2;
% end
% L2_(1,1)=1;
% L2_(n+1,n+1)=1;
% L2_(1,2)=-1;
% L2_(n+1,n-1)=-1;

for i=1:n
    for j=2:n
        for k=2:n
            d1(k-1)=1/2*r*(u0(k+1,j)-2*u0(k,j)+u0(k-1,j))+u0(k,j)-1/2*t*f(k,j);
            % d(k-1)=1/2*r*deltaxx(u,k,j)+u(k,j)-1/2*t*fexact(x(k),y(j));
        end
        u(2:n,j)=ChaseMethod(a1,b1,a1,d1);
        % u2(j,2:n)=(L2^-1)*d';
    end
    u(:,1)=u(:,2);
    u(:,n+1)=u(:,n);
    for k=2:n
        d2(1)=0;d2(n+1)=0;
        for j=2:n
            d2(j)=1/2*r*(u(k,j+1)-2*u(k,j)+u(k,j-1))+u(k,j)-1/2*t*f(k,j);
        end
        u0(k,:)=ChaseMethod(a2,b2,c2,d2);
        % u(k,:)=(L2_^-1)*dd';
    end
    et(i)=mynorm(u0,u_,h);
    if(et(i)<ep)
        break;
    end
    u_=u0;
end
u=u0;
for i=1:n+1
    for j=1:n+1
        p(i,j)=uexact(x(i),y(j));
        e(i,j)=abs(u(i,j)-p(i,j));
    end
end 
end

function [x]=ChaseMethod(a,b,c,d)
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

function [u]=uexact(x,y)
u=sin(x*pi)*cos(y*pi);
end

function [f]=fexact(x,y)
f=sin(x*pi)*cos(y*pi)*-pi*pi*2;
end

function [r]=deltaxx(U,x,y)
r=U(x+1,y)-2*U(x,y)+U(x-1,y);
end

function [r]=deltayy(U,x,y)
r=U(x,y+1)-2*U(x,y)+U(x,y-1);
end

function [r]=delta(U,x,y)
r=deltaxx(U,x,y)+deltayy(U,x,y);
end

function [res]=mynorm(A,B,h)
T=(A-B).*(A-B);
res=sqrt(sum(sum(T)))*h;
end