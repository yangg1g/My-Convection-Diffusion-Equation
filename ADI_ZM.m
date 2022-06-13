clear
tic
[p,e,x,y,u,k,et]=Solve(256,1e-6);
toc
surf(x,y,u);
mynorm(u,p,1/128)
k
n=8;
i=1;
while n<=128
    [p,e,x,y,u,k,et]=Solve(n,1e-6);
    E(i)=mynorm(e,0,1/n);
    K(i)=k;
    n=n*2;
    i=i+1;
end
i=2;
fprintf("%d & %d & & %d\\\\\n",2^(i+1),E(i-1),K(i-1));
while i<=6
    mye(i-1)=log(E(i-1)/E(i))/log(2);
    fprintf("%d & %d & %d & %d\\\\\n",2^(i+2),E(i),mye(i-1),K(i));
    i=i+1;
end
function [p,e,x,y,u,k,et]=Solve(n,ep)
h=1/n;
t=h*h/2/sin(pi*h);
%t=0.015;
r=t/h^2;

x=0:h:1;
y=0:h:1;
u0=zeros(n+1,n+1);
u_=zeros(n+1,n+1);
f=zeros(n+1,n+1);
u=zeros(n+1,n+1);
p=zeros(n+1,n+1);
e=zeros(n+1,n+1);
for i=1:n+1
    for j=1:n+1
        f(i,j)=fexact(x(i),y(j));
        p(i,j)=uexact(x(i),y(j));
        % u0(i,j)=uexact(x(i),y(j));
    end
end
a1=-r*ones(1,n-2);
b1=(1+2*r)*ones(1,n-1);
d=zeros(1,n-2);

for i=1:100000
    for j=2:n
        for k=2:n
            d(k-1)=r*(u0(k,j+1)-2*u0(k,j)+u0(k,j-1))+u0(k,j)+t*f(k,j);
        end
        u(2:n,j)=ChaseMethod(a1,b1,a1,d);
    end
    for k=2:n
        for j=2:n
            d(j-1)=r*(u(k+1,j)-2*u(k,j)+u(k-1,j))+u(k,j)+t*f(k,j);
        end
        u0(k,2:n)=ChaseMethod(a1,b1,a1,d);
    end
    et(i)=mynorm(u0,u_,h);
    if(et(i)<ep)
        break;
    end
    u_=u0;
end
u=u0;
k=i;
end

function [x]=ChaseMethod(a,b,c,d)
r=size(a);
m=r(2);
r=size(b);
n=r(2);
sd=size(d);
if size(a)~=size(c)|m~=n-1|n~=sd(2)
    error('');
end
l=zeros(1,n-1);
u=zeros(1,n);
y=zeros(1,n);
x=zeros(sd);
u(1)=b(1);
for i=2:n
    l(i-1)=a(i-1)/u(i-1);
    u(i)=b(i)-l(i-1)*c(i-1);
end
for j=1:sd(1)
    y(1)=d(j,1);
    for i=2:n
        y(i)=d(j,i)-l(i-1)*y(i-1);
    end
    x(j,n)=y(n)/u(n);
    for i=n-1:-1:1
        x(j,i)=y(i)/u(i);
        x(j,i)=(y(i)-c(i)*x(j,i+1))/u(i);
    end
end
end

function [u]=uexact(x,y)
u=exp(pi*(x+y))*sin(x*pi)*sin(y*pi);
% u=sin(x*pi)*sin(y*pi);
%u=exp(pi*(x + y))*sin(2*pi*y)*sin(pi*x);
end

function [f]=fexact(x,y)
f= - 2*pi^2*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) - 2*pi^2*exp(pi*(x + y))*cos(pi*y)*sin(pi*x);
% f=-2*pi^2*sin(pi*x)*sin(pi*y);
%f=2*pi^2*exp(pi*(x + y))*cos(pi*x)*sin(2*pi*y) + 4*pi^2*exp(pi*(x + y))*cos(2*pi*y)*sin(pi*x) - 3*pi^2*exp(pi*(x + y))*sin(pi*x)*sin(2*pi*y);
end

function [r]=deltaxx(U,x,y,h)
r=(U(x+1,y)-2*U(x,y)+U(x-1,y))/h/h;
end

function [r]=deltayy(U,x,y,h)
r=(U(x,y+1)-2*U(x,y)+U(x,y-1))/h/h;
end

function [r]=delta(U,x,y,h)
r=deltaxx(U,x,y,h)+deltayy(U,x,y,h);
end

function [res]=mynorm(A,B,h)
T=(A-B).*(A-B);
res=sqrt(sum(sum(T)))*h;
end