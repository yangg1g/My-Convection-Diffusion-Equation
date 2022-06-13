clear
tic
[p,e,x,y,u,k,et]=Solve(128,1e-6);
toc
mynorm(e,0,1/128)
surf(x,y,u);
n=4;
i=1;
while n<=256
    [p,e,x,y,u,k(i),et]=Solve(n,1e-8);
    E(i)=mynorm(e,0,1/n);
    n=n*2;
    i=i+1;
end
i=2;
while i<=7
    mye(i-1)=log(E(i-1)/E(i))/log(2);
    i=i+1;
end

function [p,e,x,y,u,i,et]=Solve(n,ep)
h=1/n;
t=100;
x=0:h:1;
y=0:h:1;
a=zeros(n+1,n+1);
b=zeros(n+1,n+1);
f=zeros(n+1,n+1);
u0=ones(n+1,n+1)*1e-5;
u_=zeros(n+1,n+1);
e=1;
for i=1:n+1
    for j=1:n+1
        a(i,j)=aexact(x(i),y(j));
        b(i,j)=bexact(x(i),y(j));
        f(i,j)=fexact(x(i),y(j));
        %u0(i,j)=uexact(x(i),y(j));
        p(i,j)=uexact(x(i),y(j));
    end
end
Y=zeros(1,n-1);
Yt=zeros(1,n-1);
u=zeros(n+1,n+1);

fi_=zeros(n-1,n-1);
A11_=zeros(n-1,1);
A12_=zeros(n-1,1);
A13_=zeros(n-1,1);
A21_=zeros(n-1,1);
A22_=zeros(n-1,1);
A23_=zeros(n-1,1);

syms sx
for i=1:10000
    for k=2:n
        for j=2:n
            A1x_=-1-h*u0(j,k)/2;
            A2x_=2;
            A3x_=-1+h*u0(j,k)/2;
            
            A11_(j-1)=t*A1x_;
            A12_(j-1)=1+t*A2x_;
            A13_(j-1)=t*A3x_;
            
            A1y=-1-h*u0(j,k)/2;
            A2y=2;
            A3y=-1+h*u0(j,k)/2;
            Y(j-1)=u0(j,k)-t*(A1y*u0(j,k-1)+A2y*u0(j,k)+A3y*u0(j,k+1))+t*h*h*f(j,k);
        end
        u(2:n,k)=zhuiganfa(A11_,A12_',A13_,Y);
    end
    for j=2:n
        for k=2:n
            A1y_=-1-h*(u(j,k))/2;
            A2y_=2;
            A3y_=-1+h*(u(j,k))/2;
            
            A21_(k-1)=t*A1y_;
            A22_(k-1)=1+t*A2y_;
            A23_(k-1)=t*A3y_;
            
            A1x=-1-h*(u(j,k))/2;
            A2x=2;
            A3x=-1+h*(u(j,k))/2;
            Y(k-1)=u(j,k)-t*(A1x*u(j-1,k)+A2x*u(j,k)+A3x*u(j+1,k))+t*h*h*f(j,k);
        end
        u0(j,2:n)=zhuiganfa(A21_,A22_',A23_,Y);
    end
    et(i)=mynorm(u0,u_,h);
    if(et(i)<ep)
        break;
    end
    u_=u0;
end
u=u0;
e=abs(u-p);
end

function [r]=aexact(x,y)
%r=exp(x+y);
r=x+y+1;
%r=x+1;
end

function [r]=bexact(x,y)
%r=exp(x+y);
%r=exp(3*x+y);
r=x+y+1;
%r=x+y+1;
end

function [r]=fexact(x,y)
%r=10*pi^2*sin(3*pi*x)*sin(pi*y) + 3*pi*cos(3*pi*x)*sin(3*pi*x)*sin(pi*y)^2 + pi*cos(pi*y)*sin(3*pi*x)^2*sin(pi*y);
%r=2*pi^2*sin(pi*x)*sin(pi*y) + pi*cos(pi*x)*sin(pi*x)*sin(pi*y)^2 + pi*cos(pi*y)*sin(pi*x)^2*sin(pi*y);
%r=x*y*(x*y*(x - 1) + x*(x - 1)*(y - 1))*(x - 1)*(y - 1) - 2*y*(y - 1) - 2*x*(x - 1) + x*y*(x*y*(y - 1) + y*(x - 1)*(y - 1))*(x - 1)*(y - 1);
%r=2*sin(pi*y) + x*sin(pi*y)*(x*sin(pi*y) + sin(pi*y)*(x - 1))*(x - 1) - x*pi^2*sin(pi*y)*(x - 1) + x^2*pi*cos(pi*y)*sin(pi*y)*(x - 1)^2;
%r=exp(pi*(x + y))*sin(pi*x)*sin(pi*y)*(pi*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) + pi*exp(pi*(x + y))*sin(pi*x)*sin(pi*y)) - 2*pi^2*exp(pi*(x + y))*cos(pi*y)*sin(pi*x) - 2*pi^2*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) + exp(pi*(x + y))*sin(pi*x)*sin(pi*y)*(pi*exp(pi*(x + y))*cos(pi*y)*sin(pi*x) + pi*exp(pi*(x + y))*sin(pi*x)*sin(pi*y));
r=(exp(x)*(exp(x) - 1)*(exp(y) - 1)*(exp(y) - 3060513257434037/1125899906842624) + exp(x)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624))*(exp(x) - 1)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624) - 2*exp(2*x)*(exp(y) - 1)*(exp(y) - 3060513257434037/1125899906842624) - exp(x)*(exp(x) - 1)*(exp(y) - 1)*(exp(y) - 3060513257434037/1125899906842624) - exp(y)*(exp(x) - 1)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624) - exp(x)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624) - exp(y)*(exp(x) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624) - 2*exp(2*y)*(exp(x) - 1)*(exp(x) - 3060513257434037/1125899906842624) + (exp(y)*(exp(x) - 1)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624) + exp(y)*(exp(x) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624))*(exp(x) - 1)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624);
end

function [r]=uexact(x,y)
%r=sin(3*pi*x)*sin(pi*y);
%r=sin(pi*x)*sin(pi*y);
%r=x*(1-x)*y*(1-y);
%r=-x*sin(pi*y)*(x - 1);
%r=exp(pi*(x+y))*sin(x*pi)*sin(y*pi);
r=(exp(x)-1)*(exp(x)-exp(1))*(exp(y)-1)*(exp(y)-exp(1));
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