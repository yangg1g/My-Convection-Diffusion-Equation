clear

[p,e,x,y,u,k,et]=Solve(128,1e-10);
surf(x,y,u);
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
function [p,e,x,y,u0,k,et]=Solve(n,ep)
h=1/n;
t=h*h/2/sin(pi*h);
t=0.015;
r=t/h^2;

x=0:h:1;
y=0:h:1;
u0=zeros(n+1,n+1);
u_=zeros(n+1,n+1);
f=zeros(n+1,n+1);
for i=1:n+1
    for j=1:n+1
        f(i,j)=fexact(x(i),y(j));
        % u0(i,j)=uexact(x(i),y(j));
    end
end
u=zeros(n+1,n+1);
a1=-r*ones(1,n-2);
b1=(1+2*r)*ones(1,n-1);
d=zeros(1,n-2);
L2=zeros(n-1,n-1);
L2(1,1)=(1-2*r);
for i=2:n-1
    L2(i,i)=(1-2*r);
    L2(i,i-1)=r;
    L2(i-1,i)=r;
end
L1=zeros(n-1,n-1);
L1(1,1)=(1+2*r);
for i=2:n-1
    L1(i,i)=(1+2*r);
    L1(i,i-1)=-r;
    L1(i-1,i)=-r;
end
% rou=sqrt(2)-1;
% p=round(log(2/pi/h)/log((1+rou)/(1-rou)));
% a=sin(pi*h/2)*sin(pi*h/2);
% t1=rand(n-1,1);
% t2=(L1^-1)*(L2*t1);
% mynorm(t1,0,h)
% mynorm(t2,0,h)
% A=(L1^-1)*(L2);
% x=ones(n-1,1);
% max(abs(eig(A)))
% norm(A^1000*x)
% norm(A^1000)*norm(x)
for i=1:100000
    %k=mod(i,p)+1;
    %t=h*h/4/a*power((1-rou)/(1+rou),2*k-1);
    for j=2:n
        for k=2:n
            d(k-1)=r*(u0(k,j+1)-2*u0(k,j)+u0(k,j-1))+u0(k,j)-t*f(k,j);
            % d(k-1)=1/2*r*deltaxx(u,k,j)+u(k,j)-1/2*t*fexact(x(k),y(j));
        end
        % u(2:n,j)=ChaseMethod(a1,b1,a1,d);
        % u(2:n,j)=(L1^-1)*d';
    end
    u(2:n,2:n)=(L1^-1)*(L2*u0(2:n,2:n)'-t*f(2:n,2:n)')';
    
    for k=2:n
        for j=2:n
            d(j-1)=r*(u(k+1,j)-2*u(k,j)+u(k-1,j))+u(k,j)-t*f(k,j);
        end
        % u0(k,2:n)=ChaseMethod(a1,b1,a1,d);
        % u0(k,2:n)=(L1^-1)*d';
    end
    u0(2:n,2:n)=((L1^-1)*(L2*u(2:n,2:n)-t*f(2:n,2:n))')';
    et(i)=mynorm(u0,u_,h);
    if(et(i)<ep)
        break;
    end
    u_=u0;
    uu(i,:,:)=u0;
end
u=u0;
k=i;
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
% u=exp(pi*(x+y))*sin(x*pi)*sin(y*pi);
% u=sin(x*pi)*sin(y*pi);
u=exp(pi*(x + y))*sin(2*pi*y)*sin(pi*x);
end

function [f]=fexact(x,y)
% f=2*pi^2*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) + 2*pi^2*exp(pi*(x + y))*cos(pi*y)*sin(pi*x);
% f=-2*pi^2*sin(pi*x)*sin(pi*y);
f=2*pi^2*exp(pi*(x + y))*cos(pi*x)*sin(2*pi*y) + 4*pi^2*exp(pi*(x + y))*cos(2*pi*y)*sin(pi*x) - 3*pi^2*exp(pi*(x + y))*sin(pi*x)*sin(2*pi*y);
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
res=norm(A-B,'fro')*h;
end