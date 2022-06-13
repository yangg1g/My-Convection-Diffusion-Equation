clear
% [p,el,x,y,u,i,et]=Solve(128,1e-6,100,10000);
% surf(x,y,u);
% tic
% [p,e,x,y,u,k,et]=Solve(128,1e-6);
% toc
% e
% surf(x,y,u);
n=256;
vl=1;
vr=1e+10;
[p,el,x,y,u,i,et]=Solve(n,1e-6,vl,10);
[p,er,x,y,u,i,et]=Solve(n,1e-6,vr,10);
while(1)
    vm=exp((log(vl)+log(vr))/2);
    [p,em,x,y,u,i,et]=Solve(n,1e-6,vm,10);
    if(el<er)
        vr=vm;
        er=em;
    else
        vl=vm;
        el=em;
    end
    if(vr-vl<1)
        break;
    end
end

[p,e,x,y,u,k,et]=Solve(n,1e-6,vm,10000);
e
vm
k
i=1;
while n<=128
    [p,e,x,y,u,k(i),et]=Solve(n,1e-6,1,10000);
    E(i)=mynorm(e,0,1/n);
    n=n*2;
    i=i+1;
end
i=2;
fprintf("|%d|%d||\n",2^(i+1),E(i-1));
while i<=6
    mye(i-1)=log(E(i-1)/E(i))/log(2);
    fprintf("|%d|%d|%d|\n",2^(i+2),E(i),mye(i-1));
    i=i+1;
end

function [p,e,x,y,u,i,et]=Solve(n,ep,t,k)
h=1/n;
% t=100;
x=0:h:1;
y=0:h:1;
a=zeros(n+1,n+1);
b=zeros(n+1,n+1);
f=zeros(n+1,n+1);
u0=ones(n+1,n+1)*1e-5;
u_=zeros(n+1,n+1);
ux=zeros(n+1,n+1);
uy=zeros(n+1,n+1);
p=zeros(n+1,n+1);
e=0.01;
for i=1:n+1
    for j=1:n+1
        a(i,j)=aexact(x(i),y(j));
        b(i,j)=bexact(x(i),y(j));
        f(i,j)=fexact(x(i),y(j));
        %u0(i,j)=uexact(x(i),y(j));
        ux(i,j)=uxexact(x(i),y(j));
        uy(i,j)=uyexact(x(i),y(j));
        p(i,j)=uexact(x(i),y(j));
    end
end
for i=1:n+1
    u0(1,i)=p(1,i);
    u0(n+1,i)=p(n+1,i);
    u0(i,1)=p(i,1);
    u0(i,n+1)=p(i,n+1);
end
u=u0;
Y=zeros(1,n-1);

A1=zeros(1,n-1);
A2=zeros(1,n-1);
A3=zeros(1,n-1);

for i=1:k
    for k=2:n
        for j=2:n
            A1x=-e-h*(a(j,k))/2;
            A2x=2*e;
            A3x=-e+h*(a(j,k))/2;
        
            A1y=-e-h*(b(j,k))/2;
            A2y=2*e;
            A3y=-e+h*(b(j,k))/2;
            
%             A1x=-1-h*(u0(j,k))/2;
%             A2x=2;
%             A3x=-1+h*(u0(j,k))/2;
%         
%             A1y=-1-h*(u0(j,k))/2;
%             A2y=2;
%             A3y=-1+h*(u0(j,k))/2;
            fi=h*h*f(j,k);
            %fi=h*h*(f(j,k)-u0(j,k)*(ux(j,k)));
            %fi=h*h*(f(j,k)-u0(j,k)*(u0(j+1,k)-u0(j-1,k))/2/h);
            
            Y(j-1)=u0(j,k)-t*(A1y*u0(j,k-1)+A2y*u0(j,k)+A3y*u0(j,k+1))+t*fi;
            A1(j-1)=t*A1x;
            A2(j-1)=1+t*A2x;
            A3(j-1)=t*A3x;
            
        end
        u(2:n,k)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y);
    end
    mynorm(u,u0,1/128);
    for j=2:n
        for k=2:n
            A1x=-e-h*(a(j,k))/2;
            A2x=2*e;
            A3x=-e+h*(a(j,k))/2;
        
            A1y=-e-h*(b(j,k))/2;
            A2y=2*e;
            A3y=-e+h*(b(j,k))/2;
            
%             A1x=-1-h*(u(j,k))/2;
%             A2x=2;
%             A3x=-1+h*(u(j,k))/2;
%         
%             A1y=-1-h*(u(j,k))/2;
%             A2y=2;
%             A3y=-1+h*(u(j,k))/2;
            fi=h*h*f(j,k);
            %fi=h*h*(f(j,k)-u(j,k)*(uy(j,k)));
            %fi=h*h*(f(j,k)-u(j,k)*(u(j,k+1)-u(j,k-1))/2/h);
            
            Y(k-1)=u(j,k)-t*(A1x*u(j-1,k)+A2x*u(j,k)+A3x*u(j+1,k))+t*fi;
            A1(k-1)=t*A1y;
            A2(k-1)=1+t*A2y;
            A3(k-1)=t*A3y;
        end
        u0(j,2:n)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y);
    end
    et(i)=mynorm(u0,u_,h);
    if(et(i)<ep)
        break;
    end
    if(isnan(et(i)))
        u0=u_;
        break;
    end
    u_=u0;
end
u=u0;
e=mynorm(u,p,h);
end

function [r]=aexact(x,y)
%r=exp(x+y);
r=x+1;
%r=x+y+1;
end

function [r]=bexact(x,y)
%r=exp(x+y);
%r=exp(3*x+y);
r=x+2*y+1;
%r=x+y+1;
end

function [r]=fexact(x,y)
r=(exp(x) + (x + 1)^101/1267650600228229401496703205376)*(exp(y) + (101*(y + 1)^100)/1267650600228229401496703205376)*(x + 2*y + 1) - ((exp(x) + (2525*(x + 1)^99)/316912650057057350374175801344)*(exp(y) + (y + 1)^101/1267650600228229401496703205376))/100 - ((exp(x) + (x + 1)^101/1267650600228229401496703205376)*(exp(y) + (2525*(y + 1)^99)/316912650057057350374175801344))/100 + (exp(x) + (101*(x + 1)^100)/1267650600228229401496703205376)*(exp(y) + (y + 1)^101/1267650600228229401496703205376)*(x + 1);
%r=10*pi^2*sin(3*pi*x)*sin(pi*y) + pi*cos(pi*y)*sin(3*pi*x)*(x + 2*y + 1) + 3*pi*cos(3*pi*x)*sin(pi*y)*(x + 1);
%r=10*pi^2*sin(3*pi*x)*sin(pi*y) + 3*pi*cos(3*pi*x)*sin(3*pi*x)*sin(pi*y)^2 + pi*cos(pi*y)*sin(3*pi*x)^2*sin(pi*y);
%r=2*pi^2*sin(pi*x)*sin(pi*y) + pi*cos(pi*x)*sin(pi*x)*sin(pi*y)^2 + pi*cos(pi*y)*sin(pi*x)^2*sin(pi*y);
%r=x*y*(x*y*(x - 1) + x*(x - 1)*(y - 1))*(x - 1)*(y - 1) - 2*y*(y - 1) - 2*x*(x - 1) + x*y*(x*y*(y - 1) + y*(x - 1)*(y - 1))*(x - 1)*(y - 1);
%r=2*sin(pi*y) + x*sin(pi*y)*(x*sin(pi*y) + sin(pi*y)*(x - 1))*(x - 1) - x*pi^2*sin(pi*y)*(x - 1) + x^2*pi*cos(pi*y)*sin(pi*y)*(x - 1)^2;
%r=exp(pi*(x + y))*sin(pi*x)*sin(pi*y)*(pi*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) + pi*exp(pi*(x + y))*sin(pi*x)*sin(pi*y)) - 2*pi^2*exp(pi*(x + y))*cos(pi*y)*sin(pi*x) - 2*pi^2*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) + exp(pi*(x + y))*sin(pi*x)*sin(pi*y)*(pi*exp(pi*(x + y))*cos(pi*y)*sin(pi*x) + pi*exp(pi*(x + y))*sin(pi*x)*sin(pi*y));
%r=(exp(x)*(exp(x) - 1)*(exp(y) - 1)*(exp(y) - 3060513257434037/1125899906842624) + exp(x)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624))*(exp(x) - 1)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624) - 2*exp(2*x)*(exp(y) - 1)*(exp(y) - 3060513257434037/1125899906842624) - exp(x)*(exp(x) - 1)*(exp(y) - 1)*(exp(y) - 3060513257434037/1125899906842624) - exp(y)*(exp(x) - 1)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624) - exp(x)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624) - exp(y)*(exp(x) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624) - 2*exp(2*y)*(exp(x) - 1)*(exp(x) - 3060513257434037/1125899906842624) + (exp(y)*(exp(x) - 1)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624) + exp(y)*(exp(x) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624))*(exp(x) - 1)*(exp(y) - 1)*(exp(x) - 3060513257434037/1125899906842624)*(exp(y) - 3060513257434037/1125899906842624);
end

function [r]=uexact(x,y)
r=(exp(x)+2^(-100)*(1+x)^(1+100))*(exp(y)+2^(-100)*(1+y)^(1+100));
%r=sin(3*pi*x)*sin(pi*y);
%r=sin(3*pi*x)*sin(pi*y);
%r=x*(1-x)*y*(1-y);
%r=-x*sin(pi*y)*(x - 1);
%r=exp(pi*(x+y))*sin(x*pi)*sin(y*pi);
%r=(exp(x)-1)*(exp(x)-exp(1))*(exp(y)-1)*(exp(y)-exp(1));
end

function [r]=uxexact(x,y)
r=- x*sin(pi*y) - sin(pi*y)*(x - 1);
end

function [r]=uyexact(x,y)
r=-x*pi*cos(pi*y)*(x - 1);
end

function [x]=zhuiganfa(a,b,c,d)
r=size(a);
m=r(2);
r=size(b);
n=r(2);
if size(a)~=size(c)|m~=n-1|size(b)~=size(d)
    error('');
end
u=zeros(1,n);
y=zeros(1,n);
x=zeros(1,n);
u(1)=b(1);
l=zeros(1,n-1);
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