clear
n=128;
[p1,p2,em,x,y,u,v,k,et]=Solve(n,1e-6,10,10000);
em
k
vl=1;
vr=1000;
nl=0;
nr=0;
[p1,p2,el,x,y,u,v,k,et]=Solve(n,1e-6,vl,10);
[p1,p2,er,x,y,u,v,k,et]=Solve(n,1e-6,vr,10);
while(1)
    vm=exp((log(vl)+log(vr))/2);
    %vm=(vl+vr)/2;
    [p1,p2,em,x,y,u,v,k,et]=Solve(n,1e-6,vm,10);
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

[p1,p2,e,x,y,u,v,i,et]=Solve(n,1e-6,vm,10000);
e
vm
i
% surf(x,y,u);


function [p1,p2,e,x,y,u,v,i,et]=Solve(n,ep,t,k)
h=1/n;
% t=100;
x=0:h:1;
y=0:h:1;
a=zeros(n+1,n+1);
b=zeros(n+1,n+1);
fu=zeros(n+1,n+1);
fv=zeros(n+1,n+1);
u0=ones(n+1,n+1)*1e-5;
u_=zeros(n+1,n+1);
v0=ones(n+1,n+1)*1e-5;
v_=zeros(n+1,n+1);
ux=zeros(n+1,n+1);
uy=zeros(n+1,n+1);
p1=zeros(n+1,n+1);
p2=zeros(n+1,n+1);
e=1;
for i=1:n+1
    for j=1:n+1
        a(i,j)=aexact(x(i),y(j));
        b(i,j)=bexact(x(i),y(j));
        fu(i,j)=fuexact(x(i),y(j));
        fv(i,j)=fvexact(x(i),y(j));
        %u0(i,j)=uexact(x(i),y(j));
        %v0(i,j)=vexact(x(i),y(j));
        ux(i,j)=uxexact(x(i),y(j));
        uy(i,j)=uyexact(x(i),y(j));
        p1(i,j)=uexact(x(i),y(j));
        p2(i,j)=vexact(x(i),y(j));
    end
end

for i=1:n+1
    u0(1,i)=p1(1,i);
    u0(n+1,i)=p1(n+1,i);
    u0(i,1)=p1(i,1);
    u0(i,n+1)=p1(i,n+1);
    
    v0(1,i)=p2(1,i);
    v0(n+1,i)=p2(n+1,i);
    v0(i,1)=p2(i,1);
    v0(i,n+1)=p2(i,n+1);
end

Y1=zeros(1,n-1);
Y2=zeros(1,n-1);
u=u0;
v=v0;

A1=zeros(1,n-1);
A2=zeros(1,n-1);
A3=zeros(1,n-1);

for i=1:k
    for k=2:n
        for j=2:n
            aix=u0(j,k);           
            A1x=-e-h*aix/2;
            A2x=2*e;
            A3x=-e+h*aix/2;
            
            aiy=v0(j,k);         
            A1y=-e-h*aiy/2;
            A2y=2*e;
            A3y=-e+h*aiy/2;

            fiu=h*h*fu(j,k);
            fiv=h*h*fv(j,k);

            Y1(j-1)=u0(j,k)-t*(A1y*u0(j,k-1)+A2y*u0(j,k)+A3y*u0(j,k+1))+t*fiu;
            Y2(j-1)=v0(j,k)-t*(A1y*v0(j,k-1)+A2y*v0(j,k)+A3y*v0(j,k+1))+t*fiv;
            
            A1(j-1)=t*A1x;
            A2(j-1)=1+t*A2x;
            A3(j-1)=t*A3x;
            
        end
        Y1(1)=Y1(1)-A1(1)*u0(1,k);
        Y1(n-1)=Y1(n-1)-A3(n-1)*u0(n+1,k);
        Y2(1)=Y2(1)-A1(1)*v0(1,k);
        Y2(n-1)=Y2(n-1)-A3(n-1)*v0(n+1,k);
        
        u(2:n,k)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y1);
        v(2:n,k)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y2);
    end
%     mynorm(u,u0,1/128);
    for j=2:n
        for k=2:n
            aix=u(j,k);          
            A1x=-e-h*aix/2;
            A2x=2*e;
            A3x=-e+h*aix/2;
            
            aiy=v(j,k);            
            A1y=-e-h*aiy/2;
            A2y=2*e;
            A3y=-e+h*aiy/2;

            fiu=h*h*fu(j,k);
            fiv=h*h*fv(j,k);
            
            Y1(k-1)=u(j,k)-t*(A1x*u(j-1,k)+A2x*u(j,k)+A3x*u(j+1,k))+t*fiu;
            Y2(k-1)=v(j,k)-t*(A1x*v(j-1,k)+A2x*v(j,k)+A3x*v(j+1,k))+t*fiv;
            
            A1(k-1)=t*A1y;
            A2(k-1)=1+t*A2y;
            A3(k-1)=t*A3y;
        end
        Y1(1)=Y1(1)-A1(1)*u(j,1);
        Y1(n-1)=Y1(n-1)-A3(n-1)*u(j,n+1);
        Y2(1)=Y2(1)-A1(1)*v(j,1);
        Y2(n-1)=Y2(n-1)-A3(n-1)*v(j,n+1);
        
        u0(j,2:n)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y1);
        v0(j,2:n)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y2);
    end
    et(i)=(mynorm(u0,u_,h)+mynorm(v0,v_,h))/2;
    if(et(i)<ep)
        break;
    end
    if(isnan(et(i)))
        u0=u_;
        break;
    end
    u_=u0;
    v_=v0;
end
u=u0;
v=v0;
e=(mynorm(u,p1,h)+mynorm(v,p2,h))/2;
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

function [r]=fuexact(x,y)
r=(exp(y - x) + (101*(y + 1)^100)/1267650600228229401496703205376)*(exp(x) + (x + 1)^101/1267650600228229401496703205376)*(exp(y) + (y + 1)^101/1267650600228229401496703205376) - (101*(y + 1)^99)/1267650600228229401496703205376 - exp(y - x)*(exp(y - x) + (y + 1)^101/1267650600228229401496703205376) - exp(y - x)/50;
end

function [r]=fvexact(x,y)
r=(exp(y - x) + (y + 1)^101/1267650600228229401496703205376)*(exp(x) + (101*(x + 1)^100)/1267650600228229401496703205376)*(exp(y) + (y + 1)^101/1267650600228229401496703205376) - ((exp(x) + (2525*(x + 1)^99)/316912650057057350374175801344)*(exp(y) + (y + 1)^101/1267650600228229401496703205376))/100 - ((exp(x) + (x + 1)^101/1267650600228229401496703205376)*(exp(y) + (2525*(y + 1)^99)/316912650057057350374175801344))/100 + (exp(x) + (x + 1)^101/1267650600228229401496703205376)^2*(exp(y) + (y + 1)^101/1267650600228229401496703205376)*(exp(y) + (101*(y + 1)^100)/1267650600228229401496703205376);
end

function [r]=uexact(x,y)
r=exp(y - x) + (y + 1)^101/1267650600228229401496703205376;
end

function [r]=vexact(x,y)
r=(exp(x)+2^(-100)*(1+x)^(1+100))*(exp(y)+2^(-100)*(1+y)^(1+100));
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