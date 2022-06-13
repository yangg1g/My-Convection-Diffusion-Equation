clear
tic
[p,e,u,x,y,k,et]=Solve(128,100000,1e-6);
toc
s=size(et);
myet=et(2:s(2))./et(1:s(2)-1);
surf(x,y,u);

function [p,e,u,x,y,k,et]=Solve(n,kmax,ep)
    h=1/n;
    u=zeros(n+1,n+1);
    u0=zeros(n+1,n+1);
    x=0:h:1;
    y=0:h:1;
%     for i=1:n+1
%         for j=1:n+1
%             u0(i,j)=uexact(x(i),y(j));
%         end
%     end
    a=zeros(n+1,n+1);
    b=zeros(n+1,n+1);
    f=zeros(n+1,n+1);
    for i=1:n+1
        for j=1:n+1
            a(i,j)=aexact(x(i),y(j));
            b(i,j)=bexact(x(i),y(j));
            f(i,j)=fexact(x(i),y(j));
            p(i,j)=uexact(x(i),y(j));
        end
    end
    for k=1:kmax
        for i=2:n
            for j=2:n
                A1=-u0(i,j+1)-u0(i,j-1)-u0(i+1,j)-u0(i-1,j)+4*u0(i,j);
                A2=h*a(i,j)*(u0(i+1,j)-u0(i-1,j))/2;
                A3=h*b(i,j)*(u0(i,j+1)-u0(i,j-1))/2;
                u(i,j)=u0(i,j)-1/4*(A1+A2+A3-h*h*f(i,j));
            end
        end
%         u(:,1)=u(:,2);
%         u(:,n+1)=u(:,n);
        et(k)=mynorm(u,p,h);
        if(et(k)<ep)
            break;
        end
        u0=u;
    end
    e=abs(u-p);
end

function [r]=aexact(x,y)
r=x+y+1;
%r=1;
end

function [r]=bexact(x,y)
r=x+y+1;
%r=1;
end

function [u]=uexact(x,y)
%u=x*(1-x)*y*(1-y);
% u=sin(x*pi)*cos(y*pi);
r=exp(x)+2^(-1/1e-2)*(1+x)^(1+1/1e-2);
end

function [f]=fexact(x,y)
% f=-2*pi^2*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) - 2*pi^2*exp(pi*(x + y))*cos(pi*y)*sin(pi*x);
% f=2*pi*pi*sin(x*pi)*cos(y*pi);
%f=2*x^3*y - x^3 + 4*x^2*y^2 - 3*x^2*y - 2*x^2 + 2*x*y^3 - 3*x*y^2 - 2*x*y + 3*x - y^3 - 2*y^2 + 3*y;
r=(exp(x) + (2525*(x + 1)^99)/316912650057057350374175801344)/(x + 1) - (9999*(x + 1)^98)/1267650600228229401496703205376 - (exp(x) + (101*(x + 1)^100)/1267650600228229401496703205376)/(x + 1)^2 - exp(x)/100;
end

function [res]=mynorm(A,B,h)
T=(A-B).*(A-B);
res=sqrt(sum(sum(T)))*h;
end