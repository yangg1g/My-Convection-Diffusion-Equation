clear
clc
[p,e,u,x,y,k,et]=Solve(32,10000,1e-6);
s=size(et);
myet=et(2:s(2))./et(1:s(2)-1);
surf(x,y,u);

function [p,e,u,x,y,k,et]=Solve(n,kmax,ep)
    h=1/n;
    u=zeros(n+1,n+1);
    u0=zeros(n+1,n+1);
    Au0=zeros(n+1,n+1);
    r=zeros(n+1,n+1);
    p=zeros(n+1,n+1);
    Ar=zeros(n+1,n+1);
    Ap=zeros(n+1,n+1);
    x=0:h:1;
    y=0:h:1;
%     for i=1:n+1
%         for j=1:n+1
%             u0(i,j)=uexact(x(i),y(j));
%         end
%     end
    for i=1:n+1
        for j=1:n+1
            f(i,j)=fexact(x(i),y(j));
        end
    end
    for i=2:n
        for j=2:n
            Au0(i,j)=-(u0(i,j+1)+u0(i,j-1)+u0(i+1,j)+u0(i-1,j)-4*u0(i,j))/h/h;
        end
    end
    r0(2:n,2:n)=f(2:n,2:n)-Au0(2:n,2:n);
    p(2:n,2:n)=r0(2:n,2:n);
    for k=1:kmax
        for i=2:n
            for j=2:n
                Ap(i,j)=-(p(i,j+1)+p(i,j-1)+p(i+1,j)+p(i-1,j)-4*p(i,j))/h/h;
            end
        end
        a=sum(sum(r0(2:n,2:n).*r0(2:n,2:n)))/sum(sum(Ap(2:n,2:n).*p(2:n,2:n)));
        u(2:n,2:n)=u0(2:n,2:n)+a*p(2:n,2:n);
        r(2:n,2:n)=r0(2:n,2:n)-a*Ap(2:n,2:n);
        b=sum(sum(r(2:n,2:n).*r(2:n,2:n)))/sum(sum(r0(2:n,2:n).*r0(2:n,2:n)));
        p(2:n,2:n)=r(2:n,2:n)+b*p(2:n,2:n);
%         u(:,1)=u(:,2);
%         u(:,n+1)=u(:,n);
        et(k)=mynorm(u,u0,h);
        if(et(k)<ep)
            break;
        end
        u0=u;
        r0=r;
    end
    for i=1:n+1
        for j=1:n+1
            p(i,j)=uexact(x(i),y(j));
            e(i,j)=abs(u(i,j)-p(i,j));
        end
    end
end

function [u]=uexact(x,y)
u=exp(pi*(x+y))*sin(x*pi)*sin(y*pi);
% u=sin(x*pi)*cos(y*pi);
end

function [f]=fexact(x,y)
f=-2*pi^2*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) - 2*pi^2*exp(pi*(x + y))*cos(pi*y)*sin(pi*x);
% f=2*pi*pi*sin(x*pi)*cos(y*pi);
end

function [res]=mynorm(A,B,h)
T=(A-B).*(A-B);
res=sqrt(sum(sum(T)))*h;
end
