clear
clc
[p,e,u,x,y,k,et]=Solve(64,10000,1e-6);
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
%         u(i,1)=sin(pi*y(i));
%         u(i,m+1)=exp(1)*exp(1)*sin(pi*y(i));
%     end
    for i=1:n
        for j=1:n
            f(i,j)=fexact(x(i),y(j));
        end
    end
    for k=1:kmax
        t=0;
%         for i=2:n
%             for j=2:n
%                 u(i,j)=h*h*f(i,j)/4+(u(i,j+1)+u(i,j-1)+u(i+1,j)+u(i-1,j))/4;
%             end
%         end
        ij=[1,1];
        while(ij(2)<n)
            i=ij(1)+1;
            j=ij(2)+1;
            u(i,j)=h*h*f(i,j)/4+(u0(i,j+1)+u(i,j-1)+u0(i+1,j)+u(i-1,j))/4;
            ij=nextij(ij,n-1);
        end
%         u(:,1)=u(:,2);
%         u(:,n+1)=u(:,n);
        et(k)=mynorm(u,u0,h);
        if(et(k)<ep)
            break;
        end
        u0=u;
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

function ij=nextij(ij,n)
ij(1)=ij(1)-1;
ij(2)=ij(2)+1;
if(ij(2)>n)
    ij(2)=ij(1)+ij(2)+1-n;
    ij(1)=n;
elseif(ij(1)==0)
    ij(1)=ij(2);
    ij(2)=1;
else 
end
end

function [res]=mynorm(A,B,h)
T=(A-B).*(A-B);
res=sqrt(sum(sum(T)))*h;
end