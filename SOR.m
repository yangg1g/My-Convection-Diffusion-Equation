tic
[p,e,u,x,y,k,et]=Solve(256,10000,1e-6);
toc
% s=size(et);
% myet=et(2:s(2))./et(1:s(2)-1);
e
k
surf(x,y,p-u);
n=8;
i=1;
while n<=128
    [p,e,x,y,u,k,et]=Solve(n,10000,1e-6);
    E(i)=e;
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
    w=1;
    w=2/(1+sqrt(1-cos(pi*h)*cos(pi*h)));
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
            u(i,j)=u0(i,j)+w*(h*h*f(i,j)+u0(i,j+1)+u(i,j-1)+u0(i+1,j)+u(i-1,j)-4*u0(i,j))/4;
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
        end
    end
    e=mynorm(u,p,h);
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