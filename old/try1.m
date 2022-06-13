clear

n=32;
h=1/n;
t=1/n;
x=0:h:1;
y=0:h:1;
a=zeros(n+1,n+1);
b=zeros(n+1,n+1);
f=zeros(n+1,n+1);
u0=zeros(n+1,n+1);
u_=zeros(n+1,n+1);
e=1;
for i=1:n+1
    for j=1:n+1
        a(i,j)=aexact(x(i),y(j));
        b(i,j)=bexact(x(i),y(j));
        f(i,j)=fexact(x(i),y(j));
        %u0(i,j)=uexact(x(i),y(j));
        u_(i,j)=u0(i,j);
    end
end
A1=zeros(1,n-1);
A2=zeros(1,n-1);
A3=zeros(1,n-1);
Y=zeros(1,n-1);
u=zeros(n+1,n+1);
for i=1:n
    for j=2:n
        for k=2:n
            ai=a(j,k);
            ri=h*ai/e;
            fi=f(j,k)+e*(u0(j+1,k)+u0(j-1,k)-2*u0(j,k))/h/h-b(j,k)*(u0(j+1,k)-u0(j-1,k))/2/h;
            A1(k-1)=A1_1(ri,ai,e);
            A2(k-1)=A2_1(ri,ai,e);
            A3(k-1)=A3_1(ri,ai,e);
            Y(k-1)=Y0_1(ri,ai,e)*fi;
        end
        
        u(2:n,j)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y);
    end
    for k=2:n
        for j=2:n
            ai=b(j,k);
            ri=h*ai/e;
            fi=f(j,k)+e*(u(j,k+1)+u(j,k-1)-2*u(j,k))/h/h-a(j,k)*(u(j,k+1)-u(j,k-1))/2/h;
            A1(j-1)=A1_1(ri,ai,e);
            A2(j-1)=A2_1(ri,ai,e);
            A3(j-1)=A3_1(ri,ai,e);
            Y(j-1)=Y0_1(ri,ai,e)*fi;
        end
        u0(k,2:n)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y);
    end
end
u=u0;
surf(x,y,u);

function [r]=aexact(x,y)
%r=x+y+1;
r=1;
end

function [r]=bexact(x,y)
%r=x+y+1;
r=1;
end

function [r]=fexact(x,y)
% r=2*x^3*y - x^3 + 4*x^2*y^2 - 3*x^2*y - 2*x^2 + 2*x*y^3 - 3*x*y^2 - 2*x*y + 3*x - y^3 - 2*y^2 + 3*y;
r=2*x^2*y - 3*x^2 + 2*x*y^2 - 4*x*y + 3*x - 3*y^2 + 3*y;
end

function [r]=uexact(x,y)
r=x*(1-x)*y*(1-y);
end

function u=Solve1_(e,a,f,n)
h=1/n;
for i=1:n-1
    ai=a(i);
    ri=h*ai/e;
    fi=f(i);
    A1(i)=A1_1(ri,ai,e);
    A2(i)=A2_1(ri,ai,e);
    A3(i)=A3_1(ri,ai,e);
    Y(i)=Y0_1(ri,ai,e)*fi;
end
A1(1)=[];
A3(n-1)=[];
u=[0,zhuiganfa(A1,A2,A3,Y),0];
end

function res=P_h_1_i(ri,ai,e)
res=e/ai*(1-exp(-ri));
end

function res=P_p_1_i(ri,ai,e)
res=e/ai*(exp(-2*ri)-exp(-ri));
end

function res=Q_h_1_0i(ri,ai,e)
res=-e/ai^2*(1-exp(-ri)-ri*exp(-ri));
end

function res=Q_p_1_0i(ri,ai,e)
res=-e/ai^2*(exp(-2*ri)-exp(-ri)+ri*exp(-ri));
end

function res=A1_1(ri,ai,e)
res=-exp(-ri)*(P_h_1_i(ri,ai,e));
end

function res=A2_1(ri,ai,e)
res=exp(-ri)*(P_h_1_i(ri,ai,e)-P_p_1_i(ri,ai,e));
end

function res=A3_1(ri,ai,e)
res=exp(-ri)*(P_p_1_i(ri,ai,e));
end

function res=Y0_1(ri,ai,e)
res=P_p_1_i(ri,ai,e)*Q_h_1_0i(ri,ai,e)-P_h_1_i(ri,ai,e)*Q_p_1_0i(ri,ai,e);
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