clear
tic
[p,e,x,y,u,k,et]=Solve(128,1e-6);
toc
surf(x,y,u);
n=4;
i=1;
while n<=128
    [p,e,x,y,u,k(i),et]=Solve(n,1e-6);
    E(i)=mynorm(e,0,1/n);
    n=n*2;
    i=i+1;
end
i=2;
while i<=6
    mye(i-1)=log(E(i-1)/E(i))/log(2);
    i=i+1;
end

function [p,e,x,y,u,i,et]=Solve(n,ep)
h=1/n;
t=1/n/1e-6;
x=0:h:1;
y=0:h:1;
a=zeros(n+1,n+1);
b=zeros(n+1,n+1);
ax=zeros(n+1,n+1);
by=zeros(n+1,n+1);
f=zeros(n+1,n+1);
fx=zeros(n+1,n+1);
fy=zeros(n+1,n+1);
u0=zeros(n+1,n+1);
u_=zeros(n+1,n+1);
e=1;
for i=1:n+1
    for j=1:n+1
        a(i,j)=aexact(x(i),y(j));
        b(i,j)=bexact(x(i),y(j));
        ax(i,j)=axexact(x(i),y(j));
        by(i,j)=byexact(x(i),y(j));
        f(i,j)=fexact(x(i),y(j));
        fx(i,j)=fxexact(x(i),y(j));
        fy(i,j)=fyexact(x(i),y(j));
        %u0(i,j)=uexact(x(i),y(j));
        p(i,j)=uexact(x(i),y(j));
    end
end
A1=zeros(1,n-1);
A2=zeros(1,n-1);
A3=zeros(1,n-1);
Y=zeros(1,n-1);
u=zeros(n+1,n+1);
for i=1:n
    for k=2:n
        for j=2:n
            aix=u0(j,k)+1;
            a1ix=(u0(j+1,k)-u0(j-1,k))/2/h;
            rix=h*aix/e;
            A1x=A1_2(rix,aix,a1ix,e);
            A2x=A2_2(rix,aix,a1ix,e);
            A3x=A3_2(rix,aix,a1ix,e);
            
            aiy=u0(j,k)+1;
            a1iy=(u0(j,k+1)-u0(j,k-1))/2/h;
            riy=h*aiy/e;
            A1y=A1_2(riy,aiy,a1iy,e);
            A2y=A2_2(riy,aiy,a1iy,e);
            A3y=A3_2(riy,aiy,a1iy,e);
            
            fi=Y0_2(rix,aix,a1ix,e)*f(j,k)/2+Y1_2(rix,aix,a1ix,e)*fx(j,k)/2+Y0_2(riy,aiy,a1iy,e)*f(j,k)/2+Y1_2(riy,aiy,a1iy,e)*fx(j,k)/2;
            
%             A1x*u0(j-1,k)+A2x*u0(j,k)+A3x*u0(j+1,k)+A1y*u0(j,k-1)+A2y*u0(j,k)+A3y*u0(j,k+1)
%             fi
            A1(j-1)=t*A1x;
            A2(j-1)=1+t*A2x;
            A3(j-1)=t*A3x;
            Y(j-1)=u0(j,k)-t*(A1y*u0(j,k-1)+A2y*u0(j,k)+A3y*u0(j,k+1))+t*fi;
        end
        
        u(2:n,k)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y);
    end
    for j=2:n
        for k=2:n
            aix=u0(j,k)+1;
            a1ix=(u0(j+1,k)-u0(j-1,k))/2/h;
            rix=h*aix/e;
            A1x=A1_2(rix,aix,a1ix,e);
            A2x=A2_2(rix,aix,a1ix,e);
            A3x=A3_2(rix,aix,a1ix,e);
            
            aiy=u0(j,k)+1;
            a1iy=(u0(j,k+1)-u0(j,k-1))/2/h;
            riy=h*aiy/e;
            A1y=A1_2(riy,aiy,a1iy,e);
            A2y=A2_2(riy,aiy,a1iy,e);
            A3y=A3_2(riy,aiy,a1iy,e);
            
            fi=Y0_2(rix,aix,a1ix,e)*f(j,k)/2+Y1_2(rix,aix,a1ix,e)*fx(j,k)/2+Y0_2(riy,aiy,a1iy,e)*f(j,k)/2+Y1_2(riy,aiy,a1iy,e)*fx(j,k)/2;
            
            A1(k-1)=t*A1y;
            A2(k-1)=1+t*A2y;
            A3(k-1)=t*A3y;
            Y(k-1)=u(j,k)-t*(A1x*u(j-1,k)+A2x*u(j,k)+A3x*u(j+1,k))+t*fi;
        end
        u0(j,2:n)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y);
    end
    et(i)=mynorm(u0,p,h);
    if(et(i)<ep)
        break;
    end
    u_=u0;
end
u=u0;
e=abs(u-p);
end

function [r]=aexact(x,y)
r=x+y+1;
%r=1;
end

function [r]=axexact(x,y)
r=1;
end

function [r]=bexact(x,y)
r=x+y+1;
%r=1;
end

function [r]=byexact(x,y)
r=1;
end

function [r]=fexact(x,y)
r=(x*y*(x - 1) + x*(x - 1)*(y - 1))*(x*y*(x - 1)*(y - 1) + 1) - 2*y*(y - 1) - 2*x*(x - 1) + (x*y*(y - 1) + y*(x - 1)*(y - 1))*(x*y*(x - 1)*(y - 1) + 1);
%r=2*x^3*y - x^3 + 4*x^2*y^2 - 3*x^2*y - 2*x^2 + 2*x*y^3 - 3*x*y^2 - 2*x*y + 3*x - y^3 - 2*y^2 + 3*y;
%r=2*x^2*y - 3*x^2 + 2*x*y^2 - 4*x*y + 3*x - 3*y^2 + 3*y;
end

function [r]=fxexact(x,y)
r=(x*y*(y - 1) + y*(x - 1)*(y - 1))^2 - 4*x + (x*y*(x - 1) + x*(x - 1)*(y - 1))*(x*y*(y - 1) + y*(x - 1)*(y - 1)) + (x*y*(x - 1)*(y - 1) + 1)*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)) + 2*y*(x*y*(x - 1)*(y - 1) + 1)*(y - 1) + 2;
end

function [r]=fyexact(x,y)
r=(x*y*(x - 1) + x*(x - 1)*(y - 1))^2 - 4*y + (x*y*(x - 1) + x*(x - 1)*(y - 1))*(x*y*(y - 1) + y*(x - 1)*(y - 1)) + (x*y*(x - 1)*(y - 1) + 1)*((x - 1)*(y - 1) + x*y + x*(y - 1) + y*(x - 1)) + 2*x*(x*y*(x - 1)*(y - 1) + 1)*(x - 1) + 2;
end

function [r]=uexact(x,y)
r=x*(1-x)*y*(1-y);
end

function u=Solve2_(e,a,f,n,u0)
h=1/n;
df=diff(f);
a1=diff(a);
for i=1:n-1
    ai=subs(a,i/n);
    a1i=subs(a1,i/n);
    ri=h*ai/e;
    fi=subs(f,i/n);
    dfi=subs(df,i/n);
    A1(i)=A1_2(ri,ai,a1i,e);
    A2(i)=A2_2(ri,ai,a1i,e);
    A3(i)=A3_2(ri,ai,a1i,e);
    Y(i)=Y0_2(ri,ai,a1i,e)*fi+Y1_2(ri,ai,a1i,e)*dfi;
end
Y(1)=Y(1)-A1(1)*u0(1);
A1(1)=[];
Y(n-1)=Y(n-1)-A3(n-1)*u0(2);
A3(n-1)=[];
u=[u0(1),double(zhuiganfa(A1,A2,A3,Y)),u0(2)];
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

function res=P_h_2_i(ri,ai,a1i,e)
res=a1i*e^2/2/ai^3*((ri^2-2*ri+2)-2*exp(-ri));
end

function res=P_p_2_i(ri,ai,a1i,e)
res=a1i*e^2/2/ai^3*(exp(-2*ri)*(ri^2+2*ri+2)-2*exp(-ri));
end

function res=Q_h_2_0i(ri,ai,a1i,e)
res=-a1i*e^2/2/ai^4*((1-exp(-ri)-ri*exp(-ri))*(ri^2-2*ri)+ri^3*exp(-ri));
end

function res=Q_p_2_0i(ri,ai,a1i,e)
res=-a1i*e^2/2/ai^4*((exp(-2*ri)-exp(-ri)+ri*exp(-ri))*(ri^2+2*ri)-ri^3*exp(-ri));
end

function res=Q_h_2_1i(ri,ai,e)
res=-e^2/ai^3*(1-exp(-ri)-ri*exp(-ri)-ri^2/2*exp(-ri));
end

function res=Q_p_2_1i(ri,ai,e)
res=-e^2/ai^3*(exp(-2*ri)-exp(-ri)+ri*exp(-ri)-ri^2/2*exp(-ri));
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

function res=A1_2(ri,ai,a1i,e)
res=-exp(-ri)*(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e));
end

function res=A2_2(ri,ai,a1i,e)
res=exp(-ri)*(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))-exp(-ri)*(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e));
end

function res=A3_2(ri,ai,a1i,e)
res=exp(-ri)*(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e));
end

function res=Y0_2(ri,ai,a1i,e)
res=(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e))*(Q_h_1_0i(ri,ai,e)+Q_h_2_0i(ri,ai,a1i,e))-(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))*(Q_p_1_0i(ri,ai,e)+Q_p_2_0i(ri,ai,a1i,e));
end

function res=Y1_2(ri,ai,a1i,e)
res=(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e))*Q_h_2_1i(ri,ai,e)-(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))*Q_p_2_1i(ri,ai,e);
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