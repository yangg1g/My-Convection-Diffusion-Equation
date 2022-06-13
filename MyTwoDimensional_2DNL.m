clear
tic
[p1,p2,x,y,u,v,k,et]=Solve(128,1e-6);
toc
surf(x,y,u);
et(k)


function [p1,p2,x,y,u,v,i,et]=Solve(n,ep)
div1=0.5;
h=1/n;
t=1e+4;
x=0:h:1;
y=0:h:1;
a=zeros(n+1,n+1);
b=zeros(n+1,n+1);
fu=zeros(n+1,n+1);
fv=zeros(n+1,n+1);
u0=ones(n+1,n+1)*1e-5;
v0=ones(n+1,n+1)*1e-5;
u_=zeros(n+1,n+1);
v_=zeros(n+1,n+1);
e=1;
for i=1:n+1
    for j=1:n+1
        fu(i,j)=fuexact(x(i),y(j));
        fv(i,j)=fvexact(x(i),y(j));
        %u0(i,j)=uexact(x(i),y(j));
        %v0(i,j)=vexact(x(i),y(j));
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
A1=zeros(1,n-1);
A2=zeros(1,n-1);
A3=zeros(1,n-1);
Y1=zeros(1,n-1);
Y2=zeros(1,n-1);
u=u0;
v=v0;
for i=1:10000
    for k=2:n
        for j=2:n
            aix=u0(j,k);
            rix=h*aix/e;
            Yx=Y0_1(rix,aix,e);
            A1x=A1_1(rix,aix,e);
            A2x=A2_1(rix,aix,e);
            A3x=A3_1(rix,aix,e);
            
            aiy=v0(j,k);
            riy=h*aiy/e;
            Yy=Y0_1(riy,aiy,e);
            A1y=A1_1(riy,aiy,e);
            A2y=A2_1(riy,aiy,e);
            A3y=A3_1(riy,aiy,e);
           
            fiu=Yy*fu(j,k);
            fiv=Yy*fv(j,k);
            
            A1(j-1)=t*A1x-t*(Yy-Yx)*(1/h/h+aix/2/h);
            A2(j-1)=1+t*A2x+t*(Yy-Yx)*2/h/h;
            A3(j-1)=t*A3x-t*(Yy-Yx)*(1/h/h-aix/2/h);
            
            Y1(j-1)=u0(j,k)-t*(A1y*u0(j,k-1)+A2y*u0(j,k)+A3y*u0(j,k+1))+t*fiu;
            Y2(j-1)=v0(j,k)-t*(A1y*v0(j,k-1)+A2y*v0(j,k)+A3y*v0(j,k+1))+t*fiv;
        end
        
        u(2:n,k)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y1);
        v(2:n,k)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y2);
    end
    for j=2:n
        for k=2:n
            aix=u(j,k);
            rix=h*aix/e;
            Yx=Y0_1(rix,aix,e);
            A1x=A1_1(rix,aix,e);
            A2x=A2_1(rix,aix,e);
            A3x=A3_1(rix,aix,e);
            
            aiy=v(j,k);
            riy=h*aiy/e;
            Yy=Y0_1(riy,aiy,e);
            A1y=A1_1(riy,aiy,e);
            A2y=A2_1(riy,aiy,e);
            A3y=A3_1(riy,aiy,e);
            
            fiu=Yx*fu(j,k);
            fiv=Yx*fv(j,k);
%            
            A1(k-1)=t*A1y-t*(Yx-Yy)*(1/h/h+aiy/2/h);
            A2(k-1)=1+t*A2y+t*(Yx-Yy)*2/h/h;
            A3(k-1)=t*A3y-t*(Yx-Yy)*(1/h/h-aiy/2/h);
            
%             A1(k-1)=t*A1y;
%             A2(k-1)=1+t*A2y;
%             A3(k-1)=t*A3y;

            Y1(k-1)=u(j,k)-t*(A1x*u(j-1,k)+A2x*u(j,k)+A3x*u(j+1,k))+t*fiu;
            Y2(k-1)=v(j,k)-t*(A1x*v(j-1,k)+A2x*v(j,k)+A3x*v(j+1,k))+t*fiv;
        end
        u0(j,2:n)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y1);
        v0(j,2:n)=zhuiganfa(A1(2:n-1),A2,A3(1:n-2),Y2);
    end
    et(i)=max(mynorm(u0,u_,h),mynorm(v0,v_,h));
    if(et(i)<ep)
        break;
    end
    u_=u0;
    v_=v0;
end
u=u0;
v=v0;
e=max(mynorm(u0,p1,h),mynorm(v0,p2,h));
end

function [r]=aexact(x,y)
%r=x+y+1;
r=x+1;
end

function [r]=bexact(x,y)
%r=x+y+1;
r=x+2*y+1;
end

function [r]=fuexact(x,y)
r=10*pi^2*sin(3*pi*x)*sin(pi*y) + 3*pi*cos(3*pi*x)*sin(3*pi*x)*sin(pi*y)^2 - x*pi*cos(pi*y)*sin(3*pi*x)*sin(pi*y)*(x - 1);
end

function [r]=fvexact(x,y)
r=2*sin(pi*y) - sin(3*pi*x)*sin(pi*y)*(x*sin(pi*y) + sin(pi*y)*(x - 1)) - x*pi^2*sin(pi*y)*(x - 1) + x^2*pi*cos(pi*y)*sin(pi*y)*(x - 1)^2;
end

function [r]=uexact(x,y)
r=sin(3*pi*x)*sin(pi*y);
end

function [r]=vexact(x,y)
r=-x*sin(pi*y)*(x - 1);
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

function [res]=mynorm(A,B,h)
T=(A-B).*(A-B);
res=sqrt(sum(sum(T)))*h;
end