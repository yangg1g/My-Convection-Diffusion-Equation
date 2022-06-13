clear

% [p,e,x,y,u,i,et]=Solve(128,1e-6,0.1,10000);
% surf(x,y,u);

n=256;
vl=0;
vr=2;
[p,el,x,y,u,i,et]=Solve(n,1e-6,vl,10);
[p,er,x,y,u,i,et]=Solve(n,1e-6,vr,10);
while(1)
    %vm=exp((log(vl)+log(vr))/2);
    vm=(vl+vr)/2;
    [p,em,x,y,u,i,et]=Solve(n,1e-6,vm,10);
    if(el<er)
        vr=vm;
        er=em;
    else
        vl=vm;
        el=em;
    end
    if(vr-vl<0.1)
        break;
    end
end
[p,e,x,y,u,i,et]=Solve(n,1e-6,vm,10000);
surf(x,y,u);
e
vm
i
n=4;
i=1;
while n<=128
    [p,e,x,y,u,k(i),et]=Solve(n,1e-6);
    E(i)=mynorm(e,0,1/n);
    n=n*2;
    i=i+1;
end
i=2;
while i<=5
    mye(i-1)=log(E(i-1)/E(i))/log(2);
    i=i+1;
end

function [p,e,x,y,u,i,et]=Solve(n,ep,t,k)
div1=0.5;
h=1/n;
%t=1e+5;
%t=2e+3;
x=0:h:1;
y=0:h:1;
a=zeros(n+1,n+1);
b=zeros(n+1,n+1);
f=zeros(n+1,n+1);
u0=ones(n+1,n+1)*1e-5+rand(n+1,n+1)*1e-5;
u_=zeros(n+1,n+1);
ux=zeros(n+1,n+1);
uxx=zeros(n+1,n+1);
uy=zeros(n+1,n+1);
uyy=zeros(n+1,n+1);
e=0.01;
for i=1:n+1
    for j=1:n+1
        a(i,j)=aexact(x(i),y(j));
        b(i,j)=bexact(x(i),y(j));
        f(i,j)=fexact(x(i),y(j));
        %u0(i,j)=uexact(x(i),y(j))+2*rand()*0.001;
        p(i,j)=uexact(x(i),y(j));
        ux(i,j)=uxexact(x(i),y(j));
        uxx(i,j)=uxxexact(x(i),y(j));
        uy(i,j)=uyexact(x(i),y(j));
        uyy(i,j)=uyyexact(x(i),y(j));
    end
end
for i=1:n+1
    u0(1,i)=p(1,i);
    u0(n+1,i)=p(n+1,i);
    u0(i,1)=p(i,1);
    u0(i,n+1)=p(i,n+1);
end
%u0=p;
% for i=2:n
%     for j=2:n
%         u0(i,j)=u0(i,1)*(n+1-j)/(n+1)+u0(i,n+1)*(j)/(n+1);
%     end
% end
u=u0;
A1=zeros(1,n-1);
A2=zeros(1,n-1);
A3=zeros(1,n-1);
Y=zeros(1,n-1);
for i=1:k
    for k=2:n
       for j=2:n
            aix=a(j,k);
            rix=h*aix/e;
            Yx=Y0_1(rix,aix,e);
            A1x=A1_1(rix,aix,e);
            A2x=A2_1(rix,aix,e);
            A3x=A3_1(rix,aix,e);

            aiy=b(j,k);
            riy=h*aiy/e;
            Yy=Y0_1(riy,aiy,e);
            A1y=A1_1(riy,aiy,e);
            A2y=A2_1(riy,aiy,e);
            A3y=A3_1(riy,aiy,e);

            A_x_1_y=A1x-(Yy-Yx)*(e/h/h+aix/2/h);
            A_x_y=A2x+(Yy-Yx)*2*e/h/h;
            A_x1_y=A3x-(Yy-Yx)*(e/h/h-aix/2/h);
            Y_x_y=-(A1y*u(j,k-1)+A2y*u(j,k)+A3y*u(j,k+1))+Yy*f(j,k);
        
            A_x_y_1=A1y-(Yx-Yy)*(e/h/h+aiy/2/h);
            A_x_y=A_x_y+A2y+(Yx-Yy)*2*e/h/h;
            A_x_y1=A3y-(Yx-Yy)*(e/h/h-aiy/2/h);

            Y_x_y=Y_x_y-(A1x*u(j-1,k)+A2x*u(j,k)+A3x*u(j+1,k))+Yx*f(j,k);
            
            u(j,k)=u(j,k)+t*(Y_x_y-A_x_1_y*u(j-1,k)-A_x1_y*u(j+1,k)-A_x_y_1*u(j,k-1)-A_x_y1*u(j,k+1)-A_x_y*u(j,k))/A_x_y;
       end
    end
    et(i)=mynorm(u0,u,h);
    if(et(i)<ep*t)
        break;
    end

    if(isnan(et(i)))
        break;
    end
    u0=u;
end
% u=u0;
e=mynorm(u,p,h);
end

function [r]=aexact(x,y)
%r=x+y+1;
r=x+1;
end

function [r]=bexact(x,y)
%r=x+y+1;
r=x+2*y+1;
end

function [r]=fexact(x,y)
r=(exp(x) + (x + 1)^101/1267650600228229401496703205376)*(exp(y) + (101*(y + 1)^100)/1267650600228229401496703205376)*(x + 2*y + 1) - ((exp(x) + (2525*(x + 1)^99)/316912650057057350374175801344)*(exp(y) + (y + 1)^101/1267650600228229401496703205376))/100 - ((exp(x) + (x + 1)^101/1267650600228229401496703205376)*(exp(y) + (2525*(y + 1)^99)/316912650057057350374175801344))/100 + (exp(x) + (101*(x + 1)^100)/1267650600228229401496703205376)*(exp(y) + (y + 1)^101/1267650600228229401496703205376)*(x + 1);
%r=10*pi^2*sin(3*pi*x)*sin(pi*y) + pi*cos(pi*y)*sin(3*pi*x)*(x + 2*y + 1) + 3*pi*cos(3*pi*x)*sin(pi*y)*(x + 1);
%r=10*pi^2*sin(3*pi*x)*sin(pi*y) + pi*cos(pi*y)*sin(3*pi*x)*(x + 2*y + 1) + 3*pi*cos(3*pi*x)*sin(pi*y)*(x + 1);
%r=2*pi^2*sin(pi*x)*sin(pi*y) + pi*cos(pi*y)*sin(pi*x)*(x + y + 1) + pi*cos(pi*x)*sin(pi*y)*(x + 1);
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
%r=sin(pi*x)*sin(pi*y);
%r=x*(1-x)*y*(1-y);
%r=-x*sin(pi*y)*(x - 1);
%r=exp(pi*(x+y))*sin(x*pi)*sin(y*pi);
%r=(exp(x)-1)*(exp(x)-exp(1))*(exp(y)-1)*(exp(y)-exp(1));
end

function [r]=uxexact(x,y)
r=pi*cos(pi*x)*sin(pi*y);
end

function [r]=uxxexact(x,y)
r=-pi^2*sin(pi*x)*sin(pi*y);
end

function [r]=uyexact(x,y)
r=pi*cos(pi*y)*sin(pi*x);
end

function [r]=uyyexact(x,y)
r=-pi^2*sin(pi*x)*sin(pi*y);
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