clear
% myt=0;%13920
% maxe=1;
% for t=1:1:13920
%     [p,e,x,y,u,k,et]=Solve(128,1e-4,t);
%     nowe=mynorm(e,0,1/128);
%     if(nowe<maxe)
%         maxe=nowe;
%         myt=t;
%     end
% end
% myt=0;
% maxe=1;
% i=1;
% n=8;
% for t=0.1:0.001:0.5
%     [p,e,x,y,u,k,et]=Solve(n,1e-4,t);
%     nowe(i)=mynorm(e,0,1/n);
%     if(nowe(i)<maxe)
%         maxe=nowe(i);
%         myt=t;
%     end
%     i=i+1;
% end
% a=[ 2 3 4];
% b=[4 4 6 7];
% c=[ 8 8 4];
% d = [9 7 2 5];
% AA=[4 8 0 0;
%     2 4 8 0;
%     0 3 6 4;
%     0 0 4 7;]
% x=crout(b,c,a,d)
tic
[p,e,x,y,u,k,et]=Solve(128,1e-9);
toc
surf(x,y,u);
% norm(e,inf);
mynorm(e,0,1/128)
n=4;
i=1;
while n<=256
    [p,e,x,y,u,k(i),et]=Solve(n,1e-6);
    E(i)=mynorm(e,0,1/n);
    n=n*2;
    i=i+1;
end
i=2;
while i<=7
    mye(i-1)=log(E(i-1)/E(i))/log(2);
    i=i+1;
end

function [p,e,x,y,u,i,et]=Solve(n,ep)
div1=0.184;
%div1=0.267;
h=1/n;
t=13920;
x=0:h:1;
y=0:h:1;
a=zeros(n+1,n+1);
b=zeros(n+1,n+1);
f=zeros(n+1,n+1);
u0=zeros(n+1,n+1);
u_=zeros(n+1,n+1);
p=zeros(n+1,n+1);
myf=zeros(n+1,n+1);
e=1;
for i=1:n+1
    for j=1:n+1
        a(i,j)=aexact(x(i),y(j));
        b(i,j)=bexact(x(i),y(j));
        f(i,j)=fexact(x(i),y(j));
        %u0(i,j)=uexact(x(i),y(j));
        p(i,j)=uexact(x(i),y(j));
        myf(i,j)=fexact(x(i),y(j))+e*uyyexact(x(i),y(j))-bexact(x(i),y(j))*uyexact(x(i),y(j));
    end
end
Y=zeros(1,n-1);
u=zeros(n+1,n+1);

A1x_=zeros(n-1,n-1);
A2x_=zeros(n-1,n-1);
A3x_=zeros(n-1,n-1);
A1y_=zeros(n-1,n-1);
A2y_=zeros(n-1,n-1);
A3y_=zeros(n-1,n-1);
fi_=zeros(n-1,n-1);
A11_=zeros(n-1,n-1);
A12_=zeros(n-1,n-1);
A13_=zeros(n-1,n-1);
A21_=zeros(n-1,n-1);
A22_=zeros(n-1,n-1);
A23_=zeros(n-1,n-1);

myFi=zeros(n-1,n-1);
for k=2:n
    for j=2:n
        aix=a(j,k);
        rix=h*aix/e;
        A1x_(j-1,k-1)=A1_1(rix,aix,e);
        A2x_(j-1,k-1)=A2_1(rix,aix,e);
        A3x_(j-1,k-1)=A3_1(rix,aix,e);

        aiy=b(j,k);
        riy=h*aiy/e;
        A1y_(j-1,k-1)=A1_1(riy,aiy,e);
        A2y_(j-1,k-1)=A2_1(riy,aiy,e);
        A3y_(j-1,k-1)=A3_1(riy,aiy,e);
        
        fi_(j-1,k-1)=f(j,k)*div1*Y0_1(rix,aix,e)+f(j,k)*(1-div1)*Y0_1(riy,aiy,e);
        
        A11_(j-1,k-1)=t*A1_1(rix,aix,e);
        A12_(j-1,k-1)=1+t*A2_1(rix,aix,e);
        A13_(j-1,k-1)=t*A3_1(rix,aix,e);
        A21_(j-1,k-1)=t*A1_1(riy,aiy,e);
        A22_(j-1,k-1)=1+t*A2_1(riy,aiy,e);
        A23_(j-1,k-1)=t*A3_1(riy,aiy,e);
        
        myFi(j-1,k-1)=t*Y0_1(rix,aix,e)*myf(j,k);
        
        Y0x_(j-1,k-1)=Y0_1(rix,aix,e);
        Y0y_(j-1,k-1)=Y0_1(riy,aiy,e);
    end
end
% for k=2:n
%     u(2:n,k)=zhuiganfa(A11_(2:n-1,k-1),A12_(:,k-1)',A13_(1:n-2,k-1),myFi(:,k-1));
% end
for i=1:100000
    for k=2:n
        A=zeros(n-1,n+1);
        for j=2:n
            A(j-1,j-1)=A11_(j-1,k-1);
            A(j-1,j)=A12_(j-1,k-1);
            A(j-1,j+1)=A13_(j-1,k-1);
%             f_x=f(j,k)+e*(u0(j,k+1)-2*u0(j,k)+u0(j,k-1))/h/h-b(j,k)*(u0(j,k+1)-u0(j,k-1))/2/h;
%             f_y=f(j,k)+e*(u0(j+1,k)-2*u0(j,k)+u0(j-1,k))/h/h-a(j,k)*(u0(j+1,k)-u0(j-1,k))/2/h;
%             A1y_(j-1,k-1)*u0(j,k-1)+A2y_(j-1,k-1)*u0(j,k)+A3y_(j-1,k-1)*u0(j,k+1)-Y0y_(j-1,k-1)*f_y
%             A1x_(j-1,k-1)*u0(j-1,k)+A2x_(j-1,k-1)*u0(j,k)+A3x_(j-1,k-1)*u0(j+1,k)-Y0y_(j-1,k-1)*f_y
            Y(j-1)=u0(j,k)-t*(A1y_(j-1,k-1)*u0(j,k-1)+A2y_(j-1,k-1)*u0(j,k)+A3y_(j-1,k-1)*u0(j,k+1))+t*fi_(j-1,k-1);
        end
        u(2:n,k)=zhuiganfa(A11_(2:n-1,k-1),A12_(:,k-1)',A13_(1:n-2,k-1),Y);
        %u(2:n,k)=linsolve(A(:,2:n),Y');
    end
    for k=2:n
        for j=2:n
            Y(j-1)=u0(j,k)-t*(A1y_(j-1,k-1)*u0(j,k-1)+A2y_(j-1,k-1)*u0(j,k)+A3y_(j-1,k-1)*u0(j,k+1))+t*fi_(j-1,k-1);
%             u(j-1,k)*A11_(j-1,k-1)+u(j,k)*A12_(j-1,k-1)+u(j+1,k)*A13_(j-1,k-1)
%             Y(j-1)
        end
    end
    %surf(x,y,u);
    A=zeros(n-1,n+1);
    for j=2:n
        for k=2:n
            A(k-1,k-1)=A11_(j-1,k-1);
            A(k-1,k)=A12_(j-1,k-1);
            A(k-1,k+1)=A13_(j-1,k-1);
            Y(k-1)=u(j,k)-t*(A1x_(j-1,k-1)*u(j-1,k)+A2x_(j-1,k-1)*u(j,k)+A3x_(j-1,k-1)*u(j+1,k))+t*fi_(j-1,k-1);
        end 
        u0(j,2:n)=zhuiganfa(A21_(j-1,2:n-1),A22_(j-1,:),A23_(j-1,1:n-2),Y);
        %u0(j,2:n)=linsolve(A(:,2:n),Y');
    end
    mynorm(u0,p,h);
    %surf(x,y,u0);
    et(i)=mynorm(u0,u_,h);
    if(et(i)<ep)
        mynorm(u0,p,h);
        break;
    end
    u_=u0;
end
u=u0;
e=abs(u-p);
end

function [r]=aexact(x,y)
r=exp(x+y);
%r=x+y+1;
%r=x+1;
end

function [r]=bexact(x,y)
r=exp(x+y);
%r=exp(3*x+y);
%r=x+y+1;
%r=x+y+1;
end

function [r]=fexact(x,y)
r=exp(x + y)*((exp(x) - 1)*(exp(y) - 1)*(x - 1) + exp(y)*(exp(x) - 1)*(x - 1)*(y - 1)) + exp(x + y)*((exp(x) - 1)*(exp(y) - 1)*(y - 1) + exp(x)*(exp(y) - 1)*(x - 1)*(y - 1)) - 2*exp(y)*(exp(x) - 1)*(x - 1) - 2*exp(x)*(exp(y) - 1)*(y - 1) - exp(x)*(exp(y) - 1)*(x - 1)*(y - 1) - exp(y)*(exp(x) - 1)*(x - 1)*(y - 1);
%r=2*exp(x)*sin(pi*y)*(exp(2*y) - 1) - exp(x + y)*(2*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) + pi*cos(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1)) - exp(x + y)*(sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1) + exp(x)*sin(pi*y)*(exp(2*y) - 1)*(x - 1)) + 4*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) + exp(x)*sin(pi*y)*(exp(2*y) - 1)*(x - 1) - pi^2*sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1) + 4*pi*exp(2*y)*cos(pi*y)*(exp(x) - 1)*(x - 1);
%r=2*exp(x)*sin(pi*y)*(exp(2*y) - 1) - exp(x + y)*(sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1) + exp(x)*sin(pi*y)*(exp(2*y) - 1)*(x - 1)) - exp(3*x + y)*(2*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) + pi*cos(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1)) + 4*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) + exp(x)*sin(pi*y)*(exp(2*y) - 1)*(x - 1) - pi^2*sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1) + 4*pi*exp(2*y)*cos(pi*y)*(exp(x) - 1)*(x - 1);
%r=2*x^3*y - x^3 + 4*x^2*y^2 - 3*x^2*y - 2*x^2 + 2*x*y^3 - 3*x*y^2 - 2*x*y + 3*x - y^3 - 2*y^2 + 3*y;
%r=2*sin(pi*y) - (x*sin(pi*y) + sin(pi*y)*(x - 1))*(x + 1) - x*pi^2*sin(pi*y)*(x - 1) - x*pi*cos(pi*y)*(x - 1)*(x + y + 1);
%r=(x + 1)*(pi*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) + pi*exp(pi*(x + y))*sin(pi*x)*sin(pi*y)) + (pi*exp(pi*(x + y))*cos(pi*y)*sin(pi*x) + pi*exp(pi*(x + y))*sin(pi*x)*sin(pi*y))*(x + y + 1) - 2*pi^2*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) - 2*pi^2*exp(pi*(x + y))*cos(pi*y)*sin(pi*x);
end

function [r]=uexact(x,y)
r=(exp(x) - 1)*(exp(y) - 1)*(x - 1)*(y - 1);
%r=(exp(x) - 1)*(exp(y) - 1)*(x - 1)*(y - 1);
%r=-sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1);
%r=x*(1-x)*y*(1-y);
%r=-x*sin(pi*y)*(x - 1);
%r=exp(pi*(x+y))*sin(x*pi)*sin(y*pi);
end

function [r]=uyexact(x,y)
r=- 2*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) - pi*cos(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1);
end

function [r]=uyyexact(x,y)
r=pi^2*sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1) - 4*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) - 4*pi*exp(2*y)*cos(pi*y)*(exp(x) - 1)*(x - 1);
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
l=zeros(1,n-1);
u=zeros(1,n);
y=zeros(1,n);
x=zeros(1,n);
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

function [x]=crout(a,c,d,b)%数组a存储三角矩阵A的主对角线元素，c、d存储主对角线上边下边带宽为1的元素
    n=length(a);
    n1=length(c);
    n2=length(d);
    %错误检查
    if n1~=n2%存储矩阵的数组维数错误
        error('MATLAB:Crout:不是三对角矩阵，参数数组中元素个数错误.');
    elseif n~=n1+1
        error('MATLAB:Crout:不是三对角矩阵，参数数组中元素个数错误.');
    end
   
    %初始化
    L=zeros(n);%生成n*n的全零矩阵
    U=zeros(n);
    p=1:n;
    q=1:n-1;
    x=1:n;
    y=1:n;
   
    %追赶法程序主体
    p(1)=a(1);
    for i=1:n-1
        q(i)=c(i)/p(i);
        p(i+1)=a(i+1)-d(i)*q(i);%d的下标改为1到n-1
    end
    %正解y
    y(1)=b(1)/p(1);%用x存储y
    for i=2:n
        y(i)=(b(i)-d(i-1)*y(i-1))/p(i);
    end
    %倒解x
    x(n)=y(n);
    for i=(n-1):-1:1
        x(i)=y(i)-q(i)*x(i+1);
    end
    %L,U矩阵
    for i=1:n
        L(i,i)=p(i);
        U(i,i)=1;
    end
    for i=1:n-1
        L(i+1,i)=d(i);
        U(i,i+1)=q(i);
    end
end
