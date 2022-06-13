clear
tic
[p,e,x,y,u,k,et]=Solve(256,1e-8);
toc
surf(x,y,u);
mynorm(e,0,1/256)

n=4;
i=1;
while n<=256
    [p,e,x,y,u,k(i),et]=Solve(n,1e-8);
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
h=1/n;
%t=3000*h;
t=62.7;
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
        %myf(i,j)=fexact(x(i),y(j))+e*uxxexact(x(i),y(j))-aexact(x(i),y(j))*uxexact(x(i),y(j));
        myf(i,j)=fexact(x(i),y(j))+e*uyyexact(x(i),y(j))-bexact(x(i),y(j))*uyexact(x(i),y(j));
    end
end
Y=zeros(1,n-1);
u=zeros(n+1,n+1);

A1x_=zeros(n-1,n-1);
A2x_=zeros(n-1,n-1);
A3x_=zeros(n-1,n-1);
f_x=zeros(n-1,n-1);
Y0x_=zeros(n-1,n-1);
A1y_=zeros(n-1,n-1);
A2y_=zeros(n-1,n-1);
A3y_=zeros(n-1,n-1);
Y0y_=zeros(n-1,n-1);
f_y=zeros(n-1,n-1);

fix_=zeros(n-1,n-1);
fiy_=zeros(n-1,n-1);
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
        Y0x_(j-1,k-1)=Y0_1(rix,aix,e);

        aiy=b(j,k);
        riy=h*aiy/e;
        A1y_(j-1,k-1)=A1_1(riy,aiy,e);
        A2y_(j-1,k-1)=A2_1(riy,aiy,e);
        A3y_(j-1,k-1)=A3_1(riy,aiy,e);
        Y0y_(j-1,k-1)=Y0_1(riy,aiy,e);
        

        A11_(j-1,k-1)=t*A1_1(rix,aix,e);
        A12_(j-1,k-1)=1+t*A2_1(rix,aix,e);
        A13_(j-1,k-1)=t*A3_1(rix,aix,e);
        A21_(j-1,k-1)=t*A1_1(riy,aiy,e);
        A22_(j-1,k-1)=1+t*A2_1(riy,aiy,e);
        A23_(j-1,k-1)=t*A3_1(riy,aiy,e);
        
%         f_x(j-1,k-1)=f(j,k)+e*(u0(j,k+1)-2*u0(j,k)+u0(j,k-1))/h/h-b(j,k)*(u0(j,k+1)-u0(j,k-1))/2/h;
%         f_y(j-1,k-1)=f(j,k)+e*(u0(j+1,k)-2*u0(j,k)+u0(j-1,k))/h/h-a(j,k)*(u0(j+1,k)-u0(j-1,k))/2/h;
%         
        myFi(j-1,k-1)=Y0_1(riy,aiy,e)*f_y(j-1,k-1);
    end
end
% for k=2:n
%     u(k,2:n)=zhuiganfa(A1y_(k-1,2:n-1),A2y_(k-1,:),A3y_(k-1,1:n-2),myFi(k-1,:));
% end
% mynorm(u,u0,h)
A1=zeros(1,n-1);
A2=zeros(1,n-1);
A3=zeros(1,n-1);
Y=zeros(1,n-1);
f_x=zeros(1,n-1);
f_y=zeros(1,n-1);
for i=1:100000
    for k=2:n
        for j=2:n
            f_x(j-1)=f(j,k)+e*(u0(j,k+1)-2*u0(j,k)+u0(j,k-1))/h/h-b(j,k)*(u0(j,k+1)-u0(j,k-1))/2/h;
            f_y(j-1)=f(j,k)+e*(u0(j+1,k)-2*u0(j,k)+u0(j-1,k))/h/h-a(j,k)*(u0(j+1,k)-u0(j-1,k))/2/h;
        end
        for j=2:n
            
%             -e*(u0(j,k+1)-2*u0(j,k)+u0(j,k-1))/h/h+b(j,k)*(u0(j,k+1)-u0(j,k-1))/2/h-e*(u0(j+1,k)-2*u0(j,k)+u0(j-1,k))/h/h+a(j,k)*(u0(j+1,k)-u0(j-1,k))/2/h
%             f(j,k)
            
%             A1x_(j-1,k-1)*u0(j-1,k)+A2x_(j-1,k-1)*u0(j,k)+A3x_(j-1,k-1)*u0(j+1,k)
%             (b(j,k)*(u0(j,k+1)-u0(j,k-1))/2/h)
%             (b(j,k)*uyexact(x(j),y(k)))
            
%             e*(u0(j,k+1)-2*u0(j,k)+u0(j,k-1))/h/h
%             e*uyyexact(x(j),y(k))
%             
%             b(j,k)*(u0(j,k+1)-u0(j,k-1))/2/h
%             bexact(x(j),y(k))*uyexact(x(j),y(k))
%             A1y_(j-1,k-1)*u0(j,k-1)+A2y_(j-1,k-1)*u0(j,k)+A3y_(j-1,k-1)*u0(j,k+1)
%             Y0y_(j-1,k-1)*(f(j,k)+e*(u0(j+1,k)-2*u0(j,k)+u0(j-1,k))/h/h-a(j,k)*(u0(j+1,k)-u0(j-1,k))/2/h)
            %Y0y_(j-1,k-1)*(myf(j,k))
            
            %Y(k-1,j-1)=Y0_1(h*a(j,k)/e,a(j,k),e)*(f(j,k)+e*(u0(j,k+1)-2*u0(j,k)+u0(j,k-1))/h/h-b(j,k)*(u0(j,k+1)-u0(j,k-1))/2/h);
            
            A1(j-1)=t*A1x_(j-1,k-1);
            A2(j-1)=1+t*A2x_(j-1,k-1);
            A3(j-1)=t*A3x_(j-1,k-1);
            
%             A1y_(j-1,k-1)*u0(j,k-1)+A2y_(j-1,k-1)*u0(j,k)+A3y_(j-1,k-1)*u0(j,k+1)-Y0y_(j-1,k-1)*f_y
%             
%             A1x_(j-1,k-1)*u0(j-1,k)+A2x_(j-1,k-1)*u0(j,k)+A3x_(j-1,k-1)*u0(j+1,k)-Y0y_(j-1,k-1)*f_y
            
            Y(j-1)=u0(j,k)-t*(A1y_(j-1,k-1)*u0(j,k-1)+A2y_(j-1,k-1)*u0(j,k)+A3y_(j-1,k-1)*u0(j,k+1)-Y0y_(j-1,k-1)*f_y(j-1)-Y0x_(j-1,k-1)*f_x(j-1));
        end
        
        u(2:n,k)=zhuiganfa(A1(2:n-1),A2(:)',A3(1:n-2),Y(:));
    end
    for k=2:n
        for j=2:n
            A1(j-1)=t*A1x_(j-1,k-1);
            A2(j-1)=1+t*A2x_(j-1,k-1);
            A3(j-1)=t*A3x_(j-1,k-1);
            
            f_x=f(j,k)+e*(u0(j,k+1)-2*u0(j,k)+u0(j,k-1))/h/h-b(j,k)*(u0(j,k+1)-u0(j,k-1))/2/h;
            f_y=f(j,k)+e*(u0(j+1,k)-2*u0(j,k)+u0(j-1,k))/h/h-a(j,k)*(u0(j+1,k)-u0(j-1,k))/2/h;
            Y(j-1)=u0(j,k)-t*(A1y_(j-1,k-1)*u0(j,k-1)+A2y_(j-1,k-1)*u0(j,k)+A3y_(j-1,k-1)*u0(j,k+1)-Y0y_(j-1,k-1)*f_y-Y0x_(j-1,k-1)*f_x);
%             A1y_(j-1,k-1)*u0(j,k-1)+A2y_(j-1,k-1)*u0(j,k)+A3y_(j-1,k-1)*u0(j,k+1)-Y0y_(j-1,k-1)*f_y
%             
%             A1x_(j-1,k-1)*u0(j-1,k)+A2x_(j-1,k-1)*u0(j,k)+A3x_(j-1,k-1)*u0(j+1,k)-Y0y_(j-1,k-1)*f_x
%             Y(j-1)=u0(j,k)-t*(A1y_(j-1,k-1)*u0(j,k-1)+A2y_(j-1,k-1)*u0(j,k)+A3y_(j-1,k-1)*u0(j,k+1)-Y0y_(j-1,k-1)*f_y-Y0x_(j-1,k-1)*f_x);
%              u(j-1,k)-u0(j-1,k)
%             A1(j-1)*u(j-1,k)+A2(j-1)*u(j,k)+A3(j-1)*u(j+1,k)
%             Y(j-1)
%             -t*(A1y_(j-1,k-1)*u0(j,k-1)+A2y_(j-1,k-1)*u0(j,k)+A3y_(j-1,k-1)*u0(j,k+1)-Y0y_(j-1,k-1)*f_y+A1x_(j-1,k-1)*u(j-1,k)+A2x_(j-1,k-1)*u(j,k)+A3x_(j-1,k-1)*u(j+1,k)-Y0x_(j-1,k-1)*f_x)
        end
    end
    %surf(x,y,u);
    %mynorm(u,u0,h)
    for j=2:n
        for k=2:n
            %Y(j-1,k-1)=Y0_1(h*b(j,k)/e,b(j,k),e)*(f(j,k)+e*(u(j+1,k)-2*u(j,k)+u(j-1,k))/h/h-a(j,k)*(u(j+1,k)-u(j-1,k))/2/h);
            
            A1(k-1)=t*A1y_(j-1,k-1);
            A2(k-1)=1+t*A2y_(j-1,k-1);
            A3(k-1)=t*A3y_(j-1,k-1);
            
            f_x=f(j,k)+e*(u(j,k+1)-2*u(j,k)+u(j,k-1))/h/h-b(j,k)*(u(j,k+1)-u(j,k-1))/2/h;
            f_y=f(j,k)+e*(u(j+1,k)-2*u(j,k)+u(j-1,k))/h/h-a(j,k)*(u(j+1,k)-u(j-1,k))/2/h;
            
            Y(k-1)=u0(j,k)-t*(A1x_(j-1,k-1)*u0(j-1,k)+A2x_(j-1,k-1)*u0(j,k)+A3x_(j-1,k-1)*u0(j+1,k)-Y0x_(j-1,k-1)*f_x-Y0y_(j-1,k-1)*f_y);
        end 
        u0(j,2:n)=zhuiganfa(A1(2:n-1),A2(:)',A3(1:n-2),Y(:));
    end
%     mynorm(u0,p,h)
%     surf(x,y,u0);
%     pause(0.1);
    et(i)=mynorm(u0,u_,h);
%     if(et(i)<ep)
%         break;
%     end
    if(i==4000)
%         surf(x,y,u0);
%         mynorm(u0,p,h)
        i=4000;
        break
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

function [r]=axexact(x,y)
r=exp(x+y);
%r=x+y+1;
%r=x+1;
end

function [r]=bexact(x,y)
r=exp(3*x+y);
%r=x+y+1;
%r=x+y+1;
end

function [r]=byexact(x,y)
r=exp(3*x+y);
%r=x+y+1;
%r=x+y+1;
end


function [r]=fexact(x,y)
%r=(2*x - 2)*(y - 1)^2*(x + y + 1) - 2*(y - 1)^2 - 2*(x - 1)^2 + (2*y - 2)*(x - 1)^2*(x + y + 1);
r=2*exp(x)*sin(pi*y)*(exp(2*y) - 1) - exp(x + y)*(sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1) + exp(x)*sin(pi*y)*(exp(2*y) - 1)*(x - 1)) - exp(3*x + y)*(2*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) + pi*cos(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1)) + 4*exp(2*y)*sin(pi*y)*(exp(x) - 1)*(x - 1) + exp(x)*sin(pi*y)*(exp(2*y) - 1)*(x - 1) - pi^2*sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1) + 4*pi*exp(2*y)*cos(pi*y)*(exp(x) - 1)*(x - 1);
%r=2*x^3*y - x^3 + 4*x^2*y^2 - 3*x^2*y - 2*x^2 + 2*x*y^3 - 3*x*y^2 - 2*x*y + 3*x - y^3 - 2*y^2 + 3*y;
%r=2*sin(pi*y) - (x*sin(pi*y) + sin(pi*y)*(x - 1))*(x + 1) - x*pi^2*sin(pi*y)*(x - 1) - x*pi*cos(pi*y)*(x - 1)*(x + y + 1);
%r=(x + 1)*(pi*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) + pi*exp(pi*(x + y))*sin(pi*x)*sin(pi*y)) + (pi*exp(pi*(x + y))*cos(pi*y)*sin(pi*x) + pi*exp(pi*(x + y))*sin(pi*x)*sin(pi*y))*(x + y + 1) - 2*pi^2*exp(pi*(x + y))*cos(pi*x)*sin(pi*y) - 2*pi^2*exp(pi*(x + y))*cos(pi*y)*sin(pi*x);
end

function [r]=uexact(x,y)
r=-sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1)*(x - 1);
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

function [r]=uxexact(x,y)
r=- sin(pi*y)*(exp(x) - 1)*(exp(2*y) - 1) - exp(x)*sin(pi*y)*(exp(2*y) - 1)*(x - 1);
end

function [r]=uxxexact(x,y)
r=- 2*exp(x)*sin(pi*y)*(exp(2*y) - 1) - exp(x)*sin(pi*y)*(exp(2*y) - 1)*(x - 1);
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

function res=A1_2(ri,ai,a1i,e)
res=-exp(-ri)*(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e));
%res=-exp(-ri)*(P_h_1_i(ri,ai,e));
end

function res=A2_2(ri,ai,a1i,e)
res=exp(-ri)*(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))-exp(-ri)*(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e));
%res=exp(-ri)*(P_h_1_i(ri,ai,e))-exp(-ri)*(P_p_1_i(ri,ai,e));
end

function res=A3_2(ri,ai,a1i,e)
res=exp(-ri)*(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e));
%res=exp(-ri)*(P_p_1_i(ri,ai,e));
end

function res=Y0_2(ri,ai,a1i,e)
res=(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e))*(Q_h_1_0i(ri,ai,e)+Q_h_2_0i(ri,ai,a1i,e))-(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))*(Q_p_1_0i(ri,ai,e)+Q_p_2_0i(ri,ai,a1i,e));
%res=(P_p_1_i(ri,ai,e))*(Q_h_1_0i(ri,ai,e))-(P_h_1_i(ri,ai,e))*(Q_p_1_0i(ri,ai,e));
end

function res=Y1_2(ri,ai,a1i,e)
res=(P_p_1_i(ri,ai,e)+P_p_2_i(ri,ai,a1i,e))*Q_h_2_1i(ri,ai,e)-(P_h_1_i(ri,ai,e)+P_h_2_i(ri,ai,a1i,e))*Q_p_2_1i(ri,ai,e);
%res=(P_p_1_i(ri,ai,e))*Q_h_2_1i(ri,ai,e)-(P_h_1_i(ri,ai,e))*Q_p_2_1i(ri,ai,e);
%res=0;
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