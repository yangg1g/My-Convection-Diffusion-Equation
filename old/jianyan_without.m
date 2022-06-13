clear
syms x t

% e=2^(-20);
% p=2-x^2;
% q=(x+1)*(t+1);
% r=-1;
% g=10*t^2*exp(-t)*x*(1-x);

e=1;
p=1;
q=0;
r=0;
g=-1+2*x;

[x,t,u]=solution(32,40,e,p,q,r,g);
surf(x,t,u');
function [x,t,u]=solution(Mh,Mt,eps,P,Q,R,G)
h=1/Mh;
global dt
dt=2/Mt;
x=0:h:1;
t=0:dt:2;
u=zeros(Mh+1,Mt+1);
for j=2:Mt+1
    for i=1:Mh
        g(i)=double(subs(subs(G,x(i+1)),t(i+1)));
        if (j-1)*dt<=1 && j==1
            F(i)=g(i);
        elseif (j-1)*dt<=1 && j~=1
            F(i)=g(i)+1/dt*u(i,j-1);
        elseif (j-1)*dt>1 
            F(i)=-R*u(i,j-1)+g(i)+1/dt*u(i,j-1);
        end
        p(i)=double(subs(P,x(i+1)));
        q(i)=double(subs(subs(Q,x(i+1)),t(i+1)));
        a(i)=h*p(i)/eps;
        k1(i)=K_i_1(eps,p(i),a(i));
        k(i)=K_i_j(eps,p(i),q(i),a(i));
        k2(i)=K_i_2(eps,p(i),a(i));
        kw(i)=K_i_w(eps,p(i),a(i));
    end
    Y=kw.*F;
    Y(1)=Y(1)-k1(1)*u(1,j);
    Y(Mh)=Y(Mh)-k2(Mh)*u(Mh+1,j);
    k1(1)=[];
    k2(Mh)=[];
    u(2:Mh,j)=chase(k,k2,k1,Y);
end
end

% function res=P(x)
% res=2-x^2;
% end
% 
% function res=Q(x,t)
% res=(x+1)*(t+1);
% end
% 
% function res=G(x,t)
% res=10*t^2*exp(-t)*x*(1-x);
% end

function res=A_i(eps,pi,ai)
res=eps/pi*(1-exp(-ai));
end
    
function res=A_i_w(eps,pi,ai)
res=eps/pi*(exp(-2*ai)-exp(-ai));
end

function res=B_i_j(eps,pi,qij,ai)
global dt
res=(qij+1/dt)*eps/(pi^2)*(1-exp(-ai)-ai*exp(-ai));
end

function res=B_i_j_w(eps,pi,qij,ai)
global dt
res=(qij+1/dt)*eps/(pi^2)*(exp(-2*ai)-exp(-ai)+ai*exp(-ai));
end

function res=C_i(eps,pi,ai)
res=eps/(pi^2)*(1-exp(-ai)-ai*exp(-ai));
end

function res=C_i_w(eps,pi,ai)
res=eps/(pi^2)*(exp(-2*ai)-exp(-ai)+ai*exp(-ai));
end

function res=K_i_1(eps,pi,ai)
res=-A_i(eps,pi,ai)*exp(-ai);
end

function res=K_i_j(eps,pi,qij,ai)
res=A_i(eps,pi,ai)*exp(-ai)-A_i_w(eps,pi,ai)*exp(-ai)+A_i(eps,pi,ai)*B_i_j_w(eps,pi,qij,ai)-A_i_w(eps,pi,ai)*B_i_j(eps,pi,qij,ai);
end

function res=K_i_2(eps,pi,ai)
res=A_i_w(eps,pi,ai)*exp(-ai);
end

function res=K_i_w(eps,pi,ai)
res=A_i(eps,pi,ai)*C_i_w(eps,pi,ai)-A_i_w(eps,pi,ai)*C_i(eps,pi,ai);
end

%追赶法求解
function [x]=chase(a,b,c,d)%a为主对角线，c为上对角线，d为下对角线，b为右侧“常数”值
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

function [x]=ChaseMethod(a,b,c,d)
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