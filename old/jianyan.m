clear
syms x t
p=2-x^2;
q=(x+1)*(t+1);
g=10*t^2*exp(-t)*x*(1-x);
solution(32,40,2^(-20),p,q,g);
function u=solution(Mh,Mt,eps,P,Q,G)
h=1/Mh;
global dt
dt=2/Mt;
x=0:h:1;
t=0:dt:2;
u=zeros(Mh+1,Mt+1);
for j=2:Mt
    for i=1:Mh
        g(i)=subs(subs(G,x(i)),t(i));
        if j*dt<=1 && j==1
            F(i)=g(i);
        elseif j*dt<=1 && j~=1
            F(i)=g(i)+1/dt*u(i,j-1);
        elseif j*dt>1 
            F(i)=u(i,j-1/dt)+g(i)+1/dt*u(i,j-1);
        end
        p(i)=subs(P,x(i));
        r(i)=h*p(i)/eps;
        q(i)=subs(subs(Q,x(i)),t(i));
        k1(i)=K_i_1(eps,p(i),r(i));
        k(i)=K_i_j(eps,p(i),q(i),r(i));
        k2(i)=K_i_2(eps,p(i),r(i));
        kw(i)=K_i_w(eps,p(i),r(i));
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

function res=A_i(eps,pi,ri)
res=eps/pi*(exp(ri)-1);
end
    
function res=A_i_w(eps,pi,ri)
res=eps/pi*(exp(-ri)-1);
end

function res=B_i_j(eps,pi,qij,ri)
global dt
res=(qij+1/dt)*eps/(pi^2)*(exp(ri-1-ri));
end

function res=B_i_j_w(eps,pi,qij,ri)
global dt
res=(qij+1/dt)*eps/(pi^2)*(exp(-ri-1+ri));
end

function res=C_i(eps,pi,ri)
res=eps/(pi^2)*(exp(ri)-1-ri);
end

function res=C_i_w(eps,pi,ri)
res=eps/(pi^2)*(exp(-ri)-1+ri);
end

function res=K_i_1(eps,pi,ri)
res=-A_i(eps,pi,ri);
end

function res=K_i_j(eps,pi,qij,ri)
res=A_i(eps,pi,ri)-A_i_w(eps,pi,ri)+A_i(eps,pi,ri)*B_i_j_w(eps,pi,qij,ri)-A_i_w(eps,pi,ri)*B_i_j(eps,pi,qij,ri);
end

function res=K_i_2(eps,pi,ri)
res=A_i_w(eps,pi,ri);
end

function res=K_i_w(eps,pi,ri)
res=A_i(eps,pi,ri)*C_i_w(eps,pi,ri)-A_i_w(eps,pi,ri)*C_i(eps,pi,ri);
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
