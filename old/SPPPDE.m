clear
syms x t

e=2^-10;
p=2-x^2;
q=(x+1)*(t+1);
r=-1;
g=10*t^2*exp(-t)*x*(1-x);

% e=1;
% p=1;
% q=0;
% r=0;
% g=2*t - t*x - t*(x - 1) - x*(x - 1);

[x,t,u]=solution(32,40,e,p,q,r,g);
surf(x,t,u');
function [x,t,u]=solution(Mh,Mt,eps,P,Q,R,G)
h=1/Mh;
dt=2/Mt;
x=0:h:1;
t=0:dt:2;
u=zeros(Mh+1,Mt+1);
u0=zeros(Mh+1,Mt+1);
for j=1:Mt+1
    for i=1:Mh+1
        u0(i,j)=-x(i)*(x(i)-1)*t(j);
    end
end
for j=2:Mt+1
    k__1=zeros(1,Mh-1);
    k_1=zeros(1,Mh-1);
    k0=zeros(1,Mh-1);
    k=zeros(1,Mh-1);
    F=zeros(1,Mh-1);
    for i=1:Mh-1
        gi=subs(subs(G,x(i+1)),t(j));
%         if (j-1)*dt<=1 && j~=1
%             F(i)=gi+1/dt*u(i+1,j-1);
%         elseif (j-1)*dt>1 
%             F(i)=-R*u(i+1,j-1)+gi+1/dt*u(i+1,j-1);
%         end
        F(i)=-R*u(i+1,j-1)+gi+1/dt*u(i+1,j-1);
        pi=subs(P,x(i+1));
        vi=subs(subs(Q,x(i+1)),t(j))+1/dt;
        ai=h*pi/eps;
        k__1(i)=K_i__1(eps,pi,ai);
        k0(i)=K_i_0(eps,pi,vi,ai);
        k_1(i)=K_i_1(eps,pi,ai);
        k(i)=K_i(eps,pi,ai);
            
%         a1=exp(-ai)*u0(i+1,j)
%         a2=(exp(-ai)+C_i_j(eps,vi,pi,ai))*u0(i,j)+D_i(eps,pi,ai)*(u0(i+1,j)-u0(i-1,j))/2/h-E_i(eps,pi,ai)*F(i)
% 
%         b1=exp(-ai)*u0(i-1,j)
%         b2=(exp(-ai)+C_i_j_w(eps,vi,pi,ai))*u0(i,j)+D_i_w(eps,pi,ai)*(u0(i+1,j)-u0(i-1,j))/2/h-E_i_w(eps,pi,ai)*F(i)
% 
%         k(i)*F(i)
%         k__1(i)*u0(i-1,j)-k0(i)*u0(i,j)+k_1(i)*u0(i+1,j)
    end
    B=k.*F;
    B(1)=B(1)-k__1(1)*u(1,j);
    B(Mh-1)=B(Mh-1)-k_1(Mh-1)*u(Mh+1,j);
    k__1(1)=[];
    k_1(Mh-1)=[];
    u(2:Mh,j)=ChaseMethod(k__1,-k0,k_1,B)';
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

function res=C_i_j(eps,vi,pi,ai)
res=eps*vi/(pi^2)*(1-exp(-ai)-ai*exp(-ai));
end

function res=C_i_j_w(eps,vi,pi,ai)
res=eps*vi/(pi^2)*(exp(-2*ai)-exp(-ai)+ai*exp(-ai));
end

function res=D_i(eps,pi,ai)
res=eps/pi*(1-exp(-ai));
end
    
function res=D_i_w(eps,pi,ai)
res=eps/pi*(exp(-2*ai)-exp(-ai));
end

function res=E_i(eps,pi,ai)
res=eps/(pi^2)*(1-exp(-ai)-ai*exp(-ai));
end

function res=E_i_w(eps,pi,ai)
res=eps/(pi^2)*(exp(-2*ai)-exp(-ai)+ai*exp(-ai));
end

function res=K_i__1(eps,pi,ai)
res=-exp(-ai)*D_i(eps,pi,ai);
end

function res=K_i_0(eps,pi,vi,ai)
res=exp(-ai)*D_i_w(eps,pi,ai)-exp(-ai)*D_i(eps,pi,ai)+C_i_j(eps,vi,pi,ai)*D_i_w(eps,pi,ai)-C_i_j_w(eps,vi,pi,ai)*D_i(eps,pi,ai);
end

function res=K_i_1(eps,pi,ai)
res=exp(-ai)*D_i_w(eps,pi,ai);
end

function res=K_i(eps,pi,ai)
res=D_i(eps,pi,ai)*E_i_w(eps,pi,ai)-D_i_w(eps,pi,ai)*E_i(eps,pi,ai);
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