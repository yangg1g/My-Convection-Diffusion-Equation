clc
clear

m=5;
A=zeros(m,m);
nowij=[1,1];
nowt=1;
while(nowij(2)<m+1)
    A(nowij(1),nowij(2))=nowt;
    nowt=nowt+1;
    nowij=nextij(nowij,m);
end
% h=1/m;
% u=zeros(m+1,m+1);
% A=zeros((m-1)*(m-1),(m-1)*(m-1));
% for i=1:m-1
%     for j=1:m-1
%         t=ij2index([i,j],m-1);
%         A(t,t)=4;
%         if(i>1)
%             A(t,ij2index([i-1,j],m-1))=-1;
%         end
%         if(j>1)
%             A(t,ij2index([i,j-1],m-1))=-1;
%         end
%         if(i<m-1)
%             A(t,ij2index([i+1,j],m-1))=-1;
%         end
%         if(j<m-1)
%             A(t,ij2index([i,j+1],m-1))=-1;
%         end
%     end
% end

function res=ij2index(ij,m)
% n=ij(1)+ij(2)-1;
% if(n<=m)
%     if(mod(n,2)==0)
%         res=n*(n-1)/2+ij(2);
%     else
%         res=n*(n-1)/2+ij(2);
%     end
% else
%     res=m*m-ij2index([m-ij(1)+1,m-ij(2)+1],m)+1;
% end
res=(ij(2)-1)*m+ij(1);
end

function res=index2ij(index)
n=1;
while(index-(n+1)*n/2>0)
    n=n+1;
end
index=index-(n-1)*n/2;
res(1)=n-index+1;
res(2)=index;
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