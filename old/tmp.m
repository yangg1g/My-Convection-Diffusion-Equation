clear
clc
n=8;
for i=1:n+1
    for j=1:n+1
        t="syms u_"+i+j+";";
        eval(t);
        t="u0("+i+","+j+")="+"u_"+i+j+";";
        eval(t);
        t="u("+i+","+j+")="+"u_"+i+j+";";
        eval(t);
    end
end
ij_ht=java.util.Hashtable;
ij=[1,1];
index=1;
while(ij(2)<n)
    ij_ht.put(ij(1)+"_"+ij(2),index);
    i=ij(1)+1;
    j=ij(2)+1;
    u(i,j)=(u0(i,j+1)+u(i,j-1)+u0(i+1,j)+u(i-1,j))/4;
    ij=nextij(ij,n-1);
    index=index+1;
end
A=zeros(49,49);
for i=1:n-1
    for j=1:n-1
        index1=ij_ht.get(i+"_"+j);
        for ii=1:n-1
            for jj=1:n-1
                index2=ij_ht.get(ii+"_"+jj);
                c=coeffs(u(i+1,j+1),u0(ii+1,jj+1));
                s=size(c);
                if(s(2)>1)
%                     if(index1==49&&index2==2)
%                         c
%                     end
                    A(index1,index2)=c(2);
                end
            end
        end
    end
end

% ij1=[1,1];
% index1=1;
% while(ij1(2)<n)
%     ij2=[1,1];
%     index2=1;
%     while(ij2(2)<n)
%         c=coeffs(u(ij1(1)+1,ij1(2)+1),u0(ij2(1)+1,ij2(2)+1));
%         s=size(c);
%         if(s(2)>1)
%             A(index1,index2)=c(2);
%         end
%         ij2=nextij(ij2,n-1);
%         index2=index2+1;
%     end
%     ij1=nextij(ij1,n-1);
%     index1=index1+1;
% end

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