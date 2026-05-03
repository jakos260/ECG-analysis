function T=distripol(VER,LD,nodes2)
%function T=intripol(VER,ITRI,nodes2)
% 30-05-10
% computes the transfer from n2 values defined at nodes2(j=1:n2)
% of a triangulated surface to all of its nn nodes
% nodes2 contains the n2 nodes at which the values are known, i.e.,
% the righthandside 

% method: minimizing the Surface Laplacian: L
% see: J COMPUT PHYS, 80/2, 331-343, 1989

% T=-(L1'*L1)^(-1)*L1'*L2 ; with
% L1 the columns of L related to the nodes whose values need to be determined
% L2 the columns of L related to the nodes for which values are known

nn=length(VER);
n2=length(nodes2);
n1=nn-n2;
nodes1=zeros(1,n1);
index=zeros(1,nn);
index(nodes2)=1;

% % L=surflapl(VER,ITRI,0);
LD(LD>0) = (LD(LD>0).^-1.5);
L = LD - diag(sum(LD));

% create submatrices
L1=zeros(nn,nn-n2);
L2=zeros(nn,n2);
% form L1 and L2
k=0;
l=0;

for j=1:nn,
  if index(j)==0,
    k=k+1;
    L1(:,k)=L(:,j);
    nodes1(k)=j;
  else
    l=l+1;
    L2(:,l)=L(:,nodes2(l));
  end
end

%compute the interpolating matrix
INT=-pinv(L1'*L1)*L1'*L2;

% form tranfer matrix T
T=zeros(nn,n2);
T(nodes2,:)=eye(n2);
T(nodes1,:)=INT;

 	