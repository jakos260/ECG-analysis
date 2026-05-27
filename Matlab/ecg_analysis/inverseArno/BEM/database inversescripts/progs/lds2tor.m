% function T=lds2tor(filenameleadspecs,filenamesurface)
% computes the transfer from a scalar function defined at 
% nl positions (leads) on triangeles of a triangulated surface
% to all its vertices;lead positions specified by lambda and mu value of the 
% triangle that carries the lead

% 2003-04-18
% METHOD:
% first:
% include all lead positions as nodes in the geometry
% next:
% compute transfer from these nodes to all others
% by  minimizing the Surface Laplacian: L
% see: J COMPUT PHYS, 80/2, 331-343, 1989
% T=-(L1'*L1)^(-1)*L2 ; with
% L1 the columns of L related to the input-vertices
% L2 the columns of L related to the other vertices

function T=lds2tor(leadspecs,surface)
LEADS=loadmat(leadspecs);
%LEADS=sortrows(LEADSPEC,4);
% e.g. gubsme.lst
dim=size(LEADS);
[nl jdum]=size(LEADS)
[VER, ITRI]=loadtri(surface);
dim=size(VER);
nn=dim(1)
index=zeros(1,nn);
n2=0;
% inspect lead positions;
% for those close to a triangle vertex: accept vertex as lead position

% for the remaining ones close to an edge: create a new vertex; 
% adapt the local triangulation, replacing two triangles bordering the edge
% by four triangles; 

% for the remaining leads: create a new vertex; 
% split the triangle into three parts;

%NBNBNB
% it is assumed that a triangle carries at most one lead
% it is assumed that two triangles sharing an edge
% carry at most one lead

delta=0.1
for l=1:nl,
  la=LEADS(l,2);
  la=max(la,0);
  la=min(la,1);
  mu=LEADS(l,2);
  mu=max(mu,0);
  mu=min(mu,1);
  nu=1-la-mu;
  nu=max(nu,0);
  nu=min(nu,1);
  if la < delta & mu < delta,
     n2=n2+1;
     col2(n2)=ITRI(l,1);
  end
  if la > 1- delta & mu < delta,
     n2=n2+1;
     col2(n2)=ITRI(l,2);
  end
  if la < delta & mu >  1- delta,
     n2=n2+1;
     col2(n2)=ITRI(l,3);
  end
  n1=nn-n2;
  [nn n1 n2]
  L=surflapl(VER,ITRI,0);
  % create submatrices
  L1=zeros(nn,n1);
  L2=zeros(nn,n2);
  % form L1 and L2
  k=0;
  l=0;
  for j=1:nn,
    if index(j)==0,
       k=k+1;
       L1(:,k)=L(:,j);
    else
       l=l+1;
       L2(:,l)=L(:,j);
    end
  end
  %compute the interpolating matrix
  INT=-pinv(L1'*L1)*L1'*L2;
  size(INT);
  k=0;
  l=0;
  T=zeros(nn,n2);
  for i=1:nn;
     if index(i)==1,
        k=k+1;
        T(i,k)=1;
     else
        l=l+1;
        for m=1:n2,
           T(i,m)=INT(l,m);
        end
     end
  end
end