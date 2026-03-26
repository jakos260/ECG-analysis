% surflapl.m
% function SL=surflapl(VER,ITRI,mode)
% approximate the surface LAPLACIAN operator of a function defined at the
% nodes of a triangulated surface
% ITRI specifies three vertices for each triangle
% VER  are the vertex coordinates
% mode=0 surflapl as such
% mode=1 SL is multiplied by square root of the nodal surface area
%        and, when applied to a function tau with dimension time over a surface,
%        SL(tau) then has the dimension s/m

% a van oosterom
% method: see Huiskamp; J Comput Phys

function SL=surflapl(VER,ITRI,mode)
dim=size(VER);
nver=dim(1);
dim=size(ITRI);
ntri=dim(1);
[nbver,BVER,nbtri,BTRI]=voisin(nver,ITRI);

if mode==1,
   % compute, for all vertices, the area of the triangles around it
   area=nodearea(VER,ITRI);
end

% MAIN LOOP
SL=zeros(nver,nver);
dm=zeros(nver,1);
for i=1:nver,
   % compute row(i) of the SL matrix using angular weights
   BVER(i,:);
   dm(i)=0;
   for j=1:nbver(i),
       D(i,j)=norm(VER(BVER(i,j),:)-VER(i,:));
       dm(i)=dm(i)+D(i,j);
   end
   dm(i)=dm(i)/nbver(i);
   
   wtot=0;
   wgt=zeros(nbver(i),1);
   for j=1:nbver(i),
       jm=icyc(j-1,nbver(i));
       jp=icyc(j+1,nbver(i));
       r0=VER(BVER(i,j ),:)-VER(i,:);
       rm=VER(BVER(i,jm),:)-VER(i,:);
       rp=VER(BVER(i,jp),:)-VER(i,:);
       pim=rm*r0';
       pip=rp*r0';
       argm=pim/(D(i,jm)*D(i,j));
       amin=acos(argm);
       argp=pip/(D(i,j)*D(i,jp));
       apls=acos(argp);
       wgt(j)=(1-cos(amin))/sin(amin)+(1-cos(apls))/sin(apls);
       wtot=wtot+wgt(j);
	   wgt(j)=wgt(j)/D(i,j);
    end
    fac=4/(wtot*dm(i));
    if mode==1 fac=fac*sqrt(area(i)); end
    wgt=fac*wgt;
    for j=1:nbver(i),
        SL(i,BVER(i,j))=wgt(j);
    end
    SL(i,i)=-sum(wgt);
end 
