% rowforw.m; 
% function [rowb,jsing]=rowforw(VER,ITRI,obs)

% 2003-03-06; A. van Oosterom
% improved treatment of auto solid angle in case of concave parts of
% bounding surface

%  compute (part of a) row of Bmatrix:
%  distributed solid angles as seen by obs subtended
%  by triangles ITRI representing a single closed surface
%  VER are its vertices
%  called by forward.m

% functions called: dsa.m      (general distributed solid angle function)
%                   loopnode.m (identifies loop of direct vertex neighbours
%		                around node; orientation clockwise when viewed from
%		                outside) 
%                    solida.m   (solid angle function)

function [rowb,jsing]=rowforw(VER,ITRI,obs)
nver=length(VER);
ntri=length(ITRI);

jsing=[];
% identify singularity (if present)
jsing=find(VER(:,1)==obs(1)&VER(:,2)==obs(2)&VER(:,3)==obs(3));

OMEGA=zeros(nver,ntri);
rowb=zeros(1,nver);
[OMEGA,index]=dsa(VER,ITRI,obs,.1);
for j=1:ntri,
    ij=ITRI(j,:);
    rowb(ij)=rowb(ij) + OMEGA(:,j)';
end
 
if isempty(jsing)==1, jsing=0; end

if jsing~=0,% treat singularity (determine auto  solid angle)
   llist=find(ITRI(:,1)==jsing|ITRI(:,2)==jsing|ITRI(:,3)==jsing);
   lus=loopnode(VER,ITRI,jsing);
   nb=length(lus);
   for kk=1:nb,
      k=icyc(kk-1,nb); l=kk; m=icyc(kk+1,nb);
      sa(kk)=solida([VER(lus(k),:);VER(lus(l),:);VER(lus(m),:)],[1 3 2],obs);
   end
   ndpos=sum(sa>0);
   ndneg=sum(sa<0);
   
   if ndpos==nb | ndneg==nb,% use spherical cap approximation
       % but only if theta, the cone angle of the cap as viewed from the origine
       % of the approximating sphere, is less than pi/6.
       % center of gravity of direct neighbours
       center=mean(VER(lus,:));
       % [jsing nb ndneg ndpos]
       % find radius of circle through direct neighbours
       % use: rms value of distances of direct neighbours to center
       radius=norm(VER(lus,:)-ones(nb,1)*center,'fro')/sqrt(nb);
       % distance from center to node(jsing)
       node2c=norm(VER(jsing,:)-center);
       coshalft=radius/sqrt(radius^2+node2c^2);
       theta=2*acos(coshalft);
       if theta < pi/6,
          if theta < 1.e-4,
              rowb(jsing)=pi*theta/2;,
          else,
              rowb(jsing)=4*pi*(1-coshalft)/theta;
          end
          if ndneg==nb, rowb(jsing)=-rowb(jsing);end
	       rowb(lus)=rowb(lus)+(2*pi-sum(rowb))/nb;
      else,
          rowb(jsing)=2*pi-sum(rowb);
	   end
   else,
	  rowb(jsing)=2*pi-sum(rowb);
   end
end
rowb=rowb/(2*pi);
