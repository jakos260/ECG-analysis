% rowforw.m;
% function [rowb,jsing]=rowforw(VER,ITRI,obs)
% size(rowb)= 1,nver

% 2012-03-12; A. van Oosterom
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

nver=size(VER,1);
ntri=size(ITRI,1);

OMEGA=zeros(nver,ntri);
rowb=zeros(1,nver);
[OMEGA,index]=dsa(VER,ITRI,obs,.01);

for j=1:ntri,
    ij=ITRI(j,:);
    rowb(ij)=rowb(ij) + OMEGA(:,j)';
end

jsing=[];
% search singularity

jsing=find(norm3d(VER-ones(nver,1)*obs)/norm(VER,'fro')<1.e-6);

if isempty(jsing), jsing=0; end

if jsing~=0,% treat singularity (determine auto solid angle)
    
    lus=loopnode(ITRI,jsing);
    
    nb=size(lus,2);
    for kk=1:nb,
        k=icyc(kk-1,nb); l=kk; m=icyc(kk+1,nb);
        sa(kk)=solida([VER(lus(k),:);VER(lus(l),:);VER(lus(m),:)],[1 3 2],obs);
    end
    ndpos=sum(sa>0);
    ndneg=sum(sa<0);
    
    if ndpos==nb | ndneg==nb,% use spherical cap approximation
        % but only if theta, the cone angle of the cap as viewed from the origin
        % of the approximating sphere, is less than pi/6.
        % center of gravity of direct neighbours
        center=mean(VER(lus,:));
       
        % find rho: the radius of circle through direct neighbours
        % use: rms value of distances of direct neighbours to center
        rho=norm(VER(lus,:)-ones(nb,1)*center,'fro')/sqrt(nb);
        
        % distance from center to node(jsing)
        node2c=norm(VER(jsing,:)-center);
        coshalft=rho/sqrt(rho^2+node2c^2);
        theta=2*acos(coshalft);
        
        if theta <2*pi/3,
            if theta < 1.e-4,
                rowb(jsing)=pi*theta/2;
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
