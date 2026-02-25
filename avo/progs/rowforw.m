% rowforw.m;
% function [rowb,jsing,index]=rowforw(VER,ITRI,obs)
% size(rowb)= 1,nver

% 2015-07-20; A. van Oosterom

%  compute (part of a) row of Bmatrix:
%  distributed solid angles as seen by obs subtended
%  by triangles ITRI ; VER are its vertices
% ITRI MUST be CLOSED, non-self intersecting

%  called by forward.m

% functions called: dsa.m      (general distributed solid angle function)
%                   loopnode.m (identifies loop of direct vertex neighbours
%		                         around node; orientation clockwise when viewed from
%		                         outside)
%                   solida.m   (solid angle function)

% 2015-07-16 index added to output

function [rowb,jsing,index]=rowforw(VER,ITRI,obs)

nver=size(VER,1);
ntri=size(ITRI,1);

[OMEGA,index] = dsa(VER,ITRI,obs,0.01);   % NOTE: uses dsa

rowb=zeros(1,nver);
for j=1:ntri
    ij=ITRI(j,:);
    rowb(ij)=rowb(ij) + OMEGA(:,j)';
end
rowb=rowb/(2*pi);

% search singularity
jsing = find(norm3d(VER-ones(nver,1)*obs)/norm(VER,'fro')<1.e-6);

if ~isempty(jsing) % treat singularity (determines auto solid angle)
    
    loop=loopnode(ITRI,jsing);
    
    nb=size(loop,2);
    sa = zeros(nb,1);
	index = zeros(nb,1);
    for kk=1:nb
        k=icyc(kk-1,nb); 
        l=kk; 
        m=icyc(kk+1,nb);
        [sa(kk),index(kk)]=solida([VER(loop(k),:);VER(loop(l),:);VER(loop(m),:)],[1 3 2],obs);
    end
    
    ndpos=sum(sa>=0);
    ndneg=sum(sa<=0);
    si = sum(index>0);
    
    
    if si~=nb && (ndpos==nb || ndneg==nb) % use spherical cap approximation
        
        if si>0
            a=1;
        end
        % but only if theta, the cone angle of the cap as viewed from the origin
        % of the approximating sphere, is less than pi/6.
        % center of gravity of direct neighbours
        center=mean(VER(loop,:));
       
        % find rho: the radius of circle through direct neighbours
        % use: rms value of distances of direct neighbours to center
        rho = norm(VER(loop,:) - ones(nb,1)*center,'fro')/sqrt(nb);
        
        % distance from center to node(jsing)
        node2c = norm(VER(jsing,:)-center);
        coshalft = rho/sqrt(rho^2 + node2c^2);
        theta    = 2.0 * acos(coshalft);
        
        if theta > pi/6% 2*pi/3
            rowb(jsing) = 1 - sum(rowb);
        else
            if theta > 1.e-4
                rowb(jsing) = 2 * ( 1.0 - coshalft )/theta;
            else
                rowb(jsing) = theta/4;  %/2 avo???
            end
            if ndneg==nb 
                rowb(jsing)=-rowb(jsing);
            end
            rowb(loop) = rowb(loop) + (1-sum(rowb))/nb;           
        end
    else
        rowb(jsing)=1-sum(rowb);
    end
else
    jsing = 0;
end





%%


% rowforw.m;
% function [rowb,jsing,index]=rowforw(VER,ITRI,obs)
% size(rowb)= 1,nver

% 2015-07-20; A. van Oosterom

%  compute (part of a) row of Bmatrix:
%  distributed solid angles as seen by obs subtended
%  by triangles ITRI ; VER are its vertices
% ITRI MUST be CLOSED, non-self intersecting

%  called by forward.m

% functions called: dsa.m      (general distributed solid angle function)
%                   loopnode.m (identifies loop of direct vertex neighbours
%		                         around node; orientation clockwise when viewed from
%		                         outside)
%                   solida.m   (solid angle function)

% 2015-07-16 index added to output

% function [rowb,jsing,index]=rowforw(VER,ITRI,obs)
% 
% nver=size(VER,1);
% ntri=size(ITRI,1);
% 
% % OMEGA=zeros(nver,ntri);
% rowb=zeros(1,nver);
% 
% [OMEGA,index]=dsa(VER,ITRI,obs,.01);   % NOTE: uses dsa
% 
% for j=1:ntri
%     ij=ITRI(j,:);
%     rowb(ij)=rowb(ij) + OMEGA(:,j)';
% end
% 
% % jsing=[];
% 
% % search singularity
% 
% jsing = find(norm3d(VER-ones(nver,1)*obs)/norm(VER,'fro')<1.e-6);
% 
% if isempty(jsing), jsing=0; end
% 
% 
% if jsing == 280
%     stop=1;
% end
% if jsing~=0 % treat singularity (determines auto solid angle)
%     
%     loop=loopnode(ITRI,jsing);
%     
%     nb=size(loop,2);
%     sa = zeros(nb,1);
% 	index = zeros(nb,1);
%     for kk=1:nb
%         k=icyc(kk-1,nb); l=kk; m=icyc(kk+1,nb);
%         [sa(kk),index(kk)]=solida([VER(loop(k),:);VER(loop(l),:);VER(loop(m),:)],[1 3 2],obs);
%     end
%     
%     ndpos=sum(sa>=0);
%     ndneg=sum(sa<=0);
%     
%     if ndpos==nb || ndneg==nb % use spherical cap approximation
%         % but only if theta, the cone angle of the cap as viewed from the origin
%         % of the approximating sphere, is less than pi/6.
%         % center of gravity of direct neighbours
%         center=mean(VER(loop,:));
%        
%         % find rho: the radius of circle through direct neighbours
%         % use: rms value of distances of direct neighbours to center
%         rho = norm(VER(loop,:) - ones(nb,1)*center,'fro')/sqrt(nb);
%         
%         % distance from center to node(jsing)
%         node2c=norm(VER(jsing,:)-center);
%         coshalft=rho/sqrt(rho^2+node2c^2);
%         theta=2*acos(coshalft);
%         
%         if theta <= pi/6% 2*pi/3
%             if theta <= 1.e-4
%                 rowb(jsing)=pi*theta/2;  %/4???
%             else
%                 rowb(jsing)=4*pi*(1-coshalft)/theta;
%             end
%             if ndneg==nb 
%                 rowb(jsing)=-rowb(jsing);
%             end
%             rowb(loop) = rowb(loop) + (2*pi-sum(rowb))/nb;
%         else
%             rowb(jsing)=2*pi-sum(rowb);
%         end
%     else
%         rowb(jsing)=2*pi-sum(rowb);
%     end
% end
% 
% rowb=rowb/(2*pi);


