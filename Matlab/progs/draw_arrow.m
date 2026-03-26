% draw_arrow.m
% [VER,ITRI]=draw(arrow(tail,head,frl,frr);
% draw an arrow from row-vectors tail to head
% frl is the fraction of the total length occupied by the length of the arrowhead
% frw is the fraction of the total length of the span of the arrowhead
% defaults: frl=0.1; frr=0.05
% if size(tail,2)==size(head,2)>2, a 3D version is constructed
% the latter is based on on  nphi equidistant phi values and nz z values;

% defaults:[0.1 3 1 0.01 3 16 ]
% 2012-01-11; A. van Oosterom


function [VER,ITRI]=draw_arrow(tail,head,frl,frr)

if nargin<3
   frl=0.2;
   frr=0.02;
end
ITRI=[];

if size(head,2)>2 && size(tail,2)>2
    ll=norm(head-tail); %length of the arrow 
    phi=atan2(head(2)-tail(2),head(1)-tail(2))/pi;
    theta=acos((head(3)-tail(3))/ll)/pi;
    
   % 3D arrow
   % parameters used for creation 3D arrow 
   % plot using triplot
   
   nphi=16;
   nz=2;
   [VER,ITRI]=make_cylinder(nphi,nz);
   nver=size(VER,1);
   VER(2:nver-1,3)=VER(2:nver-1,3)*(1-frl)*ll;
   VER(nver,3)=ll;
   VER(2:nver-1,1:2)=VER(2:nver-1,1:2)*frr*ll;
   VER=rotash(VER,[phi,theta,0],tail);
else
    % 2D ; plot results VER
    ll=norm(head(1:2)-tail(1:2));
    phi=0;
    theta=acos((head(2)-tail(2))/ll)/pi;
    VER=[0 0 0;
         0 ll 0;
        -ll*frr (1-frl)*ll 0;
         0  (1-frl/2)*ll  0 ;
        ll*frr (1-frl)*ll 0;
        0 ll 0];
        VER=rotash(VER,[0 0 -theta],[tail 0]);
    
end
    
    
      
      
     
         
        
    
    
    
    
    
  

   
   




