% make_diabolo.m
% variant of dolly.m; dedicated to axial-symmetric shapes
% aiming at equal node density
% triangulation based on Fourier representation  (SINE TRANSFORM) of
% the distance to the origin expressed in cylinder coordinates: r=r(z)
% for diabolo: r(z)=a0-a3*cos(3*(pi/2)*z/a0); z=[-a0 .... a0]; 
% z levels selected forcing equal increments dels along the s-curve 
% the number of phi values along the circles rho(z)=sqrt(r^2(z)-z^2) of cross-section
% at z-level is set such that rho(z)*2*pi/lr~dels
% A. van Oosterom 2006_11_6

clear all
% 
%   a0=30; 
%   fra3=0.6
%   a3=fra3*a0;
%   c=3*pi/2/a0;
%   dels=0.566 % used for diabolo; a0=30; fra3=0.6 ; 301 levels;
%              nver=83689; ntris=167374

  % search for z values producing equal s increments: dels
  % tune dels for obtaining the desired number of segments dels along s

% 
%   a0=30; 
%   fra3=0
%   a3=fra3*a0;
%   c=3*pi/2/a0;
%   dels=0.469;  % used for sphere; a0=30; 201 levels
%   % nver=51402; ntri=102800 
  
  
  
  a0=30; 
  fra3=0
  a3=fra3*a0;
  c=3*pi/2/a0;
  dels=2.9; 
  % used for sphere; a0=3;  33 levels
  
  
  

% for alternate shift over phi/2 of nodes at subsequent crossections use:  nehs=1;
  nehs=0;

% use finetune=1 for fine tuning; use only if symmetry about xy plane is
% involved, else, put: finetune=0
finetune=1;



s(1)=0;
k=1;
zk(1)=-a0;
rhok(1)=0;
epsilon=1.e-6;
zl=-a0;
zr=-a0+dels;
incr(1)=dels;

while zk(k)<=a0,
    % find x such that zk(k)+x produces the desired increment dels
    % uses bi_section algorithm
    i=1;
    while i<50,
       i=i+1;
       z=(zl+zr)/2;
       % function call
       rz=a0-a3*cos(c*z);
       rho=sqrt(rz^2-z^2);
       increment=sqrt((rho-rhok(k))^2+(z-zk(k))^2); %function; value should be: dels
       test_incr=increment-dels;
       if abs(test_incr)<epsilon,
           k=k+1;
           zk(k)=z;
           rhok(k)=rho;
           zl=z; zr=z+dels;
           incr(k)=increment;
           break
       end
       if test_incr>0,
           zr=z;
       else,
           zl=z;
       end
    end
    if i==50, break, end 
end 

zk
nsegs=size(zk,2)

'length of s curve'
slength=sum(incr)

'adjust dels, and rerun,  if final zk value shown is not close to a0; or for setting desired n_zlevels '
pause

% plots to check geometry and segmentation
delz=2*a0/10000; % merely use for plotting
z=[-a0:delz:a0]; 
nz=length(z);
rz=a0-a3*cos(c*z);
rhoz=sqrt(rz.^2-z.^2);

figure(1)
clf
plot(z,rz)
xlabel('z');
ylabel(' blue: r(z); green: rho(z)')
title(' blue: r(z); green: rho(z)')
hold on

if finetune==1,
   nzk=length(zk);
   zk=-a0+2*a0/(zk(nzk)+a0)*(zk+a0);
   zk=((zk'-reverse(zk)')/2)';
   % recompute rhok
   rk=a0-a3*cos(c*zk);
   rhok=sqrt(rk.^2-zk.^2);
end
  
plot(z,rhoz,'g')
plot(zk,rhok,'r*')
title(' blue: r(z); green: rho(z); asterisks at equal distances along the curve')

% compute nodes
   kr=length(zk);
   lrk=round(2*pi*rhok/dels);
   %lrk=max(lrk,lrk(2));
   lrk(1)=1; lrk(kr)=1;   
   % lrk(2)=6; lrk(kr-1)=6;
   % ranges indices vertices at the different z levels  
   levelend=cumsum(lrk);
   levelbeg=[1 levelend(1:kr-1)+1];
   %levelbeg
   %savemat(['levelindices_0' num2str(10*fra3) '.lst'],levelbeg)
   %pause
   VER(1,1:3)=[0 0 zk(1)];
   ind=1;
   for iz=2:kr-1,
     lr=lrk(iz);
     fphi=2*pi/lr;
     z(iz)=zk(iz);
     rho(iz)=rhok(iz);
     jfi=1:lr;
     jfi=jfi-rem(iz,2)*nehs/2;
     phi=fphi*jfi';
     VER(ind+1:ind+lr,1:3)=[rho(iz)*cos(phi) rho(iz)*sin(phi) z(iz)*ones(lr,1)];
     ind=ind+lr;
   end

   rho(kr)=0;
   ind=ind+1;
   VER(ind,1:3)=[0 0 zk(kr)];
   nver=ind;

   % figure(3)
   % clf

' creating  triangulation'
itri=0;

'creating lower polecap'
k=1;
listb=[levelbeg(k+1):levelend(k+1) levelbeg(k+1)];
nb=length(listb)-1;
ITRI(1:nb,1:3)=[ones(nb,1) listb(1:nb)' listb(2:nb+1)'];
ntri=nb;

% triangulate peels(k), k=2:kr-2

for k=2:kr-2,
    ['creating triangulation of peel(' num2str(k) ')']
    lista=levelbeg(k):levelend(k);
    listb=levelbeg(k+1):levelend(k+1);
    na1=length(lista);
    nb1=length(listb);
    if na1>nb1,
       temp=listb;
       listb=lista; lista=temp;
    end
    na=length(lista);
    nb=length(listb);
    lista=[lista lista(1)];
     
    %identify nearest neighbour of lista(1)    
    d=norm3d(ones(nb,1)*VER(lista(1),:)-VER(listb,:));
    [mi imi]=min(d);
    if imi>1, listb=listb([imi:nb 1:imi-1]);end
    listb=[listb listb(1)];
    
    while length(lista)>1,
       while length(listb)>1,
             d1=norm3d(VER(lista(1),:)-VER(listb(2),:));
             if length(lista)>1, 
                 d2=norm3d(VER(lista(2),:)-VER(listb(1),:));,
             else,
                 d2=d1;
             end
             %pause
             ntri=ntri+1;
             
             if d1<=d2,
                if na1<=na,
                     ITRI(ntri,1:3)=[lista(1) listb(1) listb(2)];
                 else,
                     ITRI(ntri,1:3)=[lista(1) listb(2) listb(1)];
                 end
                 listb(1)=[];
              else,
                 if na1<=na,
                     ITRI(ntri,1:3)=[lista(2) lista(1) listb(1)];
                 else,
                     ITRI(ntri,1:3)=[lista(2) listb(1) lista(1)];
                 end
                 lista(1)=[];
             end

             if length(lista)+length(listb)==3,
                trits_ntri=ITRI(ntri,:);
                ntri=ntri+1;
                trits=[lista listb];
               
                % test if order of the indices needs to be reversed
                for i=1:3,
                    test=trits_ntri([i icyc(i+1,3) icyc(i+2,3)]);
                    if sum(test==trits)==2,
                        trits([2 1])=trits(1:2);
                        break
                    end
                end      
                 ITRI(ntri,:)=trits;        
                    
                lista=[];
                listb=[];
             end
             if length(listb)+length(lista)<3, break, end
         end
     end
 end
 
'creating north pole cap'
listb=[levelbeg(kr-1):levelend(kr-1) levelbeg(kr-1)];
nb=length(listb)-1;
ADD(1:nb,1:3)=[nver*ones(nb,1) listb(2:nb+1)' listb(1:nb)'];
ITRI=[ITRI;ADD];

figure(3)
clf

area=nodearea(VER,ITRI); 
VALS=area;
triplot
sumarea=sum(area)
[ma ima]=max(area)
[mi imi]=min(area)
'[mean(area) std(area) std/mean     mi        ma         mi/ma]'
quality=[mean(area) std(area) std(area)/mean(area) mi ma mi/ma]
% use triareas(VER,ITRI) if desired








