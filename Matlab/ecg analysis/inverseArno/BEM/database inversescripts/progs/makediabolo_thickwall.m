% makediabolo.m
% variant of dolly.m; dedicated to axial-symmetric shapes
% aiming at equal node density
% triangulation based on Fourier representation  (SINE TRANSFORM) of
% the distance to the origin expressed in cylinder coordinates :r=r(z)
% r(z)=a0-a3*cos(3*(pi/2)*z/a0); z=[-a0 .... a0]; 
% w is the desired wall thicknes,
% resulting in an expanded version of the surface specified by a0 and a3

% rho(z)=sqrt(r^2(z)-z^2)

% specify reconstruction parms: kr  for the desired z- levels           
%                               lr  for the number of phi values
% set nehs=1 for alternate shift over phi/2 of nodes at subsequent crossections

% 
clear
  kr=29; nehs=1;
  
  a0=18; % note: outer layer 
  a3=.3*a0;
  %a3=0;
  w=2;
  % note: outer layer has 'a0'=a0+w

% prepare; search for z values producing triangles having roughly equal
% size
figure(1)

clf
z=-a0:0.0001*a0:a0;
z=max(z,-a0); z=min(z,a0);
c=3*pi/2/a0;

nz=length(z);
rz=a0-a3*cos(c*z);
rhoz=sqrt(rz.^2-z.^2);

if w > 0,
  drhozdz=diffrows(rhoz)*(nz-1)/(2*a0);
  alpha=atan(drhozdz); % the angle of the local slope to the rhoz curve
  z=z-w*sin(alpha);
  rhoz=rhoz+w*cos(alpha);
  rz=sqrt(rhoz.^2+z.^2);
end


plot(z,rz)
xlabel('z');
ylabel(' blue: r(z); green: rho(z)')
title(' blue: r(z); green: rho(z)')
hold on
plot(z,rhoz,'g')
[ma ima]=max(rhoz);
[z(ima) rhoz(ima)];
pause

figure(2)
clf
  
   nz=size(z,2);
   % compute the length along the rhoz curve
   s(1)=0;
   for i=2:nz,
       s(i)=sqrt((z(i)-z(i-1))^2+(rhoz(i)-rhoz(i-1))^2);
   end
   scumula=cumsum(s);
   plot(z,5*scumula,'k')
   hold on

zk(1)=z(1); zk(kr)=z(nz) ;alphak(1)=pi/2; alphak(kr)=-pi/2;
rhok(1)=0;
rhok(kr)=0;
levs=scumula(nz)*(0:1/(kr-1):1);

for k=2:kr-1,
   j=find(scumula(1:nz-1)<=levs(k)& scumula(2:nz)>levs(k));
   plot([z(1) z(j) z(j)], 5*[scumula(j) scumula(j) 0],':')
   zk(k)=z(j);
   rhok(k)=rhoz(j);
end
title('distance along the curve as a function of z')
xlabel('z;  (mm)')
ylabel('distance laong the curve (mm)')
pause

figure(1)
hold on
plot(zk,rhok,'r*')
title(' blue: r(z); green: rho(z); asterisks at equal distances along the curve')

pause

zk(2)=zk(2)*0.993 % small tuning applied to reduce tha max/min area ratio
zk=(zk-reverse(zk))/2;

next=0;
if next==1,
   lrk=[1 6 8 14 18 24 29 30 32 33 32 29 26 23 23];
   nlkh=size(lrk,2)-1;
   lrk=[lrk reverse(lrk(1:nlkh))];
end

next=1;
if next==1,
   lrk=round(rhok/max(rhok)*35);
   lrk(1)=1; lrk(kr)=1;
   lrk(2)=6; lrk(kr-1)=6;
end

figure(3)
clf
  
VER(1,1:3)=[0 0 zk(1)];
ind=1;
  
for iz=2:kr-1,
    lr=lrk(iz);
    fphi=2*pi/lr;
    z(iz)=zk(iz);
    rho(iz)=rhok(iz);
    heven=((rem(iz+1,2)-1)/2)*nehs;
    for jfi=1:lr,
        afi=jfi-heven;
        phi=fphi*afi;
        ind=ind+1;
        VER(ind,1:3)=[rho(iz)*cos(phi) rho(iz)*sin(phi) z(iz)];
    end
end
rho(kr)=0;
ind=ind+1;
VER(ind,1:3)=[0 0 zk(kr)];
nver=ind
plot3(VER(:,1),VER(:,2),VER(:,3),'+')

itri=0;

levelend=cumsum(lrk);
levelbeg=[1 levelend(1:kr-1)+1];

% create cap
k=1;

listb=[levelbeg(k+1):levelend(k+1) levelbeg(k+1)];
nb=length(listb)-1;
ITRI(1:nb,1:3)=[ones(nb,1) listb(1:nb)' listb(2:nb+1)'];
ntri=nb;

    for k=2:kr-2,
  
   % triangulate peel(k)
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
     %[k na1 nb1 na nb]
     lista=[lista lista(1)];
     %pause
     
     %identify nearest neighbour of lista(1)    
     d=norm3d(ones(nb,1)*VER(lista(1),:)-VER(listb,:));
     [mi imi]=min(d);
     
     if imi>1, listb=listb([imi:nb 1:imi-1]);end
     listb=[listb listb(1)];
     %pause
     
     
     while length(lista)>1,
         while length(listb)>1,
             d1=norm3d(VER(lista(1),:)-VER(listb(2),:));
             d2=norm3d(VER(lista(2),:)-VER(listb(1),:));
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
             if length(listb)<=1, break, end
             if length(lista)<=1, break, end
             %pause
         end
        if length(listb)<=1, break, end
        if length(lista)<=1, break, end
     end
     
     ntri=ntri+1;
     ITRI(ntri,1:3)=[lista reverse(listb)];
     
     % create cap
     listb=[levelbeg(kr-1):levelend(kr-1) levelbeg(kr-1)];
     nb=length(listb)-1;
     ADD(1:nb,1:3)=[nver*ones(nb,1) listb(2:nb+1)' listb(1:nb)'];
     ITRI=[ITRI;ADD];
 end
 
figure(3)
clf

area=nodearea(VER,ITRI);
VALS=area;
triplot
sum(area)
[mean(area) std(area) sum(area) std(area)/mean(area) max(area)/min(area)]









