% make_snake.m
% generates triangulated mesh for a closed 3D object having local cylinder symmetry 
% local symmetry axis defined by a number of nodes along the CORE (size nx,3) of the
% snake-like geometry in 3D space
% nphi: number of equidistant phi values ;
% lfractions: positions of the rings circles) of the circles supporting the mesh
% simple closure at 'tim and tail' not yet implemented; function type implementation
% not yet included
% 2013-11-19; A. van Oosterom; 


fractions=[0:1.5:36]/36  % NB: subsets are the inter-electode distances 2:1.5 and 4*1.5
nf=size(fractions,2);

nphi=12;

CORE=loadmat('posi.mat');

figure(1)
clf
CORE=CORE(1:3:end,:);
nc=size(CORE,1);
plot3(CORE(:,1),CORE(:,2),CORE(:,3),'k-+')
axis equal
axis square
grid on

hold on
plot3(CORE(1,1),CORE(1,2),CORE(1,3),'r*')

d=norm3d(CORE(2:nc,:)-CORE(1:nc-1,:))

cd1=[0;cumsum(d)]

CORE_F=spline(cd1(1:9),CORE(1:9,:)',cd1(1):0.05:cd1(9))';

plot3(CORE_F(:,1),CORE_F(:,2),CORE_F(:,3),'m')

ncc=size(CORE_F,1);
d=norm3d(CORE_F(2:ncc,:)-CORE_F(1:ncc-1,:));
cd2=[0;cumsum(d)];
cd2(end)

pause

s=fractions*36 %resampled distance along the spline curve

CENTERS=spline(cd2(1:ncc),CORE_F(1:ncc,:)',s')';

d3=norm3d(CENTERS(2:nf,:)-CENTERS(1:nf-1,:));
cd3=[0; cumsum(d3)];

plot3(CENTERS(1:nf,1),CENTERS(1:nf,2),CENTERS(1:nf,3),'g*');
plot3(CENTERS(1,1),CENTERS(1,2),CENTERS(1,3),'r*');


% NOTE the original nine positions do not comply with the supposed 36 mm along
% the catheter

% specify basic (closed)circle
phi=(0:nphi)/nphi*2*pi;
x=cos(phi); y=sin(phi); z=zeros(1,nphi+1); % normal along z axis
radius=0.5;
VERC=radius*[x' y' z'];

% shift to centers and line up there normals to local path


% tail:
normal=CENTERS(2,:)-CENTERS(1,:);
normal=-normal/norm(normal);
phi=atan2(normal(2),normal(1));
theta=acos(normal(3));
VER=rotash(VERC,[phi theta 0]/pi,CENTERS(1,:));
plot3(VER(:,1),VER(:,2),VER(:,3),'k+-')
pause

ALLVER=VER(1:nphi,:);
for i=2:nf-1,
     %body
    normal=(CENTERS(i+1,:)-CENTERS(i-1,:))/2;
    normal=-normal/norm(normal);
    phi=atan2(normal(2),normal(1));
    theta=acos(normal(3));
    VER=rotash(VERC,[phi theta 0]/pi,CENTERS(i,:));
    ALLVER=[ALLVER;VER(1:nphi,:)];
    plot3(VER(:,1),VER(:,2),VER(:,3),'k+-')
end
    
%tip
normal=CENTERS(nf,:)-CENTERS(nf-1,:);
normal=-normal/norm(normal);
phi=atan2(normal(2),normal(1));
theta=acos(normal(3));
VER=rotash(VERC,[phi theta 0]/pi,CENTERS(nf,:));
plot3(VER(:,1),VER(:,2),VER(:,3),'k+-')
pause

ALLVER=[ALLVER;VER;];

VER=ALLVER;

nver=size(ALLVER)
nverwall=nphi*nf

% create wall
ITRI=[];

listb=1:nphi,

for i=1:nz-1,
    lista=listb
    listb=lista+nphi
    ADD=make_peel(VER,lista,listb);
    ITRI=[ITRI;ADD];
end

figure(2)
clf
center
VALS=VER;
triplot

% if desired: both can be closed 




