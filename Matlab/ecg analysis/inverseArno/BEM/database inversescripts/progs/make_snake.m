function [VER, ITRI,CENTERS] =make_snake(CORE, id, nphi,radius)
% make_snake.m
% generates triangulated mesh for a closed 3D object having local cylinder symmetry 
% local symmetry axis defined by a number of nodes along the CORE (size nx,3) of the
% snake-like geometry in 3D space, reconstructed at distances id (vector)
% nphi: number of equidistant phi values ;
% CENTERS: new interpolated center points.
% lfractions: positions of the rings circles) of the circles supporting the mesh
% simple closure at 'tim and tail' not yet implemented; function type implementation
% not yet included
% 2013-11-19; A. van Oosterom; 
% 2013-11-27 oostep1 Function


% id=[0:1.5:36]/36  % NB: subsets are the inter-electode distances 2:1.5 and 4*1.5
if ~exist('nphi') || isempty(nphi)
    nphi=12;
end
nf=length(id);

% CORE=CORE(1:3:end,:);
nc=size(CORE,1);

d=norm3d(CORE(2:nc,:)-CORE(1:nc-1,:));

cd1=[0;cumsum(d)];

CORE_F=spline(cd1(1:9),CORE(1:9,:)',cd1(1):0.05:cd1(9))';


ncc=size(CORE_F,1);
d=norm3d(CORE_F(2:ncc,:)-CORE_F(1:ncc-1,:));
cd2=[0;cumsum(d)];
cd2(end);


% id=id*36; %resampled distance along the spline curve

CENTERS=spline(cd2(1:ncc),CORE_F(1:ncc,:)',id')';

d3=norm3d(CENTERS(2:nf,:)-CENTERS(1:nf-1,:));
cd3=[0; cumsum(d3)];


% NOTE the original nine positions do not comply with the supposed 36 mm along
% the catheter

% specify basic (closed)circle
phi=(0:nphi)/nphi*2*pi;
x=cos(phi); y=sin(phi); z=zeros(1,nphi+1); % normal along z axis
% radius=0.5;
VERC=radius*[x' y' z'];

% shift to centers and line up there normals to local path


% tail:
normal=CENTERS(2,:)-CENTERS(1,:);
normal=-normal/norm(normal);
phi=atan2(normal(2),normal(1));
theta=acos(normal(3));
VER=rotash(VERC,[phi theta 0]/pi,CENTERS(1,:));

ALLVER=VER(1:nphi,:);
for i=2:nf-1,
     %body
    normal=(CENTERS(i+1,:)-CENTERS(i-1,:))/2;
    normal=-normal/norm(normal);
    phi=atan2(normal(2),normal(1));
    theta=acos(normal(3));
    VER=rotash(VERC,[phi theta 0]/pi,CENTERS(i,:));
    ALLVER=[ALLVER;VER(1:nphi,:)];
%     plot3(VER(:,1),VER(:,2),VER(:,3),'k+-')
end
    
%tip
normal=CENTERS(nf,:)-CENTERS(nf-1,:);
normal=-normal/norm(normal);
phi=atan2(normal(2),normal(1));
theta=acos(normal(3));
VER=rotash(VERC,[phi theta 0]/pi,CENTERS(nf,:));

ALLVER=[ALLVER;VER;];

VER=ALLVER;

nver=size(ALLVER);
nverwall=nphi*nf;

% create wall
ITRI=[];

listb=1:nphi,

for i=1:nf-1, % nz-1
    lista=listb;
    listb=lista+nphi;
    ADD=make_peel(VER,lista,listb);
    ITRI=[ITRI;ADD];
end

% if desired: both can be closed 




