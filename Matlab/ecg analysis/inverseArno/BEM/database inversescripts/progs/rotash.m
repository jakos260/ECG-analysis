% rotash.m
% function [NEW,ROTA]=rotash(OLD,[phi theta gamma],shift)
% specify ANGELS IN UNITS PI  !!!!!!!!
% rotation is followed by shift of coordinate system
% ROTA: the applied rotation matrix; : NEW=OLD*ROTA'
% rotation x,y,z  over theta, phi and gamma
% into x'y'z' followed by shift xyzs
%	Phi,Theta = angle between z' axis and z axis
%	Gamma=angle beween y' axis and line normal to z' axis in xy plane
% see: Jeffreys & Jeffreys: Meth of Math Phys; page 122-123
% if GAMMA=-PHI: rotate THETA around axis normal to zz' plane
% for positive alpha: lefthanded; counter clockwise rotation: 
% for cc rotating over alpha holding the  +x-axis take   PHI=0.5 THETA=-alpha  GAMMA=-0.5
% for cc rotating over alpha holding the +y-axis take    PHI=0   THETA= alpha  GAMMA= 0.
% for cc rotating over alpha holding the +z-axis take    PHI=0   THETA=0       GAMMA= alpha.
% ORIG=(NEW-ones(n,1)*shift)*ROTA reverts the operation ( n=size(OLD,1) )
% 2012-12-02; A. van Oosterom

function [NEW,ROTA]=rotash(OLD,angles,shift)
persistent previousa 
persistent ROT

if isempty(previousa), previousa=[inf inf inf]; end
if sum(angles==previousa)~=3,
	a=angles*pi;
	cp=cos(a(1));
	sp=sin(a(1));
	ct=cos(a(2));
	st=sin(a(2));
	cg=cos(a(3));
	sg=sin(a(3));
    ROT=ones(3,3);
    ROT(1,1)= cp*ct*cg-sp*sg;  ROT(1,2)=-cp*ct*sg-sp*cg;  ROT(1,3)= cp*st;
    ROT(2,1)= sp*ct*cg+cp*sg;  ROT(2,2)=-sp*ct*sg+cp*cg;  ROT(2,3)= sp*st;
    ROT(3,1)=-st*cg;           ROT(3,2)= st*sg;           ROT(3,3)= ct;
    previousa=angles;
end
ROTA=ROT;
NEW=OLD*ROT';
if nargin<3, return,end
ndim=size(OLD,1);
NEW=NEW+ones(ndim,1)*shift;
