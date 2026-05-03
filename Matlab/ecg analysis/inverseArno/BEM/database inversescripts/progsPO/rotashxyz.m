% rotashxyz.m
% function [NEW,ROT]=rotash(OLD,[phi  theta psi],shift)
% angles are to be expressed in units pi
% rotation and shift of coordinate system
% ROT: the applied rotation matrix; (inv(ROT)*NEW')' reverts the rotation  
% rotation x,y,z  over theta, phi and gamma
% into x'y'z' followed by shift xyzs
% Subsequent rotation phi over x, theta over y' and psi over z'', followed
% by translation shift

% 20120821 oostep1. Adapted from rotash 2005-02-24; A. van Oosterom

function [NEW,ROTA]=rotashxyz(OLD,angles,shift)

% %block for testing purposes
% global USEZYZ
% if USEZYZ
%     [NEW,ROTA]=rotash(OLD,angles,shift);
%     return
% end
if max(angles)>2 || min(angles)<-2
    warning('large angle detected. Remember angles are in units pi.');end

persistent previousa 
persistent ROT
if isempty(previousa)==1, previousa=[inf inf inf]; end
if sum(angles==previousa)~=3,
	a=angles*pi;
	cph=cos(a(1));
	sph=sin(a(1));
	ct=cos(a(2));
	st=sin(a(2));
	cps=cos(a(3));
	sps=sin(a(3));
    ROT=ones(3,3);
    ROT(1,1)= ct*cps;
    ROT(1,2)=-cph*sps+sph*st*cps;
    ROT(1,3)= sph*sps+cph*st*cps;
    ROT(2,1)= ct*sps;
    ROT(2,2)= cph*cps+sph*st*sps;
    ROT(2,3)=-sph*cps+cph*st*sps;
    ROT(3,1)=-st;
    ROT(3,2)= sph*ct;
    ROT(3,3)= cph*ct;
    previousa=angles;
end
ROTA=ROT;
NEW=OLD*ROT';
dim=size(OLD);
NEW=NEW+ones(dim(1),1)*shift;
