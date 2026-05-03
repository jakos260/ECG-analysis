% rotashpre.m
% function [NEW,ROT]=rotash(OLD,[phi theta gamma],shift, [preshift])
% angles are to be expressed in units pi ( i.e. degrees/180)
% (preshift,) rotation and shift of coordinate system
% 20120611 oostep1: adapted from rotash.m

function [NEW,ROTA]=rotashpre(OLD,angles,shift,preshift)
persistent previousa 
persistent ROT
if ~exist('preshift')
    preshift=zeros(1,3);
end

if isempty(previousa)==1, previousa=[inf inf inf]; end
if sum(angles==previousa)~=3,
	a=angles*pi;
	cp=cos(a(1));
	sp=sin(a(1));
	ct=cos(a(2));
	st=sin(a(2));
	cg=cos(a(3));
	sg=sin(a(3));
    ROT=ones(3,3);
    ROT(1,1)= cp*ct*cg-sp*sg;
    ROT(1,2)=-cp*ct*sg-sp*cg;
    ROT(1,3)= cp*st;
    ROT(2,1)= sp*ct*cg+cp*sg;
    ROT(2,2)=-sp*ct*sg+cp*cg;
    ROT(2,3)= sp*st;
    ROT(3,1)=-st*cg;
    ROT(3,2)= st*sg;
    ROT(3,3)= ct;
    previousa=angles;
end
ROTA=ROT;
NEW=(OLD+ones(size(OLD,1),1)*preshift)*ROT';
dim=size(OLD);
NEW=NEW+ones(dim(1),1)*shift;
