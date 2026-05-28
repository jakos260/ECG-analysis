% setobs.m
% script of electrog.m
% 200301729
% sets geom specs

% observation points for electrograms
skip=1;
if ntheta==17, skip=2; end
if ntheta==33, skip=3; end
selobs=sel(1:skip:ntheta);
OBS=robs*VERB(selobs,:);


% for simple electrog: skip this part and perhaps more
next=1;
if next==1,
z=0:.01:1.5*b;
nz=length(z);
NEEDLE=[zeros(nz,2) z'];
OBS=[OBS; NEEDLE] ;
end

xobs=OBS(:,1);
yobs=OBS(:,2);
zobs=OBS(:,3);
% critical angle for un-obscured observation from stimulus site
phi1=acos(a/b);

if r6==1 & exist('ntims'), 
    plotphi,
end