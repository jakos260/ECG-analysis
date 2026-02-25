% electrog.m
% 20030205
% computes epicardial and endocardial types of electrograms
% on thick sphere model of atria
% activation according to Huygens
% A. van Oosterom
warning off
clf
clear

% settings for wave front computations
nphi=101;
phi=[0:1/(nphi-1):1];
phi=phi*pi;

% geometry of front:
[VER,ITRI]=loadtri('annulus.tri');
nver=length(VER);

% asuming first half of these are on the inner ring:
ainit=max(max(abs(VER(1:nver/2,:))));
binit=1;
a=ainit;
b=binit;

sigp=1;
sigm=3;
sigs=1;
sigobs=1;

% specification of boundary
ntheta=17; nphib=20;
ntheta=64; nphib=20;

%ntheta=9; nphib=16;

[VERB, ITRIB]=makeglobe(ntheta,nphib); % unit sphere
% used to represent endocardium
nverb =length(VERB);
ntrib =length(ITRIB);

% compute B matrix
% use symmetry related to axial source config and spherical
% boundary with radius a; since B contains solid angles its values are 
% invariant under scaling of a.

B=zeros(nverb,nverb);
nhalf=nverb/2;
ref=nverb+1-[1:nverb];
sel=[1 2:nphib:nverb];
for i=1:nhalf,
    [B(i,1:nverb),jsing]=rowforw(VERB,ITRIB,VERB(i,:));
    B(nverb+1-i,ref)=B(i,:);
end

% observation points for electrograms
nn=9;
robs=a;

% compute T matrix; used to compute potentials at the boundary
   C=(sigm-sigp)/(sigm+sigp)*B;
   AMA=inv(eye(nverb)-C);
   T=zeros(ntheta,ntheta);
   for i=1:ntheta,
     T(i,1)=AMA(sel(i),1);
     T(i,ntheta)=AMA(sel(i),sel(ntheta));
     for j=2:ntheta-1,
       T(i,j)=sum(AMA(sel(i),sel(j):sel(j+1)-1));
     end
   end 

% uibox settings
r1=1; r2=0; r3=0; r4=0; r5=1; r6=1; r7=0; r8=1; r9=0;
inho=0;

setgeom
scalinit=.6;
scal=scalinit;
sl1=1;
ntims=101;
t=1:ntims;
tijd=0.55+[1:ntims]/ntims;

stim='endo';
gettiming
getphi
plotphi

uiboxe1=[.85 .95 .15 .05]; % set inner radius

uiboxp1=[.85 .9 .15 .05];  % getradius 

uibox1r1=[0   .65  .08 .05]; % endostim
uibox1r2=[0   .70 .08 .05]; % epistim
uibox1r3=[.9  .5  .08 .05]; % inhom
%uibox1r4=[.9  .45 .08 .05]; % partepi
%uibox1r5=[.9  .4  .08 .05]; % partall
uibox1r6=[.0  .9  .12 .05]; % showpots
uibox1r7=[.0  .85 .12 .05]; % showfronts
uibox1r8=[.0  .00 .12 .05]; % obsendo (r8)
uibox1r9=[.0  .05 .12 .05]; % obsepi

uiboxs1=[.9  .05 .03 .3];   % scalpots
uiboxs2=[.09  .25  .03 .3];   % setrobs

uiboxt1=[0    .75  .12 .05]; % stimulus
uiboxt2=[0    .95 .12 .05];  % display
uiboxt3=[0.   .2 .12 .05];   % display r-obs
uiboxt4=[0.   .1   .12 .05]; % observe
uiboxt5=[0.9  .0  .08 .05];  % scale
uiboxt6=[0.  .5  .09 .05];   % r-obs

setuis



