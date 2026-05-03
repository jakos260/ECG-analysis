% getqrs.m
% extract qrs interval from **.qst data
% old Nijmegen data in which 50 samples of QRS and 50 samples of STT 
% are combined;
% shift data to get maximum correlation with reference
% shift to zero mean

clear
allfiles=loadchars('normals.labs');
nfiles=size(allfiles);

% first: get reference data
% =r01904=gu
file=allfiles(48,:);
PHI=loadmat(file);
PHI=PHI';
%PHI=zeromean(PHI);
n=size(PHI);
nlds=n(1);
ntims=n(2);

nt=ntims/2;
%nt=ntims;

% 10 labels for extracting 12-lead signals:
% for old nim system:
% V1,V2,(2 labels for V3), V4 to V6, VR, VL and ~VF
lab=[19 26 33 34 41 48 54 1 2 56];

% 10 labels for extracting 12-lead signals:
% for ams system:
% V1,V2,(2 labels for V3), V4, V5, V6,  VR, VL and ~VF
%lab=[12 18 25 25 31 40 45 64 63 48];

% establish reference
VFWCT=-(PHI(lab(8),:)+PHI(lab(9),:));
VFest=PHI(lab(10),:);
rd=reldiff(VFWCT,VFest);
reldif=rd(4);
if reldif < .5
VF=VFWCT;
'wct established'
wct=1;
else 
VF=VFest;
'VF estimated '
wct=0;
end


for nofile=1:nfiles(1),

clf
color='blue';
funscal=.5;
t=1:nt;
tbeg=1;
tend=nt;
t=t*.8/nt;
axis([0 4  0 4])
axis manual
axis('off')
hold on

for i=1:3,
for j=1:4,
k=j+(i-1)*4;
yshift=4-i;
xshift=j-.8;
if k == 1,
sig=PHI(lab(9),tbeg:tend)-PHI(lab(8),tbeg:tend);
tekst=sprintf('%s','I');
elseif k == 2,
sig=PHI(lab(9),tbeg:tend);
if wct==1,
tekst=sprintf('%s','aVL');
else
tekst=sprintf('%s','Vla');
end
elseif k == 3,
sig=PHI(lab(1),tbeg:tend);
tekst=sprintf('%s','V1');
elseif k == 4,
sig=PHI(lab(5),tbeg:tend);
tekst=sprintf('%s','V4');
elseif k == 5,
sig=VF(1,tbeg:tend)-PHI(lab(8),tbeg:tend);
tekst=sprintf('%s','II');
elseif k == 6,
sig=PHI(lab(8),tbeg:tend)
if wct==1,
tekst=sprintf('%s','aVR');
else 
tekst=sprintf('%s','Vra');
end
elseif k == 7,
sig=PHI(lab(2),tbeg:tend);
tekst=sprintf('%s','V2');
elseif k == 8,
sig=PHI(lab(5),tbeg:tend);
tekst=sprintf('%s','V5');
elseif k == 9,
sig=VF(1,tbeg:tend)-PHI(lab(9),tbeg:tend);
tekst=sprintf('%s','III');
tekst=sprintf('%s','III');
elseif k == 10,
sig=VF(1,tbeg:tend);
if wct==1
tekst=sprintf('%s','aVF');
else
tekst=sprintf('%s','Vlf');
end
elseif k == 11,
sig=.5*PHI(lab(3),tbeg:tend)+.5*PHI(lab(4),tbeg:tend);
tekst=sprintf('%s','V3');
elseif k == 12,
sig=PHI(lab(7),tbeg:tend);
tekst=sprintf('%s','V6');
end
plot(t+xshift,funscal*sig+yshift,color)
text(xshift,yshift+.5,tekst);
end
end
set(text(0.,0.2,sprintf('%s',file)),'color','blue');
xl=[1.2 1.2 2 ];
yl=[.2+funscal .2 .2];
set(line(xl,yl),'color','k')
text(1.25,.2+0.8*funscal,sprintf('%s','1 mV'));
text(1.8,.35,sprintf('%s','100 ms'));
if wct==1,
set(text(0,0,sprintf('%s','(wct)')),'color','k');
else
set(text(0,0,sprintf('%s','(zme)')),'color','k');
end

% now: get data to be processed
file2=allfiles(nofile,:);
PSI=loadmat(file2);
PSI=PSI';
n2=size(PSI);
nlds2=n2(1);
ntims2=n2(2);

nt2=ntims2;
%nt2=ntims2/2;
%PSI=zeromean(PSI);

% establish reference
VFWCT=-(PSI(lab(8),:)+PSI(lab(9),:));
VFest=PSI(lab(10),:);
reldif=norm(VFWCT-VFest)/norm(VFest)
if reldif < .5
VF=VFWCT;
'wct established'
wct=1;
else 
VF=VFest;
'VF estimated '
wct=0;
end

best=1.e+6;
TEST=zeros(64,nt2+50);
TEST(:,25:nt2+24)=PSI(:,1:nt2);

VFTEST=zeros(1,nt2+50);
VFTEST(1,25:nt2+24)=VF(1,1:nt2);

for tt=20:30,
tbeg=tt;
tend=tbeg+nt2-1;
dif=norm(PHI(:,1:nt2)-TEST(:,tbeg:tend),'fro');
if dif < best,
best=dif;
tbest=tt;
end
end
[best tbest]
%pause
tbeg=tbest;
tend=tbeg+nt-1;
color='red';

for i=1:3,
for j=1:4,
k=j+(i-1)*4;
yshift=4-i;
xshift=j-.8;
if k == 1,
sig=TEST(lab(9),tbeg:tend)-TEST(lab(8),tbeg:tend);
tekst=sprintf('%s','I');
elseif k == 2,
sig=TEST(lab(9),tbeg:tend)-(TEST(lab(8),tbeg:tend)+VFTEST(1,tbeg:tend))/2;
elseif k == 3,
sig=TEST(lab(1),tbeg:tend);
elseif k == 4,
sig=TEST(lab(5),tbeg:tend);
elseif k == 5,
sig=VFTEST(1,tbeg:tend)-TEST(lab(8),tbeg:tend);
elseif k == 6,
sig=TEST(lab(8),tbeg:tend);
elseif k == 7,
sig=TEST(lab(2),tbeg:tend);
elseif k == 8,
sig=TEST(lab(5),tbeg:tend);
elseif k == 9,
sig=VFTEST(1,tbeg:tend)-TEST(lab(9),tbeg:tend);
elseif k == 10,
sig=VFTEST(1,tbeg:tend);
elseif k == 11,
sig=.5*TEST(lab(3),tbeg:tend)+.5*TEST(lab(4),tbeg:tend);
elseif k == 12,
sig=TEST(lab(7),tbeg:tend);
end
h=plot(t+xshift,funscal*sig+yshift,color);
end
end
reldif=norm(PHI(:,1:nt)-TEST(:,tbeg:tend),'fro')/norm(PHI(:,1:nt),'fro');
set(text(3.5,-0.2,sprintf(' shifted: %0.5g  ms',-2*(tbest-25))),'color','k');
set(text(1.5,-0.2,sprintf(' reldiff: %0.5g',reldif)),'color','k');
set(text(3.5,0.2,sprintf('%s',file2)),'color','red');
if wct==1,
set(text(3.5,0,sprintf('%s','(wct)')),'color','k');
else
set(text(3.5,0,sprintf('%s','(zme)')),'color','k');
end
newfile='normals/r00015.qrswct';
newfile(11:14)=file2(3:6)
saveasci(newfile,TEST(:,tbeg:tend))
%pause
end