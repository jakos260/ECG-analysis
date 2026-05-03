% plotlds.m
% plot standard 12 leads
% in 'old' montage'
clear
clf
file='r01904.qst';
PHI=loadmat(file);
PHI=PHI';
PHI=zeromean(PHI);

color='blue';
n=size(PHI);
nlds=n(1);
ntims=n(2);
nt=ntims/2;
% 10 labels of signals to plot the 12-leads:
% V1,V2,(2 labels for V3), V4 to V6, VR VL
%lab=loadmat('nim12lds.labels');
lab=[19 26 33 34 41 48 54 1 2 56];

% establish reference
VFWCT=-(PHI(lab(8),:)+PHI(lab(9),:));
VF=VFWCT;
VFest=PHI(56,:);
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
sig=PHI(lab(9),tbeg:tend)-(PHI(lab(8),tbeg:tend)+VF(1,tbeg:tend))/2;
tekst=sprintf('%s','aVL');
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
sig=PHI(lab(8),tbeg:tend)-(PHI(lab(9),tbeg:tend)+VF(1,tbeg:tend))/2;
tekst=sprintf('%s','aVR');
elseif k == 7,
sig=PHI(lab(2),tbeg:tend);
tekst=sprintf('%s','V2');
elseif k == 8,
sig=PHI(lab(5),tbeg:tend);
tekst=sprintf('%s','V5');
elseif k == 9,
sig=VF(1,tbeg:tend)-PHI(lab(9),tbeg:tend);
tekst=sprintf('%s','III');
elseif k == 10,
sig=VF(1,tbeg:tend)-(PHI(lab(8),tbeg:tend)+PHI(lab(9),tbeg:tend))/2;
tekst=sprintf('%s','aVF');
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
set(text(0.,0.2,sprintf('%s',file)),'color','white');
xl=[1.2 1.2 2 ];
yl=[.2+funscal .2 .2];
set(line(xl,yl),'color','white')
text(1.25,.2+0.8*funscal,sprintf('%s','1 mV'));
text(1.8,.35,sprintf('%s','100 ms'));
if wct==1,
set(text(0,0,sprintf('%s','(wct)')),'color','white');
else
set(text(0,0,sprintf('%s','(zme)')),'color','white');
end
