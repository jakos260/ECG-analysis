% makefronts_elgs.m
% 20051108
clear
a=2.6; b=4;


next=1;
if next==1;
    
ra=a;
rb=b;
nr=11;
nphi=61;

[VER,ITRI]=make_annulus(a,b,nr,nphi,0);
nver=size(VER,1);


% epicardial focus
pnt1=[0 b 0];
endo=0;

% endocardial focus
% pnt1=[0 a 0];
% endo=1;

for i=1:nver;
    VALS(i,1)=distance_annulus(pnt1,VER(i,:),a,b);
end

% ms if v=0.5m/s
% VALS=VALS*20;
 
grsw=0;

figure(1)
clf
cmap='tims.mcm';
zebra=-0.25;
lsw=0;
iview=9;
%triplot
%hold on


% identify isochrone: j
nr=11;

for i=1:nr,
   r(i)=b-(i-1)*(b-a)/(nr-1);
end

nphi=800;
for i=1:nphi+1,,
    phi(i)=(i-1)*pi/nphi;
    x(i)=sin(phi(i));
    y(i)=cos(phi(i));
    for j=1:nr,
        D(i,j)=distance_annulus(pnt1,r(j)*[x(i)+eps y(i),0],a,b);
    end
end

%figure(3)
%clf
%plot(D)

step=0.05;
figure(1)
d=step:step:max(max(D));
nd=length(d);

RES=zeros(nd,nr);
X=RES;
Y=RES;

for id=1:nd,
   for j=1:nr,
       if d(id)>=D(1,j) & d(id) < D(nphi+1,j),
          is=find(D(:,j) <=d(id));
          if isempty(is)==0, 
             ir=reverse(is);
             RES(id,j)=phi(ir(1));
             X(id,j)=r(j)*cos(RES(id,j));
             Y(id,j)=r(j)*sin(RES(id,j));
          end
       end
   end
end
%RES

% treat incomplete number of nodes identified on an isochrone; 
% force all to carry nr nodes

for id=1:nd,

    if RES(id,1)~=0 | RES(id,nr)~=0,
       rho=d(id);
       %['refine' num2str(id)]
       if endo==1,
          % endocardial focus 
          if RES(id,1)==0,
             alpha=atan2(r(nr)*sin(RES(id,nr))-pnt1(1),r(nr)*cos(RES(id,nr))-pnt1(2));
             beta=(0:nr-2)/(nr-1)*alpha;
             XY=[ cos(beta)*rho; sin(beta)*rho;];
             X(id,1:nr-1)=XY(1,:)+ones(1,nr-1)*pnt1(2);
             Y(id,1:nr-1)=XY(2,:)+ones(1,nr-1)*pnt1(1);
          end
       else,
          % epicardial focus
          if RES(id,nr)==0 & id < nd/2,,
             alpha=atan2(r(1)*sin(RES(id,1))-pnt1(1),-r(1)*cos(RES(id,1))+pnt1(2));
             beta=(nr-1:-1:0)/(nr-1)*alpha;
             XY=[-cos(beta)*rho; sin(beta)*rho;];
             X(id,1:nr)=XY(1,:)+ones(1,nr)*pnt1(2);
             Y(id,1:nr)=XY(2,:)+ones(1,nr)*pnt1(1);
          end
       end
       
       if RES(id,1)~=0 | RES(id,nr)~=0,
          rho=d(id);
          if endo==1,
             % endocardial focus 
             if RES(id,1)==0,
                alpha=atan2(r(nr)*sin(RES(id,nr))-pnt1(1),r(nr)*cos(RES(id,nr))-pnt1(2));
                beta=(0:nr-2)/(nr-1)*alpha;
                XY=[ cos(beta)*rho; sin(beta)*rho;];
                X(id,1:nr-1)=XY(1,:)+ones(1,nr-1)*pnt1(2);
                Y(id,1:nr-1)=XY(2,:)+ones(1,nr-1)*pnt1(1);
             end
          else,
             % epicardial focus
             if RES(id,nr)==0 & id < nd/2,,
                alpha=atan2(r(1)*sin(RES(id,1))-pnt1(1),-r(1)*cos(RES(id,1))+pnt1(2));
                beta=(nr-1:-1:0)/(nr-1)*alpha;
                XY=[-cos(beta)*rho; sin(beta)*rho;];
                X(id,1:nr)=XY(1,:)+ones(1,nr)*pnt1(2);
                Y(id,1:nr)=XY(2,:)+ones(1,nr)*pnt1(1);
             end
          end
       
          if RES(id,nr)==0 & id > nd/2,
             list=find(RES(id,:)==0);
             nl=list(1)-1;
             %'refine arc'
             Y(id,nl+1)=0;
             if nl>=2,
                aa=Y(id,nl-1)/Y(id,nl);
                X(id,nl+1)=(X(id,nl-1)-aa*X(id,nl))/(1-aa);
             end
             if nl==1,
                 X(id,nl+1)=-b+2*(X(id,nl)+b);
             end
             x=X(id,1:nl+1);
             xx=x(1)+(0:nr-1)/(nr-1)*(x(nl+1)-x(1));
             y=spline(x,Y(id,1:nl+1),xx);
             X(id,:)=xx;
             Y(id,:)=y;             
          end
       end
   end
end
 
nfronts=size(X,1);
for i=1:nfronts;
    plot3(-Y(i,:),X(i,:),zeros(1,nr),'W')
end

% create triangle specs wave front
nphi=24;
phi=(0:nphi-1)/nphi*2*pi;
ntri=0;
ITRI=[];
for j=1:nr-1;
    for i=1:nphi-1,
        ntri=ntri+1;
        ITRI(ntri,1:3)=[ (j-1)*nphi+i  (j-1)*nphi+i+1 j*nphi+i];
        ntri=ntri+1;
        ITRI(ntri,1:3)=[ (j-1)*nphi+i+1 j*nphi+i+1 j*nphi+i];
    end
    ntri=ntri+1;
    ITRI(ntri,1:3)=[ j*nphi (j-1)*nphi+1  (j+1)*nphi];
    ntri=ntri+1;
    ITRI(ntri,1:3)=[ (j-1)*nphi+1 j*nphi+1 (j+1)*nphi];
end
ITRIF=ITRI;

% create triangle specs cross-section myocard
ntri=0;
ITRI=[];

for j=1,
    for i=1:nphi-1,
        ntri=ntri+1;
        ITRI(ntri,1:3)=[ (j-1)*nphi+i  (j-1)*nphi+i+1 j*nphi+i];
        ntri=ntri+1;
        ITRI(ntri,1:3)=[ (j-1)*nphi+i+1 j*nphi+i+1 j*nphi+i];
    end
    ntri=ntri+1;
    ITRI(ntri,1:3)=[ j*nphi (j-1)*nphi+1  (j+1)*nphi];
    ntri=ntri+1;
    ITRI(ntri,1:3)=[ (j-1)*nphi+1 j*nphi+1 (j+1)*nphi];
end
ITRIC=ITRI;

figure(1)
hold on
iview=1;
OBS( 1:12,1:3)=[zeros(12,1) zeros(12,1) (2.2:0.2:4.4)'] ;
OBS(13:24,1:3)=[zeros(12,1) (2.2:0.2:4.4)'  zeros(12,1)] ;
OBS(25:36,1:3)=[zeros(12,1) zeros(12,1) -(2.2:0.2:4.4)'] ;
nobs=size(OBS,1);

nnphi=13;
phin=(0:nnphi-1)/(nnphi-1)*pi;
OBS(nobs+1:nobs+nnphi,1:3)=[zeros(nnphi,1) a*sin(phin') a*cos(phin') ];
nobs=size(OBS,1)

OBS(nobs+1:nobs+nnphi,1:3)=[zeros(nnphi,1) b*sin(phin') b*cos(phin') ];
nobs=size(OBS,1)

nobs=nobs+1;
OBS(nobs,1:3)=zeros(1,3);

   
for t=1:nfronts, 
    clf
    nver=0;
    VERF=[];
    for j=nr:-1:1,
        rho=Y(t,j); z=X(t,j);
        for k=1:nphi,
            nver=nver+1;
            VERF(nver,1:3)=[-rho*cos(phi(k)) rho*sin(phi(k)) z];
        end
    end
    
    edgea=[1:nphi 1];
    edgeb=[nver-nphi+1:nver nver-nphi+1];
   
    % make boundary contour
    nphi=24;
    phi=(0:nphi-1)/nphi*2*pi;
    Bb=[zeros(nphi+1,1) rb*[cos(phi) cos(phi(1))]' rb*[sin(phi) sin(phi(1))]'];
    Ba=ra/rb*Bb;
    VERC=[Bb(1:nphi,:); Ba(1:nphi,:)];
    [VER,ITRI]=addtwotris(VERC,ITRIC,VERF,ITRIF);
    VALS=[];
    lsw=1;
    zebra=[];
    boxsize=4;
    triplot
    view(90,-60)
    hold on
    plot3(Bb(:,1),Bb(:,2),Bb(:,3),'b')
    plot3(Ba(:,1),Ba(:,2),Ba(:,3),'r')
    
    set(ht,'vis','off')
    
    uitim=uicontrol('style','text');
    uiboxtim=[.18 .92 .15 .05];
    set(uitim,'units','norm','position',uiboxtim',...
    'string',['t= ' num2str(t*10*step)],'fontsize',10)
    
    %pause

    for kk=1:size(OBS,1),
        SA=dsa(VERF,ITRIF,OBS(kk,1:3),0.1);
        PHI(kk,t)=-10*sum(sum(SA))/pi;
    end
    
    plot3(OBS(1:nobs,1),OBS(1:nobs,2),OBS(1:nobs,3),'k+')
    pause(0.5)
end

nt=size(PHI,2);
nextra=210-nt;
PHI=[zeros(nobs,10) PHI zeros(nobs,nextra)];

nt=size(PHI,2);

for i=1:36,
    kshift(i)=20*icyc(i,12);
end

lshift=20*(1:13);

PHI(1:36,1:nt)=PHI(1:36,1:nt)+kshift'*ones(1,nt);
PHI(37:49,:)=PHI(37:49,:)+lshift'*ones(1,nt);
PHI(50:62,:)=PHI(50:62,:)+lshift'*ones(1,nt);
uiboxxas=[.05 .0 .85 .05];

figure(2)
clf
if endo==1,
   titel='endo sitimulus; transmural electrograms; stim site';
else
   titel='epi sitimulus; transmural electrograms; stim site';
end
hold on
for i=1:12,
   plot(PHI(i,:)')
   plot(PHI(i,:)-PHI(63,:),'r')
   plot([1 nt],[kshift(i) kshift(i)], 'k:');
   text(-10, kshift(i),num2str(i))
end
axis([0 220 0 280])
axis off
uixas2=uicontrol('style','text');
set(uixas2,'units','norm','position',uiboxxas',...
'string','vertical shift between signals: 20 mV; time bar: 100 ms','fontsize',8)
title(titel)
plot([10 110],[0 0],'k')
plot([10 10],[3 -3],'k')
plot([110 110],[3 -3],'k')

figure(3)
clf
if endo==1,
   titel='endo sitimulus; transmural electrograms; site: halfway';
else
   titel='epi sitimulus; transmural electrograms; site:halfway';
end
hold on
for j=13:24,
    i=icyc(j,12);
   plot(PHI(j,:)')
   plot(PHI(j,:)-PHI(63,:),'r')
   plot([1 nt],[kshift(i) kshift(i)], 'k:');
   text(-10, kshift(i),num2str(i))
end
axis([0 220 0 280])
axis off
uixas3=uicontrol('style','text');
set(uixas3,'units','norm','position',uiboxxas',...
'string','vertical shift between signals: 20 mV; time bar:: 100 ms','fontsize',8)
title(titel)

plot([10 110],[0 0],'k')
plot([10 10],[3 -3],'k')
plot([110 110],[3 -3],'k')


figure(4)
clf
if endo==1,
   titel='endo sitimulus; transmural electrograms; site: terminal';
else
   titel='epi sitimulus; transmural electrograms; site: terminal';
end
hold on
for j=25:36,
    i=icyc(j,12);
   plot(PHI(j,:)')
   plot(PHI(j,:)-PHI(63,:),'r')
   plot([1 nt],[kshift(i) kshift(i)], 'k:');
   text(-10, kshift(i),num2str(i))
end
axis([0 220 0 280])
axis off
uixas3=uicontrol('style','text');
set(uixas3,'units','norm','position',uiboxxas',...
'string','vertical shift between signals: 20 mV; time bar: 100 ms','fontsize',8)
title(titel)
plot([10 110],[0 0],'k')
plot([10 10],[3 -3],'k')
plot([110 110],[3 -3],'k')

figure(5)
clf
if endo==1,
   titel='endo sitimulus; electrograms; site: endo';
else
   titel='epi sitimulus; electrograms; site: endo';
end
hold on
for j=37:49,
    i=j-36;;
   plot(PHI(j,:)')
   plot(PHI(j,:)-PHI(63,:),'r')
   plot([1 nt],[lshift(i) lshift(i)], 'k:');
   text(-10, lshift(i),num2str(i))
end
axis([0 220 0 280])
axis off
uixas3=uicontrol('style','text');
set(uixas3,'units','norm','position',uiboxxas',...
'string','vertical shift between signals: 20 mV; time bar: 100 ms','fontsize',8)
title(titel)
plot([10 110],[0 0],'k')
plot([10 10],[3 -3],'k')
plot([110 110],[3 -3],'k')

figure(6)
clf
if endo==1,
   titel='endo sitimulus; electrograms; site: epi';
else
   titel='epi sitimulus; electrograms; site: epi';
end
hold on
for j=50:62,
    i=j-49;;
   plot(PHI(j,:)')
   plot(PHI(j,:)-PHI(63,:),'r')
   plot([1 nt],[lshift(i) lshift(i)], 'k:');
   text(-10, lshift(i),num2str(i))
end
axis([0 220 0 280])
axis off
uixas3=uicontrol('style','text');
set(uixas3,'units','norm','position',uiboxxas',...
'string','vertical shift between signals: 20 mV; time bar: 100 ms','fontsize',8)
title(titel)
plot([10 110],[0 0],'k')
plot([10 10],[3 -3],'k')
plot([110 110],[3 -3],'k')

figure(7)
clf
plot(PHI(nobs,:))

end



next=0;
if next==1,
   % compute geometry for treating boundary effects    
   figure(1)
   clf
   % specification of globe shaped boundary
   ntheta=13; nphib=2*(ntheta-1);
   [VER, ITRI]=makeglobe(ntheta,nphib); % unit sphere
   VERb=b*VER;
   VERa=a*VER;
   savetri('epicard.tri', VERb,ITRI);
   savetri('endocard.tri',VERa,ITRI);
   VER=VERa; 
   triplot
end

next=1;
if next==1,
    % compute matrices for treating boundary effects
    [VERa,ITRIa]=loadtri('endocard.tri');
    %[VERb,ITRIb]=loadtri('epicard.tri,VERb,ITRIb');
    nvera =length(VERa);
    ntria =length(ITRIa);

    nphia=1+sqrt(2*nvera-3);

siga=3; % cavity conductivity (relative to myocard and exterior)


% compute B matrix
% use spherical symmetry; 
% since B contains solid angles its values are 
% invariant under scaling of a.

B=zeros(nvera,nvera);
nhalf=nvera/2;
ref=nvera+1-[1:nvera];
sel=[1 2:nphia:nvera];
for i=1:nhalf,
    [B(i,1:nvera),jsing]=rowforw(VERa,ITRIa,VERa(i,:));
    B(nvera+1-i,ref)=B(i,:);
end

end

next=0;
if next==1,

% compute T matrix; used to compute potentials at the boundary from
% infinite medium potentials
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


end

