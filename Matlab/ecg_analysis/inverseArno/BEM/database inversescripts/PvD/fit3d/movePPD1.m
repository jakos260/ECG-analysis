clear; 


defname='C:\Geoms\geoms\geoms\default\PPD1';
outname='C:\Geoms\geoms\geoms\out\PPD2_655';
refname='C:\Geoms\geoms\geoms\test1.hver';
defrefname='C:\Geoms\geoms\geoms\default\PPD1_defVER.hver';

%%
Yref1=loadasci(refname);
X1=loadasci(defrefname); % X1 is the one to be moved

[AVER,AITRI]=loadtri([defname '_atria.tri']);
[LVER,LITRI]=loadtri([defname '_lendo.tri']);
[RVER,RITRI]=loadtri([defname '_rendo.tri']);
[VVER,VITRI]=loadtri([defname '_ventr.tri']);
[TVER,TITRI]=loadtri([defname '_torso.tri']);
[LLVER,LLITRI]=loadtri([defname '_llung.tri']);
[RLVER,RLITRI]=loadtri([defname '_rlung.tri']);
[HVER,HITRI]=loadtri([defname '_heart.tri']);

index=loadasci([defname '_heart.idx']);
saveindex([outname '_heart.idx'],index,'1) endo L atria, 2) endo R atria, 3) endo L ventricle, 4) endo R ventricle 5) epi atria, 6) epi  ventricles, 7) vessels');

X1(1,:)=[57.0889   46.1518  -46.9163]

apexi=find(HVER(:,1)==X1(1,1) & HVER(:,2)==X1(1,2) & HVER(:,3)==X1(1,3));


figure(10);clf
patch('Faces',HITRI,'Vertices',HVER,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','b','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
[V,I]=loadtri('C:\Geoms\geoms\geoms\test.tri');
hs=patch('Faces',I,'Vertices',V,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','r','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
view(90,0);axis off equal;
% return


% HVER=[HVER(:,2) -HVER(:,1) HVER(:,3)];

shift=Yref1(1,:)-HVER(apexi,:);
HVER=rotash(HVER,[0 0 0],shift);

tShift=Yref1(1,:);
HVER=rotash(HVER,[0 0 0],-tShift);

refframe='H:\ECG_simulation\measurements\PPD2\MRI\PPD2_2 110705\MR002674'
PPD1frame='H:/ECG_simulation/measurements/PPD1/dicom/wikr/heart (serie 11-13)/IM_01508'
A=dicominfo(refframe);
X2=A.ImageOrientationPatient(1:3)'; X2=[-X2(2) X2(1) X2(3)];
Y2=A.ImageOrientationPatient(4:end)'; Y2=[-Y2(2) Y2(1) Y2(3)];
Z2=cross(X2,Y2);
P2=A.ImagePositionPatient;

A=dicominfo(PPD1frame);
X1=A.ImageOrientationPatient(1:3)';X1=[-X1(2) X1(1) X1(3)];
Y1=A.ImageOrientationPatient(4:end)';Y1=[-Y1(2) Y1(1) Y1(3)];
Z1=cross(X1,Y1);
P1=A.ImagePositionPatient;
Z3=cross(Z1,Z2);	% rotation axis
alpha1=180*asin(norm3d(Z3))/pi;	% angle

plane1=[(100*X1-100*Y1);...
	    (100*X1+100*Y1);...
	    (-100*X1+100*Y1);...
	    (-100*X1-100*Y1)];

plane2=[(100*X2-100*Y2);...
	    (100*X2+100*Y2);...
	    (-100*X2+100*Y2);...
	    (-100*X2-100*Y2)];

ITRI=[1 2 3; 1 3 4];
figure(4);clf
hs=patch('Faces',HITRI,'Vertices',HVER,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','b','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
view(90,0);axis off equal;
patch('Faces',ITRI,'Vertices',plane1,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','b','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
patch('Faces',ITRI,'Vertices',plane2,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','m','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
% line([0 X2(:,1)],[0 X2(:,2)],[0 X2(:,3)],'Color','b','linestyle',':')
% line([0 Y2(:,1)],[0 Y2(:,2)],[0 Y2(:,3)],'Color','r','linestyle',':')
% line([0 Z2(:,1)],[0 Z2(:,2)],[0 Z2(:,3)],'Color','g','linestyle',':')


% line([0 X1(:,1)],[0 X1(:,2)],[0 X1(:,3)])
% line([0 Y1(:,1)],[0 Y1(:,2)],[0 Y1(:,3)],'Color','r')
% line([0 Z1(:,1)],[0 Z1(:,2)],[0 Z1(:,3)],'Color','g')
% 


X1=doRot(X1,Z3,alpha1);
Y1=doRot(Y1,Z3,alpha1);
Z1=doRot(Z1,Z3,alpha1);
% line([0 X1(:,1)],[0 X1(:,2)],[0 X1(:,3)])
% line([0 Y1(:,1)],[0 Y1(:,2)],[0 Y1(:,3)],'Color','r')
% line([0 Z1(:,1)],[0 Z1(:,2)],[0 Z1(:,3)],'Color','g')
plane1=doRot(plane1,Z3,alpha1);
HVER=doRot(HVER,Z3,alpha1);

[X2; Y2; Z2]-[X1; Y1; Z1]

patch('Faces',HITRI,'Vertices',HVER,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','r','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
patch('Faces',ITRI,'Vertices',plane1,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','r','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');

% clf
% line([0 pLR1(:,1)],[0 pLR1(:,2)],[0 pLR1(:,3)],'color','k')
% line([0 pLR2(:,1)],[0 pLR2(:,2)],[0 pLR2(:,3)],'Color','m')
% line([0 Z1(:,1)],[0 Z1(:,2)],[0 Z1(:,3)],'Color','c')
% line([0 Y3(:,1)],[0 Y3(:,2)],[0 Y3(:,3)],'Color','c')
pL1=X1*126 + Y1*126;pR1=X1* 38 + Y1*141 ;pLR1=pL1-pR1;pLR1=pLR1/norm3d(pLR1);
pL2=X2*152 + Y2*132;pR2=X2* 78 + Y2*108 ;pLR2=pL2-pR2;pLR2=pLR2/norm3d(pLR2);





Y3=cross(pLR1,pLR2);	% rotation axis
alpha2=asin(norm3d(Y3));	% angle
alpha2=180*alpha2/pi;


phi=180*acos(sum(pLR2.*pLR1))./pi %Z2 should be Z1 is second rotation axis
alpha2=phi;
X1=doRot(X1,Z1,alpha2);
Y1=doRot(Y1,Z1,alpha2);
Z1=doRot(Z1,Z1,alpha2);
HVER=doRot(HVER,Y3,alpha2);
plane1=doRot(plane1,Y3,alpha2);
a=170*Y1;
line([-a(1) a(1)],[-a(2) a(2)],[-a(3) a(3)])

patch('Faces',HITRI,'Vertices',HVER,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','g','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
patch('Faces',ITRI,'Vertices',plane1,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','g','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
% shift=Yref1(1,:)-HVER(apexi,:);
% HVER=rotash(HVER,[0 0 0],shift);

% patch('Faces',HITRI,'Vertices',HVER,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','m','edgecolor','k','FaceAlpha',1);

HVER=rotash(HVER,[0 0 0],tShift);

savetri([outname '_heart.tri'],HVER,HITRI);


return
PP=A.ImagePositionPatient';
% PP=[ 0 0 0];
PP=mean(Yref2);
X2=	  [(PP-100*Xp -100*Yp);...
	   (PP-100*Xp +100*Yp);...
	   (PP+100*Xp +100*Yp);...
	   (PP+100*Xp -100*Yp)];
a=X2(:,1);
X2(:,1)=-X2(:,2);
X2(:,2)=a;

a=Yref2(:,1);
Yref2(:,1)=-Yref2(:,2);
Yref2(:,2)=a;


a=X2;
X2=Yref2;
Yref2=a;
Yref2=Yref2+ones(size(Yref2,1),1)*mean(Yref2);
X2=X2+ones(size(X2,1),1)*mean(Yref2);

figure(3);clf
ITRI=[1 4 2;2 4 3;1 4 3;1 3 2];
patch('Faces',ITRI,'Vertices',X2,'FaceLighting','phong','FaceColor','g','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
patch('Faces',ITRI,'Vertices',Yref2,'FaceLighting','phong','FaceColor','b','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
patch('Faces',HITRI,'Vertices',HVER,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','r','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
axis equal off

%% scale the heart

scale=norm3d(Yref1(1,:)-Yref1(3,:))./norm3d(X1(1,:)-X1(3,:));

m=mean(HVER);
X1=(X1-ones(size(X1,1),1)*m).*scale+-ones(size(X1,1),1)*m;
AVER=(AVER-ones(size(AVER,1),1)*m).*scale+-ones(size(AVER,1),1)*m;
LVER=(LVER-ones(size(LVER,1),1)*m).*scale+-ones(size(LVER,1),1)*m;
RVER=(RVER-ones(size(RVER,1),1)*m).*scale+-ones(size(RVER,1),1)*m;
VVER=(VVER-ones(size(VVER,1),1)*m).*scale+-ones(size(VVER,1),1)*m;
HVER=(HVER-ones(size(HVER,1),1)*m).*scale+-ones(size(HVER,1),1)*m;


%% shift all geometries and default reference points
shift=(mean(Yref1)-mean(X1));
shift=Yref1(1,:)-X1(1,:);

X1=X1+ones(size(X1,1),1)*shift;
AVER=AVER+ones(size(AVER,1),1)*shift;
LVER=LVER+ones(size(LVER,1),1)*shift;
RVER=RVER+ones(size(RVER,1),1)*shift;
VVER=VVER+ones(size(VVER,1),1)*shift;
HVER=HVER+ones(size(HVER,1),1)*shift;

TVER=TVER+ones(size(TVER,1),1)*shift;
RLVER=RLVER+ones(size(RLVER,1),1)*shift;
LLVER=LLVER+ones(size(LLVER,1),1)*shift;

figure(4);

hs=patch('Faces',HITRI,'Vertices',HVER,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','m','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
view(90,0);axis off equal;


%%

shift=-Yref1(1,:)+mean(Yref2);

X1=X1+ones(size(X1,1),1)*shift;
AVER=AVER+ones(size(AVER,1),1)*shift;
LVER=LVER+ones(size(LVER,1),1)*shift;
RVER=RVER+ones(size(RVER,1),1)*shift;
VVER=VVER+ones(size(VVER,1),1)*shift;
HVER=HVER+ones(size(HVER,1),1)*shift;

TVER=TVER+ones(size(TVER,1),1)*shift;
RLVER=RLVER+ones(size(RLVER,1),1)*shift;
LLVER=LLVER+ones(size(LLVER,1),1)*shift;



% Yrefh=Yref2;Xh=X2;
% figure(1);clf;fitnodes;


shift=-shift;
X1=X1+ones(size(X1,1),1)*shift;
AVER=AVER+ones(size(AVER,1),1)*shift;
LVER=LVER+ones(size(LVER,1),1)*shift;
RVER=RVER+ones(size(RVER,1),1)*shift;
VVER=VVER+ones(size(VVER,1),1)*shift;
HVER=HVER+ones(size(HVER,1),1)*shift;

TVER=TVER+ones(size(TVER,1),1)*shift;
RLVER=RLVER+ones(size(RLVER,1),1)*shift;
LLVER=LLVER+ones(size(LLVER,1),1)*shift;


figure(4);

hs=patch('Faces',HITRI,'Vertices',HVER,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','g','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
view(90,0);axis off equal;



Yrefh=Yref1;Xh=X1;
figure(1);clf;fitnodes;

%%

savetri([outname '_atria.tri'],AVER,AITRI);
savetri([outname '_lendo.tri'],LVER,LITRI);
savetri([outname '_rendo.tri'],RVER,RITRI);
savetri([outname '_ventr.tri'],VVER,VITRI);
savetri([outname '_heart.tri'],HVER,HITRI);

savetri([outname '_torso.tri'],TVER,TITRI);
savetri([outname '_llung.tri'],LLVER,LLITRI);
savetri([outname '_rlung.tri'],RLVER,RLITRI);

hs=patch('Faces',HITRI,'Vertices',HVER,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','r','edgecolor','k','FaceAlpha',.1,'buttondownFcn','selectnode');
figure(4);

hs=patch('Faces',HITRI,'Vertices',HVER,'FaceLighting','phong','BackFaceLighting','lit','AmbientStrength',0.7,'FaceColor','r','edgecolor','k','FaceAlpha',1,'buttondownFcn','selectnode');
view(90,0);axis off equal;





