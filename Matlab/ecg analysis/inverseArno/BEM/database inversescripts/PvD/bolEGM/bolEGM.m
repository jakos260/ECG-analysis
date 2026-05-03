% 
% fn='bol4cm_2mm';
% 
% % [VER,ITRI]=loadtri([fn '.tri']);
% % DIS=loadmat([fn '.dis3d']);
% 
% col=DIS(1,:);
% clf
% colormap(loadmat('TIMS.mcm'));
% patch('Faces',ITRI,'Vertices',VER,'FaceLighting','phong',...
% 		'BackFaceLighting','lit','AmbientStrength',0.7,...
% 		'FaceVertexCData',col','FaceColor','interp','edgecolor','none','ButtondownFcn','SelectNode');
% 
% colorbar






r=0.1:0.01:20;  % uitspreidende golf-front in de tijd
phi=zeros(size(r));
for i=1:length(phi)
	r1=r(i);   % ring 1 [mm]
	r2=r(i)+1; % ring 2 [mm]
    
   x1=r1:-r1/steps:-r1;
   y1=sqrt(r1^2 -x1.^2); 
	y1=[y1 -y1(2:length(y1)-1)];
	xs=sort(x1(1:length(x1)-1));
	x1=[x1 xs(1:length(xs)-1)];	
	
	x2=r2:-r2/steps:-r2;
    y2=sqrt(r2^2 -x2.^2); 
	y2=[y2 -y2(2:length(y2)-1)];	
	xs=sort(x2(1:length(x2)-1));
	x2=[x2 xs(1:length(xs)-1)];
	
	l1=[x1; y1;    zeros(1,length(x1))]';   % ring 1
	l2=[x2; y2; 2.*ones(1,length(x2))]';    % ring 2 2 mm onder of boven
	ver=[l1;l2];
    
	np=length(x1);
	itri1=[1:np ; [np+2:2*np ,np+1]; np+1:2*np]';
	itri2=[1:np ; [2:np,1] ; [np+2:2*np ,np+1]]';
	itri=[itri1;itri2];

    % V = Vd/4pi * 0.5 * Omega
    phi(i)=40/(8*pi)*sum(solida(ver,itri,p));    
end
































% return

dt=0.1;
p=[4 0.001 -0.15 0.09 0];
S=gets_v([1:dt:240],10,190,p,3);
dS=diffrows(S);


dk=10;
az=0;
len=100;
VER=[-len 0 0; len 0 0; len dk 0; -len dk 0];
ITRI=[1 2 3; 3 4 1];
tz=-1e-5:dt:200;
templ=zeros(size(tz));
for i=1:length(tz)
	z=tz(i);
	VER=[-len 0 z; len 0 z; len dk z+az; -len dk z+az];
	templ(i)=sum(solida(VER,ITRI,[0 0 0]));
end
templ=10*templ/max(templ);

EGM=conv(templ,dS);
t=[0:dt:(length(S)-1)*dt];
clf
subplot(4,1,1)
plot([0:dt:(length(S)-1)*dt],S,'k','linewidth',2)
axis tight
subplot(4,1,2)
plot([0:dt:(length(dS)-1)*dt],dS,'r','linewidth',2);
subplot(4,1,3)
plot([0:length(templ)-1]/length(templ)-0.5,templ,'r','linewidth',2);
axis tight
grid
subplot(4,1,4)
plot([0:dt:(length(EGM)-1)*dt],EGM,'b','linewidth',2);
axis tight
grid

% hold on
% plot(S,'r')

