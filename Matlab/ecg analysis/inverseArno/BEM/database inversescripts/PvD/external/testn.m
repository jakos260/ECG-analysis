ver1=[0 0 -2;-1 1 0; 1 1 0;0.5 0 2];
itri1=[1 2 3;2 4 3];

aa=[0:0.001:1]';
VM=ones(size(aa))*ver1(2,:)+aa*(ver1(3,:)-ver1(2,:));
a1=[];a2=[];a3=[];
for i=1:length(aa)
	a1=[a1 norm3d(ver1(1,:)-VM(i,:))];
	a2=[a2 norm3d(ver1(4,:)-VM(i,:))];
	a3=[a3 a1(end)+a2(end)];
end

figure(1)
clf
subplot(3,1,1)
plot(aa,a1,'r');
subplot(3,1,2)
plot(aa,a2,'g')
subplot(3,1,3)
plot(aa,a3,'b')
disp(num2str([aa(find(a3==min(a3))) a3(end)-a3(1) a3(1)<a3(end)]))
% a3(1)-a3(6)
% a3(end)-a3(1)
figure(2);clf	
a4=graphdist(itri1,ver1,4)
% colormap(loadmat('tims.mcm'))

patch('Faces',itri1,'Vertices',ver1,'FaceColor',[0.5 0.5 0.5],...
	  'edgecolor','k','buttondownFcn','selectnode','LineWidth',4);%[.99 .99 .99]
axis off equal tight; 

view(-23,0)
vm=VM(find(a3==min(a3)),:);
line([ver1(1,1) vm(:,1) ver1(4,1)],[ver1(1,2) vm(:,2) ver1(4,2)],[ver1(1,3) vm(:,3) ver1(4,3)],'color','r','linewidth',3,'linestyle',':');

foc=1;first=[2 3];sec=4;
line(vm(1),vm(2),vm(3),'Color',[0.99 0.99 0.99],'linestyle',':','linewidth',3,'Marker','.','linewidth',20,'Markersize',50) 

line(ver1(foc,1),ver1(foc,2),ver1(foc,3),'Color','y','linestyle','none','Marker','.','linewidth',20,'Markersize',50)
line(ver1(first,1),ver1(first,2),ver1(first,3),'Color','b','linestyle','none','Marker','.','linewidth',20,'Markersize',50)
line(ver1(sec,1),ver1(sec,2),ver1(sec,3),'Color','r','linestyle','none','Marker','.','linewidth',20,'Markersize',50) 


% line([ver1(4,1) vm(:,1)],[ver1(4,2) vm(:,2)],[ver1(4,3) vm(:,3)],'color','k','linewidth',3)