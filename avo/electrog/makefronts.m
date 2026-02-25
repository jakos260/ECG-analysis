% makefronts.m
% 20051119
clear
a=2.6; b=4;
ra=a;
rb=b;

nr=11;
nphi=81;

[VER,ITRI]=make_annulus(a,b,nr,nphi,0);
nver=size(VER,1);

[min(norm3d(VER)) max(norm3d(VER))]

% epicardial focus
pnt1=[0 b 0];
endo=0;

% endocardial focus
pnt1=[0 a 0];
endo=1;

endo

for i=1:nver;
    VALS(i,1)=distance_annulus(pnt1,VER(i,:),a,b);
end

 
grsw=0;

figure(1)
clf
cmap='tims.mcm';
zebra=-0.25;
lsw=0;
iview=9;
triplot
hold on


% identify isochrone: j


for i=1:nr,
   r(i)=b-(i-1)*(b-a)/(nr-1);
end

nphi=1200;
for i=1:nphi+1,,
    phi(i)=(i-1)*pi/nphi;
    x(i)=sin(phi(i));
    y(i)=cos(phi(i));
    for j=1:nr,
        D(i,j)=distance_annulus(pnt1,r(j)*[x(i)+eps y(i) 0],a,b);        
    end
end

%figure(3)
%clf
%plot(D)


figure(1)

d=.1:0.1:max(max(D));
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

% treat incomplete number of nodes identified on an isochrone; all are
% desired to show nr nodes
    
for id=80:nd,
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
         if RES(id,nr)==0 & id > nd/2,
            list=find(RES(id,:)==0);
            nl=list(1)-1;
            %'refine arc'
            Y(id,nl+1)=0;
            if nl>=2,
               a=Y(id,nl-1)/Y(id,nl);
               X(id,nl+1)=(X(id,nl-1)-a*X(id,nl))/(1-a);
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
   
%    plot3(-Y(id,:),X(id,:),zeros(1,nr),'w+')
   uitim=uicontrol('style','text');
   uiboxtim=[.18 .92 .15 .05];
   set(uitim,'units','norm','position',uiboxtim',...
   'string',['t= ' num2str(id)],'fontsize',10)
%    pause
end
 
