% remesh_ventricles.m

% A. van Oosterom; 2013_10_06
clear 

geom1='tors00.tri'
geom2='rlung.tri';
geom3='llung.tri';
geom4='atria.tri';
geom5='ventricles.tri';
geom6='rcavs.tri';
geom7='lcavs.tri';

%[VER4,ITRI4]=loadtri(geom4);
[VER5,ITRI5]=loadtri(geom5);
% [VER6,ITRI6]=loadtri(geom6);
% [VER7,ITRI7]=loadtri(geom7);
% nver4=size(VER4,1)
% nver5=size(VER5,1)
% nver6=size(VER6,1)
% nver7=size(VER7,1)

edge1=loadmat('edgeTRIC.lst');
edge2=loadmat('edgeRVOT.lst');
edge3=loadmat('edgeMV.lst');
[types,extra]=loadmat('WK_ventricle_node.types');


figure(1)
clf

% [VERA,ITRIA]=addtwotris(VER5,ITRI5)
%,VER4,ITRI4);
%    [VERB,ITRIB]=addtwotris(VER6,ITRI6,VER7,ITRI7);
%    [VER,ITRI]=addtwotris(VERA,ITRIA,VERB,ITRIB);
%VER=VERA; ITRI=ITRIA;

VER=VER5; ITRI=ITRI5;
nver=size(VER,1);


%         edge=edge3;  %endoLV
%         seed=229;


%center
%grsw=1;

ngeoms=1;
lsw=0;
figure(1)
clf
VALS=VER;
VERE=VER;
triplot
hold on
%plot3(VER(edge,1),VER(edge,2),VER(edge,3),'k','linewidth',1.5)

%     fac=1.05
%    for i=1:nver,
%        text(fac*VER(i,1),fac*VER(i,2),fac*VER(i,3),num2str(i))
%    end


VERE=VER;
hold on
plot3(VERE(edge1,1),VERE(edge1,2),VERE(edge1,3),'k','linewidth',1)
view(-90, 0)
plot3(VERE(edge2,1),VERE(edge2,2),VERE(edge2,3),'k','linewidth',1)
plot3(VERE(edge3,1),VERE(edge3,2),VERE(edge3,3),'k','linewidth',1)

%    nodes_0=find(types==0); % epiLV
%    nodes_1=find(types==1); % epiRV
%    nodes_2=find(types==2); % septum RV
%    nodes_3=find(types==3); % free wall RV
%    nodes_4=find(types==4); % free wall LV
%    nodes_5=find(types==5); % septum LV

%    plot3(VER(nodes_0,1),VER(nodes_0,2),VER(nodes_0,3),'w*')
%    plot3(VER(nodes_1,1),VER(nodes_1,2),VER(nodes_1,3),'y*')
%    plot3(VER(nodes_2,1),VER(nodes_2,2),VER(nodes_2,3),'m*')
%    plot3(VER(nodes_3,1),VER(nodes_3,2),VER(nodes_3,3),'r*')
%    plot3(VER(nodes_4,1),VER(nodes_4,2),VER(nodes_4,3),'g*')
%    plot3(VER(nodes_5,1),VER(nodes_5,2),VER(nodes_5,3),'b*')

% node=seed;
% setnode


edge=unique([edge1 edge2 edge3]) % epi
seed=855;

[ADJ,DIST]=graphdist(ITRI,VER,1);

patchnodes=find_patch_2(VER,ITRI,edge,seed);

ntri=size(ITRI,1);
trispatch=find(sum(ismember(ITRI(:,[1 2 3]),patchnodes),2)==3*ones(ntri,1));
ITRI=ITRI(trispatch,:);

%nnodes=size(unique(ITRI),1)-size(unique(edge),2)

% 851-77=774


figure(1)
clf
%center
% zebra=-10;
%stripes=[3.5 8 12  20:10:130 135];
stripes=[ 7 12  18:7.5:130.5 135];

for i=1:1500,
VALS(i,1)=min(DIST(i,edge));
end
triplot_contour;
hold on

plot3(VERE(edge1,1),VERE(edge1,2),VERE(edge1,3),'k','linewidth',1)
view(-90, 0)
plot3(VERE(edge2,1),VERE(edge2,2),VERE(edge2,3),'k','linewidth',1)
plot3(VERE(edge3,1),VERE(edge3,2),VERE(edge3,3),'k','linewidth',1)


nstripes=size(stripes,2);
line_up_contour_segments;

% contoursegments are closed

nnew=zeros(nstripes,1);
for i=1:nstripes,
    k=find(STRINGS(:,4)==i);
    STRING=STRINGS(k,1:3);
    ns=size(STRING,1);
    plot3(STRING(:,1),STRING(:,2),STRING(:,3),'k')
    l=sum(norm3d(STRING(2:end,:)-STRING(1:end-1,:))); % length of contour i
    nnew_i=round(l/6.75); % 6.75 mm
    if l<50,
        nnew_i=5;
    end
    nnew(i)=nnew_i;
end

% tune nnew
nnew(1)=77;
nnew(3:12)=nnew(3:12)-ones(10,1);
nnew(4:11)=nnew(4:11)-ones(8,1);
nnew
sum(nnew)

% define new nodes
NEWVER=[];
for i=1:nstripes,
    k=find(STRINGS(:,4)==i);
    CONT=STRINGS(k,1:3);
    NEWC=fourier_contours(CONT,round(nnew(i)/2)-1,2,nnew(i)); % end does not replicate begin
    nc=size(NEWC,1);
    NEWVER=[NEWVER; [NEWC i*ones(nc,1)]];
    figure(1)
    plot3(NEWC(:,1),NEWC(:,2),NEWC(:,3),'*k','linewidth',1)
    pause
end




for i=1:nstripes,
   
    NEWVER=[NEWVER; [NEWC i*ones(nc,1)]];
    figure(1)
    plot3(NEWC(:,1),NEWC(:,2),NEWC(:,3),'*k','linewidth',1)
    pause
end







%[ITRI,lista,listb]=make_peel(VER,lista,listb,mode,fig)



%
%    kleur=['b' 'r' 'g'];
%    figure(2)
%    clf
%    zincr=1;
%    zlevel=39;
%    delslab=0.5;
%    crossec
%    pause
% end


