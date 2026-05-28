% view_geometry.m
% view triangulated geometry
% A. van Oosterom; 2013_10_01
clear

geom1='tors00.tri';

geom2='rlung.tri';
geom3='llung.tri';

geom4='atria.tri';
geom5='ventricles.tri';

geom6='rcavs.tri';
geom7='lcavs.tri';

[VER4,ITRI4]=loadtri(geom4);
[VER5,ITRI5]=loadtri(geom5);

[VER6,ITRI6]=loadtri(geom6);
[VER7,ITRI7]=loadtri(geom7);

nver4=size(VER4,1)
nver5=size(VER5,1)
nver6=size(VER6,1)
nver7=size(VER7,1)

edge1=loadmat('edgeTRIC.lst');
edge2=loadmat('edgeRVOT.lst');
edge3=loadmat('edgeMV.lst');
[types,extra]=loadmat('WK_ventricle_node.types');


next=1;
if next==1,
    % view atria; ventricles and cavities (incl aorta).
    figure(1)
    clf
    
    % [VERA,ITRIA]=addtwotris(VER5,ITRI5)
    %,VER4,ITRI4);
    %    [VERB,ITRIB]=addtwotris(VER6,ITRI6,VER7,ITRI7);
    %    [VER,ITRI]=addtwotris(VERA,ITRIA,VERB,ITRIB);
    %VER=VERA; ITRI=ITRIA;
    
    VER=VER5; ITRI=ITRI5;
    
    %[VER ITRI]=make_sphere(1);
    
    nver=size(VER,1);
    
    %edge=17:26;
    %seed=5;
    %         edge=edge3;  %endoLV
    %         seed=229;
    %
    edge=unique([edge1 edge2]) % endo RV
    seed=131;
    
    %   seed=855;
    %
    %    edge=unique([edge1 edge2 edge3]);
    
    center
    grsw=1;
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
    
    node=seed;
    setnode
    
    patchnodes=find_patch_2(VER,ITRI,edge,seed);
    
    ntri=size(ITRI,1);
    
    trispatch=find(sum(ismember(ITRI(:,[1 2 3]),patchnodes),2)==3*ones(ntri,1));
    
    
    
    ITRI=ITRI(trispatch,:);
    
    
    figure(1)
    clf
    triplot
    hold on
    
    plot3(VERE(edge1,1),VERE(edge1,2),VERE(edge1,3),'k','linewidth',1)
    view(-90, 0)
    plot3(VERE(edge2,1),VERE(edge2,2),VERE(edge2,3),'k','linewidth',1)
    plot3(VERE(edge3,1),VERE(edge3,2),VERE(edge3,3),'k','linewidth',1)
    
    
    
    
    
    
    plot3(VER(patchnodes,1),VER(patchnodes,2),VER(patchnodes,3),'*w')
    
    
end







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


next=0;
if next==1,
    % view lungs and myocard.
    figure(1)
    clf
    [VERA,ITRIA]=addtwotris(VER5,ITRI5,VER6,ITRI6);
    [VERB,ITRIB]=addtwotris(VER1,ITRI1,VER2,ITRI2);
    [VER,ITRI]=addtwotris(VERA,ITRIA,VERB,ITRIB);
    
    lsw=0;
    figure(1)
    clf
    VALS=VER;
    triplot
    %    ngeoms=2;
    %    kleur=['b' 'r' 'g'];
    %    figure(2)
    %    clf
    %    crossec
    pause
    
end

