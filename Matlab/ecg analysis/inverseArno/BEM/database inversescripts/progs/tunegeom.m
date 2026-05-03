% file tunegeom.m
% date: 0408012
clear

[VER,ITRI1]=LOADTRI('hunterverti.dhk');
PNTS=loadmat('hunterverti.pntr');
DIRS=loadmat('hunterverti.fibr');

next=0;
if next==1,
    PNTS=loadmat('hunter.pntr');
    PNTS=PNTS/1000;
    exx= [min(VER(:,1)) max(VER(:,1))];
    exy= [min(VER(:,2)) max(VER(:,2))];   
    exz= [min(VER(:,3)) max(VER(:,3))];
    
    exax=[min(PNTS(:,2)) max(PNTS(:,2))];
    exay=[min(PNTS(:,3)) max(PNTS(:,3))];   
    exaz=[min(PNTS(:,4)) max(PNTS(:,4))];
    
    nsel=length(SEL);
    shift=[exx(1)-exax(1);
           exy(1)-exay(1);
           exz(1)-exaz(1)+.003;]
    SEL=PNTS(sel,2:4);
    SEL=SEL+ones(nsel,1)*shift';
    %savemat('hunterverti.pntr',SEL)
end

%PNTS=loadmat('externs.mat');

npnts=size(PNTS,1);
sel=1:npnts;
SEL=PNTS(sel,:);

figure(1)
clf
ITRI=ITRI1;
triplot
figure(2)
clf
crossec
level=.023
set(sllevel,'val',level)
set(sllevel,'max',0.035)
newcross

next=0;
if next==1,
del=0.0003;
slab=find(abs(SEL(:,3)-level)<=del);
    %identify all PNTS that are not contained in myocard
    nslab=size(slab);
   for i=1:nslab,
       sa=solida(VER,ITRI,SEL(slab(i),:));
       if abs(sum(sa,2))<.1,
           plot(SEL(slab(i),2),-SEL(slab(i),1),'rO');
       end
   end
rim=[682:702 682];
line(VER(rim,2),-VER(rim,1),VER(rim,3),'color','g')  
end


next=0;
 % adapt triangulation of the base
 if next==1,
    ITRI=ITRI1(1350:1400,:);
    VER=rotash(VER,[0 .5 0 ],[0 0 0 ]);
    figure(1)
    clf
    triplot
    pause
    for tri=1:51,
        [btri,BVERS]=buurtris(ITRI,tri)
        for i=1:length(btri),
            % try to improve the pair [tri btri(i)]
            % form quadrangle
            QUAD=[ITRI(tri,:); ITRI(btri(i),:)]
            clf
            triplot
            [angs,edge]=quadrangles(VER,QUAD);
            [edge(1:4);angs]
            [ma ima]=max(angs)
            hold on
            plot3(VER(edge(ima),1),VER(edge(ima),2),VER(edge(ima),3),'r*')
            if edge(ima)~=ITRI(tri,BVERS(i,1))&edge(ima)~=ITRI(tri,BVERS(i,2)),
               'adapt'
                pause
                ITRI(tri,1:3)=[edge(ima) edge(icyc(ima+2,4)) edge(ima+1) ];
                ITRI(btri(i),1:3)=[edge(icyc(ima+3,4)) edge(icyc(ima+2,4)) edge(ima)];
                clf
                triplot
            end
        pause
        break
     end
  end
end
    
    
    


      
