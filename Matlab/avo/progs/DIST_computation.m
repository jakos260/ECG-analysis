% DIST_computation.m
% A. van Oosterom; 2015_01_28
% creation/documenting/checking of ADJ matrices and distance matrices
% mode = 0 : unit edge lengths
% mode = 1 : actual edge distances (2D)
% mode = 2 : actual distances through interior (3D)
% mode = 3 : as in mode 2  (3D), but second order neighbours
%            are included
% mode = 4 : as in mode 1  (2D), but second order neighbours
%            are included

clear all
keep_edges=1;

geom='ARVC01_ventricles.tri';
[VER,ITRI]=loadtri(geom);
nv=size(VER,1);

%menu=[3:5]

menu=7

if ismember(1,menu),
    'start 1 atria'
    beep
    tic
    [ADJ0,DIST0]= graphdist(ITRI);
    tic
    savemat('ADJ0_atri.mat',ADJ0);
    time_run
    tic
    savemat('DIST0_atri.mat',DIST0);
    time_run
    
    ADJ0=[];
    DIST0=[];
    tic
    [ADJ1,DIST1]= graphdist(ITRI,VER,1);
    savemat('ADJ1_atri.mat',ADJ1);
    savemat('DIST1_atri.mat',DIST1);
    time_run
    
    ADJ1=[];
    DIST1=[];
    
    tic
    [ADJ2,DIST2]= graphdist(ITRI,VER,2);
    savemat('ADJ2_atria.mat',ADJ2);
    savemat('DIST2_atria.mat',DIST2);
    ADJ2=[];
    DIST2=[];
    time_run
end

if ismember(2,menu),
    'start 2; atria'
    beep
    %     ratio=2; %=velo_surf/velo_transm
    %     ADJ=ADJ2;
    %     ADJ(ADJ1==0)=ratio*ADJ(ADJ1==0);
    %     DIST3=graphdist(ADJ);
    ratio=2; %=velo_surf/velo_transm
    ADJ1=loadmat('ADJ1_atri.mat');
    ADJ2=loadmat('ADJ2_atri.mat');
    ADJ=ADJ2;
    ADJ(ADJ1==0)=ADJ(ADJ1==0)*ratio;
    tic
    DIST3=graphdist(ADJ);
    time_run
    savemat(DIST2_2_atri.may,DIST3);
end

if ismember(3,menu),
    'start 2; ventr'
    beep
    [VER,ITRI]=loadtri(geom);
    tic
    
    [ADJ0,DIST0]= graphdist(ITRI); % un_weighted surface connections
    tic
    savemat('ADJ0_ventr.mat',ADJ0);
    time_run
    tic
    savemat('DIST0_ventr.mat',DIST0);
    time_run
    
    ADJ0=[];
    DIST0=[];
    
    tic
    [ADJ1,DIST1]= graphdist(ITRI,VER,1); % weighted surface connections
    savemat('ADJ1_ventr.mat',ADJ1);
    savemat('DIST1_ventr.mat',DIST1);
    time_run
    
    ADJ1=[];
    DIST1=[];
    
    tic
    [ADJ2,DIST2]= graphdist(ITRI,VER,2); % weighted intramural connections
    savemat('ADJ2_ventr.mat',ADJ2);
    savemat('DIST2_1_ventr.mat',DIST2);
    ADJ2=[];
    DIST2=[];
    time_run % time run: 14 min 26.16 s
end

if ismember(4,menu),
    'start 4; ventr'
    beep
    ADJ1=loadmat('ADJ1_ventr.mat');
    ADJ2=loadmat('ADJ2_ventr.mat');
    ratio=2; %=velo_along_surfurface/velo_transmural
    ADJ=ADJ2;
    ADJ(ADJ1==0)=ADJ(ADJ1==0)*ratio;
    tic
    DIST=graphdist(ADJ);
    time_run % 9 min 8.4 s
    savemat('DIST2_2_ventr.mat',DIST);
end

if ismember(5,menu),
    'start 5;; display'
    
    VALS=loadmat('DIST2_2_ventr.mat');
    types=loadmat('ARVC01_ventricles.typ');
    
    %     edge1=loadmat('edgeMV.lst');
    %     edge2=loadmat('edgeTRIC.lst');
    %     edge3=loadmat('edgeRVOT.lst');
    
    cmap='tims.mcm';
    figure(1)
    center
    clf
    
    zebra=-10;
    triplot_contour
    hold on
    
    %     plot3(VER(edge1,1),VER(edge1,2),VER(edge1,3),'k')
    %     plot3(VER(edge2,1),VER(edge2,2),VER(edge2,3),'k')
    %     plot3(VER(edge3,1),VER(edge3,2),VER(edge3,3),'k')
    list1=find(types==1);
    plot3(VER(list1,1),VER(list1,2),VER(list1,3),'*k')
    list2=find(types==2);
    plot3(VER(list2,1),VER(list2,2),VER(list2,3),'*w')
    list3=find(types==3);
    plot3(VER(list3,1),VER(list3,2),VER(list3,3),'*m')
    
end

if ismember(6,menu),
    'start 6;; display'
    
    %geom='ARVC01_rcav.tri';
    geom='ARVC01_ventricles.tri';
    
    [VER,ITRI]=loadtri(geom);
    VALS=loadmat('DIST1_ventr.mat');
    %VALS=VER;
    %types=loadmat('ARVC01_rcav.typ');
    types=loadmat('ARVC01_ventricles.typ');
    
    %     edge1=loadmat('edgeMV.lst');
    %     edge2=loadmat('edgeTRIC.lst');
    %     edge3=loadmat('edgeRVOT.lst');
    
    cmap='tims.mcm';
    figure(1)
    center
    clf
    
    zebra=-10;
    triplot_contour
    hold on
    %     plot3(VER(edge1,1),VER(edge1,2),VER(edge1,3),'k')
    %     plot3(VER(edge2,1),VER(edge2,2),VER(edge2,3),'k')
    %     plot3(VER(edge3,1),VER(edge3,2),VER(edge3,3),'k')
    % list1=find(types==1);
    % plot3(VER(list1,1),VER(list1,2),VER(list1,3),'*k')
    % list2=find(types==2);
    % plot3(VER(list2,1),VER(list2,2),VER(list2,3),'*w')
    % list3=find(types==3);
    % plot3(VER(list3,1),VER(list3,2),VER(list3,3),'*m')
    
    
    
    route=getroute(VER,ITRI,385,393);
    route=[route getroute(VER,ITRI,393,401)];
    route=[route getroute(VER,ITRI,401,404)];
    route=[route getroute(VER,ITRI,404,601)];
    route=[route getroute(VER,ITRI,601,382)];
    route=[route getroute(VER,ITRI,382,385)];
    route([9 22 27 30])=[]
    
    form ='%4d\n'
    file='edgeMV.lst'
    saveasciform(file,route,form)
    edge=loadmat(file);
    plot3(VER(edge,1),VER(edge,2),VER(edge,3),'*-m')
    pause
    
    
    pause
    
    
    route=getroute(VER,ITRI,428,430);
    route=[route getroute(VER,ITRI,430,436)];
    route=[route getroute(VER,ITRI,436,424)];
    route=[route getroute(VER,ITRI,424,428)];
    route([4 12 26])=[]
    
    
    file='edgeTRIC.lst'
    saveasciform(file,route,form)
    edge=loadmat(file);
    plot3(VER(edge,1),VER(edge,2),VER(edge,3),'*-m')
    pause
    
    
    
    route=getroute(VER,ITRI,934,795);
    route=[route getroute(VER,ITRI,795,413)];
    route=[route getroute(VER,ITRI,413,934)];
    route([6 11])=[]
    file='edgeRVOT.lst'
    saveasciform(file,route,form)
    edge=loadmat(file);
    plot3(VER(edge,1),VER(edge,2),VER(edge,3),'*-m')
    pause
    
    plot3(VER(route,1),VER(route,2),VER(route,3),'*-m')
    
end

if ismember(7,menu),
    'check torso geom'
    [VER1,ITRI1]=loadtri('ARVC01_thorax.tri');
    [VER2,ITRI2]=loadtri('ARVC01_thoraxNijmegen_2015.tri');
end
    




