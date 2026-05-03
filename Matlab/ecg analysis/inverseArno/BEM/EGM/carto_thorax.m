%% Script to convert data from CARTO and x-ray images into a format that
%% can be used in SCIRun. 
%
% Version 1.0: 5-jul-2015

clear all
close all
clc

% If directory name is not defined:
subject = 'NICE001'; fpath = ['/Users/arnojanssen/Documents/STW/PVCs/' subject];

% load CARTO data:
load([fpath '/CATHDATA/CARTO.mat'], 'EGM', 'CARTO');

% load mesh of ventricles:
[ventricles.node, ventricles.face]  = loadtri([fpath '/' subject '_model/' subject '_model_ventricles.tri']);
ventricles_typ                      = importdata([fpath '/' subject '_model/' subject '_model_ventricles.typ']);
ventricles_typ                      = ventricles_typ(2:end,1);

% Construct structure withendo and epi cardial surface nodes:
epi_node_nr                      	= find(ventricles_typ == 1);                                                % include all epicardial points
LV_node_nr                      	= find(ventricles_typ == 2 | ventricles_typ == 4 | ventricles_typ == 7);
RV_node_nr                       	= find(ventricles_typ == 3 | ventricles_typ == 5 | ventricles_typ == 6);

epi.node                            = ventricles.node(epi_node_nr,:);
LV.node                             = ventricles.node(LV_node_nr,:);
RV.node                             = ventricles.node(RV_node_nr,:);

% Find the triangles:
epi.face    = find(sum(ismember(ventricles.face,epi_node_nr)')' > 0);   
LV.face     = find(sum(ismember(ventricles.face,LV_node_nr)')' > 2);
RV.face     = find(sum(ismember(ventricles.face,RV_node_nr)')' > 2);

% load mesh of thorax:
[thorax.node, thorax.face] = loadtri([fpath '/' subject '_model/' subject '_model_thoraxAmsterdam.tri']);

%% load x-ray data manually:
[FileName, PathName]    = uigetfile('*','Select annotate output file.', [fpath '/Fluoroscopy'], 'MultiSelect', 'on');
cath_tip_ver            = zeros(length(FileName),1);
cath_tip_pnt            = zeros(length(FileName),3);
cath_rest_ver           = [];
cath_rest_pnt           = [];
xray_mark_label         = zeros(length(FileName),1);

for k = 1:length(FileName),
    load([PathName FileName{k}]);
    
    for m = 1:size(CARTO.markers,1), 
        msel = findstr(FileName{k},num2str(CARTO.markers(m,1)));
        if ~isempty(msel), ksel = [m,msel]; end
    end
    
    xray_mark_label(k) = str2num(FileName{k}(ksel(2):ksel(2)+size(num2str(CARTO.markers(ksel(1),1)),2)-1));
    
    manshift
    for d = 1:length(geom), selventricles(d) = strcmp(geom{d}.name, 'Ventricles'); end
    
    selcolumn           = find(selventricles == 1);
    C1                  = find(strcmp(mlabels, 'c1') == 1);
    C2                  = find(strcmp(mlabels, 'c2') == 1);
    C3                  = find(strcmp(mlabels, 'c3') == 1);
    C4                  = find(strcmp(mlabels, 'c4') == 1);
    cath_tip_ver(k)     = vernearest(C1, selcolumn);
    cath_rest_ver       = [cath_rest_ver; vernearest(C2, selcolumn); vernearest(C3, selcolumn); vernearest(C4, selcolumn)];
    cath_tip_pnt(k,:)   = coord(C1,:);
    cath_rest_pnt       = [cath_rest_pnt; coord(C2,:); coord(C3,:); coord(C4,:)];
end

% Create structure for SCIRun of x-ray points:
xray.node       = [ventricles.node(cath_tip_ver,:); ventricles.node(cath_rest_ver,:)];
xray.field      = [8 5 6 7 1 1 1 1 1 1 1 1 1 1 1 1];

xraypnt.node    = [cath_tip_pnt; cath_rest_pnt];
xraypnt.field   = [8 5 6 7 1 1 1 1 1 1 1 1 1 1 1 1];

%% File CARTO points and latencies:
for n = 1:size(CARTO.coordinates,1)/2,
    VER(n,:)    = CARTO.coordinates((n*2)-1,3:5);
    MAP(n)      = CARTO.points(n,2);
    LAT(n)      = CARTO.map_anno(n,7);
end

% Add markers for x-ray positions:
MAP(CARTO.markers(1,1)) = 5;
MAP(CARTO.markers(2,1)) = 6;
MAP(CARTO.markers(3,1)) = 7;
MAP(CARTO.markers(4,1)) = 8;

% Find 'bad' points:
delsel = find(CARTO.points(:,8) == -1);

% Remove points from VER:
VER(delsel,:)   = [];
MAP(delsel)     = [];
LAT(delsel)     = [];

% Create structure for SCIRun:
cartomap.node   = [VER(:,1) VER(:,2) VER(:,3); 0 0 0];
cartomap.field  = [MAP 10];

% Save latencies:
cartolatencies.node     = cartomap.node(1:end-1,:);
cartolatencies.field    = LAT;

% Coordinates of reference electrode of CARTO system and required rotation:
referencepoint  = [-104.9812 -26.0711 79.3922];
rotx            = 90;
roty            = 0;
rotz            = 90;

[newcartopoints.node]   = carto2model(cartomap.node, referencepoint, rotx, roty, rotz);
newcartopoints.field    = cartomap.field;

%% project the new carto points onto the ventricular mesh: 
proj_cartopoints        = [];
proj_cartofaces         = [];
proj_ver                = ventricles.node;
proj_fac                = ventricles.face;

for ncp = 1:size(newcartopoints.node,1)-1,
    
    % select the correct projection surface nodes & triangles:
    if newcartopoints.field(1,ncp) == 1 || newcartopoints.field(1,ncp) == 5 || newcartopoints.field(1,ncp) == 6 || newcartopoints.field(1,ncp) == 7,
        proj_nodes = RV.node; 
        proj_faces = proj_fac(RV.face,:);
    elseif newcartopoints.field(1,ncp) == 2,
        proj_nodes = [RV.node; LV.node; epi.node];
        proj_faces = [proj_fac(RV.face,:); proj_fac(LV.face,:); proj_fac(epi.face,:)];
    elseif newcartopoints.field(1,ncp) == 3 || newcartopoints.field(1,ncp) == 4 || newcartopoints.field(1,ncp) == 8,
        proj_nodes = [LV.node; epi.node];
        proj_faces = [proj_fac(LV.face,:); proj_fac(epi.face,:)];
    end
    
    % project on vertices:
    tdist = zeros(size(proj_nodes,1),1);
    
    for prp = 1:size(proj_nodes,1), tdist(prp) = norm(newcartopoints.node(ncp,:) - proj_nodes(prp,:)); end
    
    [prpdist, prpindex] = min(tdist);
    
    proj_cartopoints.node(ncp,:)    = proj_nodes(prpindex,:);
    proj_cartopoints.field(ncp)     = prpdist;
    
    % project on faces:
    [la, mu, dist]                  = ptridist(newcartopoints.node(ncp,:), proj_ver, proj_faces);
    [mind, imind]                   = min(abs(dist));
    lai                             = la(imind);
    mui                             = mu(imind);
    
    proj_cartofaces.node(ncp,:)     = (1 - lai - mui)*proj_ver(proj_faces(imind,1),:) + lai*proj_ver(proj_faces(imind,2),:) + mui*proj_ver(proj_faces(imind,3), :);
    proj_cartofaces.field(ncp)      = mind;
    proj_cartofaces.trinearest(ncp) = imind;
end

%% save all data:
save([fpath '/CATHDATA/cartothorax_scirun.mat'], 'cartomap', 'ventricles', 'cartolatencies', 'thorax', 'newcartopoints', 'xray', 'xraypnt', 'proj_cartopoints', 'proj_cartofaces');