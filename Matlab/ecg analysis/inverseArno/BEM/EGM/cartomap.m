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
[ventricles.node, ventricles.face] = loadtri([fpath '/' subject '_model/' subject '_model_ventricles.tri']);

% load x-ray data manually:
[FileName, PathName]    = uigetfile('*','Select annotate output file.', [fpath '/Fluoroscopy'], 'MultiSelect', 'on');
cath_tip_ver            = zeros(length(FileName),1);
cath_tip_pnt            = zeros(length(FileName),3);
cath_rest_ver           = [];
cath_rest_pnt           = [];

for k = 1:length(FileName),
    load([PathName FileName{k}]);
    for d = 1:length(geom), selventricles(d) = strcmp(geom{d}.name, 'Ventricles'); end
    
    selcolumn           = find(selventricles == 1);
    C1                  = find(strcmp(mlabels, 'c1') == 1);
    C2                  = find(strcmp(mlabels, 'c2') == 1);
    C3                  = find(strcmp(mlabels, 'c3') == 1);
    C4                  = find(strcmp(mlabels, 'c4') == 1);
    cath_tip_ver(k)     = vernearest(C1, selcolumn);
    cath_rest_ver       = [cath_rest_ver; vernearest(C2, selcolumn); vernearest(C3, selcolumn); vernearest(C4, selcolumn)];
    cath_tip_pnt(k,:)   = coord(C1,:); %pnearest(C1, selcolumn,:);
    cath_rest_pnt       = [cath_rest_pnt; coord(C2,:); coord(C3,:); coord(C4,:)];
end

% Load CARTO points and latencies:
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
cartomap.node   = [-VER(:,1) VER(:,2) -VER(:,3)];
cartomap.field  = MAP;

% Save latencies:
cartolatencies.node     = cartomap.node;
cartolatencies.field    = LAT;

% Create structure for SCIRun of x-ray points:
xray.node       = [ventricles.node(cath_tip_ver,:); ventricles.node(cath_rest_ver,:)];
xray.field      = [8 5 6 7 1 1 1 1 1 1 1 1 1 1 1 1];

xraypnt.node    = [cath_tip_pnt; cath_rest_pnt];
xraypnt.field   = [8 5 6 7 1 1 1 1 1 1 1 1 1 1 1 1];

% save all data:
save([fpath '/CATHDATA/cartomap_scirun.mat'], 'cartomap', 'ventricles', 'cartolatencies', 'xray', 'xraypnt');