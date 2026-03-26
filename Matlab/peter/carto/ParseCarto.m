% -- Script to parse points and unipolar voltage recordings
%    Carto .car file -- %

% function ParseCarto(path2use,

clear all
close all

origpath = pwd;
path2use = 'C:\Users\damp2\Documents\Data\measurements\UCLA\prosCIPS\proCIPS02\carto\Patient 2014_08_18\Study 1\Export';
% path2use = uigetdir('','Select the directory with .car file');
sprintf('\n The Selected Directory with .car file: %s\n',path2use)

%disp('Format for Input of Directory/File Name is the following: ''FileName''');
outputdatadir = 'Parsed_CarFile';
outfilename_pts = 'Points';
outfilename_fld = 'Field1';
outfilename_Upot = 'Unipolar';
outfilename_Bpot = 'Bipolar';
outfilename_LAT = 'LAT';
outfilename_ABL = 'ABL';
cd(path2use)

h1 = dir(fullfile(path2use,'*car*'));

h1=h1(4);

Len_h1 = length(h1);
h2 = dir(outputdatadir);
if length(h2) == 0
    mkdir(outputdatadir);
    h2 = dir(outputdatadir);
end;
h3 = dir(fullfile(path2use,'*.vsl'));

if Len_h1 > 1
    fname = input('More than one .car file.  Please enter the desired file. Example("blauer3mo.car")', 's');
else
    fname = h1.name;    
end
if ~isempty(strfind(fname,'car'))
    fid = fopen(fname);
    % C1-Point# C2-C4 X,Y,Z Coords C5-Unipolar C6-Bipolar C7-LAT
    CAR = textscan(fid, '%*c %*n %n %*n %n %n %n %n %n %n %n %n %n %*n %*n %n %*[^\n]');
    %CAR = textscan(fid, '%*c %*n %n %*n %n %n %n %n %n %n %n %n %n');
    fclose(fid);
end

if ~isempty(strfind(fname,'vsl'))
    VSL=1;
end
% For ABL, 1 = ablation pnt and 0 != ablation point
carfields = {'Point','Xcoord','Ycoord','Zcoord','Xdir','Ydir','Zdir','Unipolar','Bipolar','LAT','ABL'};
CAR = cell2struct(CAR,carfields,2);

maxpntnum = CAR.Point(end,1);
Len_pts = length(CAR.Point);

%Bigsort = zeros(maxpntnum,7);
Bigsort = nan(maxpntnum,10);
Bigsort(:,1) = 1:maxpntnum;
k = 0;
for k = 1:Len_pts
    node = CAR.Point(k,1);
    Bigsort(node,2) = CAR.Xcoord(k,1);
    Bigsort(node,3) = CAR.Ycoord(k,1);
    Bigsort(node,4) = CAR.Zcoord(k,1);
    Bigsort(node,5) = CAR.Xdir(k,1);
    Bigsort(node,6) = CAR.Ydir(k,1);
    Bigsort(node,7) = CAR.Zdir(k,1);
    Bigsort(node,8) = CAR.Unipolar(k,1);
    Bigsort(node,9) = CAR.Bipolar(k,1);
    Bigsort(node,10) = CAR.LAT(k,1);
    Bigsort(node,11) = CAR.ABL(k,1);
end

cd(outputdatadir);

Points(:,1) = Bigsort(:,2);
Points(:,2) = Bigsort(:,3);
Points(:,3) = Bigsort(:,4);
Vectors(:,1) = Bigsort(:,5);
Vectors(:,2) = Bigsort(:,6);
Vectors(:,3) = Bigsort(:,7);

Field1.node = Points;
Field1.field = Vectors;
save(outfilename_fld,'Field1');
dlmwrite(outfilename_pts,Points,'delimiter',' ','newline','pc');
Unipolar = Bigsort(:,8);
save(outfilename_Upot, 'Unipolar');
Bipolar = Bigsort(:,9);
save(outfilename_Bpot, 'Bipolar');
LAT = Bigsort(:,10);
save(outfilename_LAT, 'LAT');
ABL = Bigsort(:,11);
save(outfilename_ABL, 'ABL');
