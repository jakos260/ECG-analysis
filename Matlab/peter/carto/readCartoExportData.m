function D=readCartoExportData(pathname,xmlFilename)


D=[];

outpathname=fullfile(pathname,'data');
if ~exist(outpathname,'dir')
    mkdir(outpathname)
end

XMLIN=xml2struct(fullfile(pathname,xmlFilename));
for i=1:length(XMLIN.Study.Maps.TagsTable.Tag)
    DATA.TagTable.Colors(i,1:4) = str2num(XMLIN.Study.Maps.TagsTable.Tag{i}.Attributes.Color);
    DATA.TagTable.Name{i} = XMLIN.Study.Maps.TagsTable.Tag{i}.Attributes.Full_Name;
    DATA.TagTable.Id(i) = str2num(XMLIN.Study.Maps.TagsTable.Tag{i}.Attributes.ID);
end
DATA.ModelView = reshape(str2num(XMLIN.Study.Enviroment.Camera.ModelViewMatrix.Text),4,4);
DATA.ProjectionMatrix = reshape(str2num(XMLIN.Study.Enviroment.Camera.ProjectionMatrix.Text),4,4);
DATA.Camera.Aspectratio = XMLIN.Study.Enviroment.Camera.Attributes.AspectRatio;
DATA.Camera.Center= XMLIN.Study.Enviroment.Camera.Attributes.Center;
DATA.Camera.Scale = XMLIN.Study.Enviroment.Camera.Attributes.Scale;

DATAALL = DATA;


myfiles=dir(fullfile(pathname, '*_car.txt'));
k=1;
for imap=1:length(myfiles)    
    DATA = DATAALL;
    if isfield(XMLIN.Study.Maps.Map{imap}.CartoPoints,'Point')
        for i=1:length(XMLIN.Study.Maps.Map{imap}.CartoPoints.Point)
            if isfield(XMLIN.Study.Maps.Map{imap}.CartoPoints.Point{i},'Tags')
                DATA.TagPoints.Id(i) = str2num(XMLIN.Study.Maps.Map{imap}.CartoPoints.Point{i}.Tags.Text);
            else
                DATA.TagPoints.Id(i) = -1;
            end
            DATA.TagPoints.Pos3D(i,1:3) = str2num(XMLIN.Study.Maps.Map{imap}.CartoPoints.Point{i}.Attributes.Position3D);
            DATA.TagPoints.cathOrient(i,1:3) = str2num(XMLIN.Study.Maps.Map{imap}.CartoPoints.Point{i}.Attributes.CathOrientation);
        end
    end
    DATA.CAR = parseCartoTxt(pathname, myfiles(imap).name);
    if ~isempty(DATA)
        D{k} =DATA;
        save(fullfile(outpathname,[myfiles(imap).name(1:end-3) 'mat']),'DATA');
        k=k+1;
    end
end


%%

function DATA = parseCartoTxt(pathname, filename)

mydir=pwd;
meshname = [filename(1:end-8) '.mesh'];
[DATA.VER,DATA.ITRI]=meshparse(pathname,meshname);
cd(mydir);
if isempty(DATA.VER)
    DATA=[];
    return;
end
fid = fopen(fullfile(pathname,filename));
if fid>0 
    str = fgetl(fid);
    %                             1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
    CAR = textscan(fid, '%c %d %d %d %f %f %f %f %f %f %f %f %d %d %d %d %d %c %d %d %d');
    fclose(fid);
else
    warning(['could not open ' filename] );
    DATA =[];
    return 
end
DATA.name = filename(1:end-8);
DATA.ipoint = cell2mat(CAR(3));
DATA.points = cell2mat(CAR(5:7));
DATA.data = [cell2mat(CAR(8:12)) double(cell2mat(CAR(13:17)))];
DATA.annotated = cell2mat(CAR(18));

keep=[];
use = zeros(length(DATA.VER(:,1)),1);
for i=1:length(DATA.points(:,1))
    d = norm3d([DATA.points(i,1) - DATA.VER(:,1) DATA.points(i,2) - DATA.VER(:,2) DATA.points(i,3) - DATA.VER(:,3)]);
    
    di= find(d==min(d(use==0)));di=di(1);
    use(di)=1;
    keep = [keep; [di i]];
end
if ~isempty(keep)
    DATA.mapping = keep;
    DATA.T=intripol(DATA.VER,DATA.ITRI,DATA.mapping(:,1));

    for i=1:length(DATA.ipoint)
        [ECG,DATA.channels] = parseCartoECG(fullfile(pathname,[ filename(1:end-8) '_P' num2str(DATA.ipoint(i)) '_ECG_Export.txt']));
        DATA.ECG_all{i} = ECG;
        DATA.ECGM1(i,:) = ECG(1,:);
        DATA.ECGM2(i,:) = ECG(2,:);    
        DATA.ECGM3(i,:) = ECG(2,:);    
        DATA.ECGM4(i,:) = ECG(2,:);    
    end
end

%%
function [ECG, channels]= parseCartoECG(filename)

fid=fopen(filename);

str=fgetl(fid);
str=fgetl(fid);
if strfind(str,'Raw ECG to MV (gain) =')
    gain= str2num(str(length('Raw ECG to MV (gain) = '):end));
end
str=fgetl(fid);
str=fgetl(fid);
channels=[];
i=1;

parseString=[];
while ~isempty(str)
    [channel,str]=strtok(str,' ');
    if ~isempty(channel)
        channels{i} = channel;
        parseString=[parseString '%d ']; 
        i=i+1;
    end
end
CAR = textscan(fid,parseString );

ECG= double(cell2mat(CAR(1:end)))*gain;
ECG= ECG';

fclose(fid);

%%


function [VER,ITRI]=meshparse(path2use,filename)
% clear all
% close all

origpath = pwd;

% path2use = uigetdir('','Select the directory with .mesh file');
sprintf('\n The Selected Directory with .mesh file: %s\n',path2use)

% h1 = dir(fullfile(path2use,'LAtrim.mesh'));
% h1 = dir(fullfile(path2use,'*.mesh'));
%h1 = dir(fullfile(path2use,'*.ashx'));

VER=[];
ITRI=[];

% cd(path2use)

outfilename_Upot = 'Unipolar';
outfilename_Bpot = 'Bipolar';
outfilename_LAT = 'LAT';
outfilename_fld = 'Field1';
outputdatadir = 'Parsed_MeshFile';
% h2 = dir(['../',outputdatadir]);
% if isempty(h2) == 1
%     mkdir(['../',outputdatadir]);
%     h2 = dir(['../',outputdatadir]);
% end

% if Len_h1 > 1
%     fname = input('More than one .mesh file.  Please enter the desired file. Example("blauer3mo.mesh")', 's');
% else
%     fname = h1.name;
% end

% for iname =1:length(h1)
fname = fullfile(path2use,filename);%  h1(iname).name;
% outfname = fname(1:end-5);
% cd(path2use)

if ~isempty(strfind(fname,'mesh'))
    fid = fopen(fname);
    if fid<=0
        warning('cannot read mesh file');
        return;
    end
    
    mapmesh = textscan(fid,'%s');
    fclose(fid);
    mapmesh = mapmesh{:};
end

if -isempty(strfind(fname,'mesh'))
    mapmesh=1;
end

normcells = strmatch('NormalZ',mapmesh);
%finds the normal cells associated with the VerticesSection and the
%TrianglesSection data

datacells12 = (normcells+4);
%identifies the cells where the data begins for the first two sections

vertcells = strmatch('[Vertices',mapmesh);
%the second and third points found are points of reference for the second
%two sections of data

vert1 = vertcells(2)+27;
%provides the start point for the Vertices Color Section

%vert2 = vertcells(3)+37;
%provides the start point for the Vertices Attributes Section

s = length(mapmesh)+1;
%provides the end point for the Vertices Attribute Section


%%Isolation of Vertices Section%%

data1 = (datacells12(1,1):9:datacells12(2,1))';
%provides the start point for each row of data in the Vertices Section

L = length(data1)-2;
%provides the correct number of iterations for the outer for loop

for k = 1:L
    
    for m = 1:6
        %there are six iterations becauset there are six columns of data
        
        norm(k,m) = str2num(mapmesh{data1(k,1)+(m-1),1});
        %str2num is necessary because the number values are in string form
    end
end

%Isolation of Triangles Section

data2 = (datacells12(2,1):9:vertcells(2,1))';
%provides the start point for each row of data in the Triangles Section

n = length(data2);
%provides the correct number of iterations for the outer loop

for k = 1:n
    
    for m = 1:6
        %there are six iterations becauset there are six columns of data
        
        tri(k,m) = str2num(mapmesh{data2(k,1)+(m-1),1});
        
    end
end

%Isolation of Vertices Color Section

% data3 = (vert1:12:vertcells(3,1))';
%
% p = length(data3);
% cd ..
% for k = 1:p
%
%     for m = 1:3
%
%         tricolor(k,m) = str2num(mapmesh{data3(k,1)+(m-1),1});
%
%     end
% end

%Isolation of Vertices Attributes Section

% data4 = (vert2:4:s)';
%
% q = length(data4);
%
% for k = 1:p
%
%     for m = 1:2
%
%         vertattrib(k,m) = str2num(mapmesh{data4(k,1)+(m-1),1});
%
%     end
% end

% Unipolar = tricolor(:,1);
% Bipolar = tricolor(:,2);
% LAT = tricolor(:,3);
% cd(outputdatadir);

%     dlmwrite('Mesh.pts', norm(:,1:3), 'delimiter', '\t');
%     dlmwrite('Mesh.itri', tri(:,1:3), 'delimiter', '\t');
%     Field1.pts = norm(:,1:3);
%     Field1.tri = tri(:,1:3);
%     save(outfilename_fld,'Field1');
VER=norm(:,1:3);
ITRI=tri(:,[2 1 3])+1;

[VER,ITRI]=cleanTRIs(VER,ITRI);

% savetri([outfname '_mesh.tri'],VER,ITRI);
% save(outfilename_Upot,'Unipolar');
% save(outfilename_Bpot,'Bipolar');
% save(outfilename_LAT,'LAT');
% end
cd(origpath)