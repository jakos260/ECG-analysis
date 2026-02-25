
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

cd(path2use)

outfilename_Upot = 'Unipolar';
outfilename_Bpot = 'Bipolar';
outfilename_LAT = 'LAT';
outfilename_fld = 'Field1';
outputdatadir = 'Parsed_MeshFile';
h2 = dir(['../',outputdatadir]);
if isempty(h2) == 1
    mkdir(['../',outputdatadir]);
    h2 = dir(['../',outputdatadir]);
end

% if Len_h1 > 1
%     fname = input('More than one .mesh file.  Please enter the desired file. Example("blauer3mo.mesh")', 's');
% else
%     fname = h1.name;
% end

% for iname =1:length(h1)
fname = fullfile(path2use,filename);%  h1(iname).name;
outfname = fname(1:end-5);
cd(path2use)

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
cd ..
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
cd(outputdatadir);

%     dlmwrite('Mesh.pts', norm(:,1:3), 'delimiter', '\t');
%     dlmwrite('Mesh.itri', tri(:,1:3), 'delimiter', '\t');
%     Field1.pts = norm(:,1:3);
%     Field1.tri = tri(:,1:3);
%     save(outfilename_fld,'Field1');
VER=norm(:,1:3);
ITRI=tri(:,[2 1 3])+1;

[VER,ITRI]=cleanTRIs(VER,ITRI);

savetri([outfname '_mesh.tri'],VER,ITRI);
% save(outfilename_Upot,'Unipolar');
% save(outfilename_Bpot,'Bipolar');
% save(outfilename_LAT,'LAT');
% end
cd(origpath)