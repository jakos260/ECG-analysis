function [VER,ITRI]=loadAsciStl(filename)
% function to import ASCII STL files into matlab.  This update uses the
% textscan function instead of reading the file line by line, which can
% significantly reduce import times for very large STL files.
%--------------------------------------------------------------
%
% inputs:       filename:   string
%
%-----------------------------------------------------------------

%open file
fid=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.
if fid == -1
    error('File could not be opened, check name or path.')
end

%=====================================================

fid=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.

fmt = '%*s %*s %f32 %f32 %f32 %*c %*s %*s %*c %*s %f32 %f32 %f32 %*c  %*s %f32 %f32 %f32 %*c %*s %f32 %f32 %f32 %*c  %*s %*c  %*s %*c ';
C=textscan(fid, fmt, 'HeaderLines', 1);
fclose(fid);
%extract normal vectors and vertices
tnorm = cell2mat(C(1:3));
tnorm = tnorm(1:end-1,:); %strip off junk from last line

v1 = cell2mat(C(4:6));
v2 = cell2mat(C(7:9));
v3 = cell2mat(C(10:12));

if isnan(C{4}(end))
    v1 = v1(1:end-1,:); %strip off junk from last line
    v2 = v2(1:end-1,:); %strip off junk from last line
    v3 = v3(1:end-1,:); %strip off junk from last line
end

v_temp = [v1 v2 v3]';
v = zeros(3,numel(v_temp)/3);

v(:) = v_temp(:);
v = v';

[VER,ITRI]=fv2pt(v,length(v)/3);%gets points and triangles

%%
function [p,t]=fv2pt(v,fnum)

%gets points and triangle indexes given vertex and facet number
c=size(v,1);

%triangles with vertex id data
t=zeros(3,fnum);
t(:)=1:c;

%now we have to keep unique points fro vertex
[p,i,j]=unique(v,'rows'); %now v=p(j) p(i)=v;
t(:)=j(t(:));
t=t';

