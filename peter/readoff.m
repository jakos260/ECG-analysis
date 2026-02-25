function [pnt,dhk]=readoff(fn)

% OFF
% 4  2  0 
% 0.0 0.0 0.0
% 1.0 0.0 0.0
% 1.0 1.0 0.0
% 0.0 1.0 0.0
% 3   0 1 2
% 3   0 2 3

% First the file has the letters "OFF" in the header. 
% Next, the "4" indicates that there are four floating 
% point vertices to come. The "2" indicates that there 
% are two polygons. The "0" is a space for the number 
% of edges, but we won't use that. The "3" before each 
% polygon indicates that there are three vertices 
% (forming a triangle). The integers are pointers into 
% the vertex list, with the first point having index "0".


fn=fn(fn~=' ');

fid = fopen(fn, 'rt');
if fid~=-1
  tmp=fgetl(fid);
  % read the vertex points
  Npnt = fscanf(fid, '%d', 3);
  Ndhk=Npnt(2);
  Nedge=Npnt(3);
  Npnt=Npnt(1);
  pnt  = fscanf(fid, '%f', [3, Npnt])';
%   pnt  = pnt(2:4,:)';

  % if present, read the triangles
  if (~(feof(fid)))
    dhk = fscanf(fid, '%d',[4 Ndhk]);
    dhk = dhk(2:4,:)'+1;
  else
    dhk = [];
  end
  fclose(fid);

else
  error('unable to open file');
end

