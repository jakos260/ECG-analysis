function [pnt, dhk] = loadtri_insbruck(fn);
% avo 2005-06-17
fn=fn(fn~=' ');

fid = fopen(fn, 'rt');
if fid~=-1

  % read the vertex points
  Npnt = fscanf(fid, '%d', 1);
  pnt  = fscanf(fid, '%f', [3, Npnt]);
  pnt  = pnt(1:3,:)';

  % if present, read the triangles
  if (~(feof(fid)))
    [Ndhk, count] = fscanf(fid, '%d', 1);
    if (count ~= 0)
      dhk = fscanf(fid, '%d', [3, Ndhk]);
      dhk = dhk(1:3,:)';
    end
  else
    dhk = [];
  end
  fclose(fid);

else
  error('unable to open file');
end


