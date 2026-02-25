function [el, lab] = readelec(fn);

% READ_ELEC reads "la/mu" electrode parameters from a MBFYS electrode file
%	which are used to position them on a triangulated surface
%
%	[el, lab] = readelec(filename), 
% 	  where el = [dhk, la, mu]
%	  and lab contains the electrode labels (if present)
%
% See also READ_TRI, TRANSFER_ELEC

% (c) Robert Oostenveld, 1998

fid = fopen(fn, 'rt');
if fid~=-1

  % read the number of electrodes
  Nel = sscanf(fgetl(fid), '%d'); 
 
  % read the electrode triangle, lambda and mu
  for i=1:Nel
    str = fgetl(fid);
    el(i,:)  = sscanf(str, '%f %f %f')';
    indx = find(str=='!');
    if (indx)
      lab(i,:) = sprintf('%6s', str((indx+1):length(str)));
    end
  end
  fclose(fid);

else
  error('unable to open file');
end

