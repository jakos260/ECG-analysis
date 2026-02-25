function savemat(name, M, extra)
% function savemat(name, M, extra);
% SAVEMAT	Save a MFBF matrix file.
%
%		Usage: savemat('file', M);
%		   or  savemat('file', M, extra);
%
% SAVEMAT('file',M) saves matrix M into 'file' in MFBF
%		binary format. 
% SAVEMAT('file',M, extra) appends the
%		string <extra> to the bottom of that file.
%
%		See also LOADMAT.
%
%		Thom Oostendorp, MF&BF Radboud University, Nijmegen, the Netherlands
% size(M)
M=M';
f=fopen(name, 'wb');
n=size(M);
fwrite(f, n(2), 'long');
fwrite(f, n(1), 'long');
fwrite(f, M, 'float');

if exist('extra')
  [m,~]=size(extra);
  for i=1:m
    fwrite(f, extra(i,:));
    fwrite(f, char(10));
  end
end
fclose(f);
