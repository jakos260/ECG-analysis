function [result]=loadbsm(name);
% LOADBSM	Load ruwe HOEKEMA ALBUM 
%
%		Usage: m         = loadmat('file');
%		   or  [m,extra] = loadmat('file');
%
%		LOADMAT('file') returns the matrix stored in 'file' and
%		the extra information stored at the bottom of that file.
%		LOADMAT works for binary as well as asci matrix files.
%

f=fopen(name);
if (f==-1)
  fprintf('\nCannot open %s\n\n', name);
  result=0;
  return;
end

N=fread(f,2,'long');
M=fread(f,[N(1),N(2)],'short');

fclose(f);
S=sprintf('\n%s contains %d rows and %d columns\n', name, N(1), N(2));
disp(S);
result=M*0.0001831054;
