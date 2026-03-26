function [result, extraresult]=loadasci(name);
% LOADASCI	Load an asci matrix file.
% version 20021128

%
%		Usage: m         = loadasci('file');
%
%		LOADASCI('file') returns the matrix stored in 'file'.
%
%		See also SAVEASCI.
%
%		Thom Oostendorp, MF&BF University of Nijmegen, the Netherlands
%		(Adriaan's version)

f=fopen(name, 'r');
if (f==-1)
  fprintf('\nCannot open %s\n\n', name);
  result=0;
  return;
end

[N,nr]=fscanf(f,'%d',2);
if (nr~=2),
  fclose(f);
  fprintf('Error reading %s\n\n', name);
  result=0;
  return;
end;

M=fscanf(f,'%f',[N(2) N(1)]);
fclose(f);
S=sprintf('\n%s contains %d rows and %d columns\n', name, N(1), N(2));
%disp(S);
result=M';
