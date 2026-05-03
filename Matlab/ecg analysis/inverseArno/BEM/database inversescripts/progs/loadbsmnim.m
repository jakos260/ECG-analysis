function [result,M]=loadbsmnim(name);
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

sel=[64  63  53  54  55  56  57  58  59   1 
      2   3   4   5   6   7   8   9  10  11 ...
     12  13  14  15  16  17  18  19  20  21 ...
     22  23  24  25  26  27  28  29  30  31 ...
     32  33  34  35  36  37  38  39  40  41 ...
     42  43  44  45  46  47  48  60  61  62 ...
     49  50  51  52];
result=-M(sel,:)*0.0001831054;
