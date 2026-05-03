function [result, extraresult]=loadmat(name);
% LOADMAT	Load a MFBF matrix file.
%
%		Usage: m         = loadmat('file');
%		   or  [m,extra] = loadmat('file');
%
%		LOADMAT('file') returns the matrix stored in 'file' and
%		the extra information stored at the bottom of that file.
%		LOADMAT works for binary as well as asci matrix files.
%
%		See also SAVEMAT.
%
%		Thom Oostendorp, MF&BF University of Nijmegen, the Netherlands

f=fopen(name);
if (f==-1)
  fprintf('\nCannot open %s\n\n', name);
  result=0;
  extraresult='';
  return;
end

[N,nr]=fscanf(f,'%d',2);
if (nr~=2)
  fclose(f);
  f=fopen(name);
  [magic ,nr]=fread(f,8,'char');
  if (char(magic')==';;mbfmat')
    fread(f,1,'char');
    hs=fread(f,1,'long');
    fread(f,1,'char');
    fread(f,1,'char');
    fread(f,1,'char');
    N=fread(f,2,'long');
    M=fread(f,[N(2),N(1)],'double');
  else
    fclose(f);
    f=fopen(name);
    N=fread(f,2,'long');
    posit=ftell(f);fseek(f,0,'eof');
    eofpos=ftell(f);fseek(f,posit,'bof');
    if (eofpos-posit)/8 >= N(2)*N(1)
        M=fread(f,[N(2),N(1)],'double');
    else
        M=fread(f,[N(2),N(1)],'float');
    end
  end
else
  M=fscanf(f,'%f',[N(2) N(1)]);
end
[extra,nextra]=fread(f,1000,'char');
fclose(f);
% S=sprintf('\n%s contains %d rows and %d columns\n', name, N(1), N(2));
% disp(S);
if (nextra~=0)
  S=sprintf('%s contains the following extra information:\n', name);
  disp(S);
  disp(setstr(extra'));
end
result=M';
extraresult=setstr(extra');
