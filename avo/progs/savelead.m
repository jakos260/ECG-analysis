function savelead(name, matrix)
% SAVEASCI	Save data in an ascii file.
%
%		Usage: saveasci('file',M);
%                  or: saveasci('file',M,comm);
%
%               comm can be used for comment lines
%		See also LOADASCI.
%
%		Thom Oostendorp, MF&BF University of Nijmegen, the Netherlands
%		(AvO's version)

f=fopen(name, 'w');

N=size(matrix);
fprintf(f, '%d %d\n', N(1), N(2)+1);

for i=1:N(1)
  k=0;
  fprintf(f, '%d ', i);
  for j=1:N(2)
    fprintf(f, '%f ', matrix(i,j));
    k=k+1;
  end
  if k~=0
    fprintf(f, '\n');
  end
end
fclose(f);

