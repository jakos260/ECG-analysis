function saveasci(name, matrix, tekst);
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
fprintf(f, '%d %d\n', N(1), N(2));

for i=1:N(1)
  k=0;
  for j=1:N(2)
    fprintf(f, '%f ', matrix(i,j));
    k=k+1;
%     if k==5,
%       k=0;
%       fprintf(f,'\n');
%     end
  end
  if k~=0
    fprintf(f, '\n');
  end
end

fprintf('\nmatrix written to file: %s\n\n',name);
if exist('tekst','var')
   M=size(tekst);
   fprintf(f, '%d %d\n', M(2), M(1));
   for i=1:M(1)
      fprintf(f, '%s\n',tekst(i,:));
   end
end
fclose(f);

