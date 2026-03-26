function savevertices(name, VER, tekst);
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

N=size(VER,1);
fprintf(f, '%d %d\n', N(1), 4);

for i=1:N
    fprintf(f, '%i %f %f %f\n', i,VER(i,1),VER(i,2),VER(i,3));

end

fprintf('\nmatrix written to file: %s\n\n',name);
if exist('tekst','var'),
   M=size(tekst);
   fprintf(f, '%d %d\n', M(2), M(1));
   for i=1:M(1),
      fprintf(f, '%s\n',tekst(i,:));
   end
end
fclose(f);

