function saveasciform(name, M, format, comm);
% Save data in a formatted ascii file.
%		Usage: saveasciform('file',M,format);
%          or: saveasciform('file',M,format,comm);
%               comm being comment lines
% example:
% form ='%3d %6.4f %4d %4d %6.4f %4d %4d %6.4f %6.4f %6.4f\n'
% saveasciform('listing.lst',RES2,form)
% 20061017 A. van Oosterom


f=fopen(name, 'w');

N=size(M)
fprintf(f, '%d %d \n', N(1), N(2));

for i=1:N(1),
   fprintf(f,format ,M(i,1:N(2))); 
end

fprintf(f, '\n');

fprintf('\nmatrix written to file: %s\n\n',name);

if nargin > 3,
   m=size(comm);
   fprintf(f, '%d %d\n', m(2), m(1));
   for i=1:m(1),
       fprintf(f, '%s\n',comm(i,:));
   end
end
fclose(f);
