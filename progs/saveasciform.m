function saveasciform(name, M, format, COMM);
% Save data in a formatted ascii file.
%		Usage: saveasciform('file',M,format);
%          or: saveasciform('file',M,format,COMM);
%               
% example:
% form ='%3d %6.4f %4d %4d %6.4f %4d %4d %6.4f %6.4f %6.4f\n'
% saveasciform('listing.lst',RES2,form)
% COMM: matrix containing comment lines: example: COMM=['comment';
%                                                       'lines   ']
% 2013_11_03 A. van Oosterom


f=fopen(name, 'w');

N=size(M)
fprintf(f, '%d %d \n', N(1), N(2));

for i=1:N(1),
   fprintf(f,format ,M(i,1:N(2))); 
end

fprintf(f, '\n');

fprintf('\nmatrix written to file: %s\n\n',name);

if nargin > 3,
   m=size(COMM);
   fprintf(f, '%d %d\n', m(2), m(1));
   for i=1:m(1),
       fprintf(f, '%s\n',COMM(i,:));
   end
end


fclose(f);













