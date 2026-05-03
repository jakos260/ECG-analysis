function [LIST]=loadchar(name);
% LOADCHAR Load characters from a file file.
%
%		Usage: m         = loadchar('file');
%
%		LOADCHAR('file') returns in a matrix
%               the character stored in 'file'.
%               the number of rows,  and length string should be heading the data

f=fopen(name, 'r');
if (f==-1)
  fprintf('\nCannot open %s\n\n', name);
  LIST=0;
  return;
end

[N,nr]=fscanf(f,'%d',2);
if (nr~=2)
   fclose(f);
   fprintf('Error reading %s\n\n', name);
   LIST=0;
  return;
end;
bl=blanks(N(1))';
LIST=bl(:,ones(1,N(2)));
line=fgetl(f);
for i=1:N(1),
    line=fgetl(f);
    LIST(i,:)=line(1:N(2));
end
fclose(f);
S=sprintf('\n%s contains %d rows and %d columns\n', name, N(1), N(2));
disp(S);
