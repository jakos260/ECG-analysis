function savechar(name,LIJST)
%function savechar(name,LIJST)

f=fopen(name, 'w');
m=size(LIJST)
fprintf(f, '%d %d\n', m(1), m(2));
for i=1:m(1),
  fprintf(f, '%s\n',LIJST(i,:));
end
fclose(f);
