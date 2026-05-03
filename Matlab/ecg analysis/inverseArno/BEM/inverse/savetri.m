function savetri(fn,VER,ITRI);
% savetri stores vertices and triangles from a MBFYS triangulation in file
% savetri(fn,VER,ITRI);
% nodes now stored in 8.6f format


f = fopen(fn, 'w');
[nver dim]=size(VER);
fprintf(f,'%d\n ',nver);
for i=1:nver;
      fprintf(f,'%5d %8.6f %8.6f %8.6f\n',i,VER(i,1:3));
end

[ntri dim]=size(ITRI);
fprintf(f,'%d\n ',ntri);
for i=1:ntri;
      fprintf(f,'%5d %5d %5d %5d\n',i,ITRI(i,1:3));
end
sprintf('\ntriangle specs written to file: %s\n',fn); 
fclose(f);
