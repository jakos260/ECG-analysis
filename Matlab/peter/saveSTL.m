% /** @brief stl file format handling
%  *  facet normal 0.0 0.0 1.0
%  *      outer loop
%  *          vertex  1.0  1.0  0.0
%  *          vertex -1.0  1.0  0.0
%  *          vertex  0.0 -1.0  0.0
%  *      endloop
%  *  endfacet
%  * Binary STL files consist of a 80 byte header line that can be
%  * interpreted as a comment string. The following 4 bytes interpreted
%  * as a long integer give the total number of facets. What follows is
%  * a normal and 3 vertices for each facet, each coordinate represented
%  * as a 4 byte floating point number (12 bytes in all). There is a 2
%  * byte spacer between each facet. The result is that each facet is represented
%  * by 50 bytes, 12 for the normal, 36 for the 3 vertices, and 2 for the spacer.
%  * @param str the name of the file
%  */
function saveSTL( filename,VER,ITRI )

f=fopen(filename, 'wb');

if (f==-1)
    fprintf('\nCannot open %s\n\n', filename);
    return;
end
header = [];
for i=1:80
    header= [header ' '];
end
count = fwrite(f, header, 'uchar');
fwrite(f,size(ITRI,1),'uint32');
normals = cross(VER(ITRI(:,2),:) -VER(ITRI(:,1),:), VER(ITRI(:,3),:) -VER(ITRI(:,1),:));
for i=1:size(ITRI,1)    
    fwrite(f,normals(i,:),'float32');% normals
    fwrite(f,VER(ITRI(i,1),:),'float32');
    fwrite(f,VER(ITRI(i,2),:),'float32');
    fwrite(f,VER(ITRI(i,3),:),'float32');
    fwrite(f,50,'short');
end
fclose(f);
