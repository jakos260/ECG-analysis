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
function [VER,ITRI] = loadSTL( filename )


f=fopen(filename);
if (f==-1)
  error('\nCannot open %s\n\n', filename);
  VER=[];
  ITRI=[];
  return;
end
str = fread(f,80,'uchar');
if  length(str) == 80    
    Nfaces=fread(f,1,'uint32');
    ITRI=zeros(Nfaces,3);
    VER = [];
    for i=1:Nfaces
   
        a  = fread(f,3,'float32');% normals
        v1 = fread(f,3,'float32');
        v2 = fread(f,3,'float32');
        v3 = fread(f,3,'float32');
        n = fread(f,1,'short');  
        if n~=0
            stop = 1;
        end
        if isempty(VER)
            i1=[];
            i2=[];
            i3=[];
        else            
            i1 = find(VER(:,1) == v1(1) & VER(:,2) == v1(2) & VER(:,3) == v1(3));
            i2 = find(VER(:,1) == v2(1) & VER(:,2) == v2(2) & VER(:,3) == v2(3));
            i3 = find(VER(:,1) == v3(1) & VER(:,2) == v3(2) & VER(:,3) == v3(3));
        end
        if isempty(i1)
            VER=[VER; v1'];
            i1 = size(VER,1);
        end
        if isempty(i2)
            VER=[VER; v2'];
            i2 = size(VER,1);
        end
        if isempty(i3)
            VER=[VER; v3'];
            i3 = size(VER,1);
        end
        ITRI(i,:) = [i1 i2 i3];  
    end
end


function readSTLtext( filename )
% STL2tri a .STL file to noid.tri file with
% vertices and triangles, a MBFYS triangulation file
%

  [raw_nodes,normals] = rw_stl(filename);

  done = zeros(1,length(raw_nodes));
  pnt = zeros(size(raw_nodes));
  nr_vertex = 0;
	
  for j = 1:length(raw_nodes)
		a1=find(pnt(:,1) == raw_nodes(j,1));
		a2=find(pnt(:,2) == raw_nodes(j,2));
		a3=find(pnt(:,3) == raw_nodes(j,3));
	  if isempty(a1) | isempty(a2) | isempty(a3)
		  nr_vertex = nr_vertex+1;
		  pnt(nr_vertex,1:3)= raw_nodes(j,:);
		end;
	end;
	pnt = pnt(1:nr_vertex,:);
  
  ndhk = length(raw_nodes)/3
  dhk = zeros(ndhk,3);
	nt = 0;
  for j = 1:3:length(raw_nodes)
    nt = nt +1;	
		for i=0:2
			a1=find(pnt(:,1) == raw_nodes(j+i,1));
			a2=find(pnt(:,2) == raw_nodes(j+i,2));
			a3=find(pnt(:,3) == raw_nodes(j+i,3));
			if ~isempty(a1) & ~isempty(a2) & ~isempty(a3)
				dhk(nt,i+1)=a3(find(a1(find(a1==a2))==a3));
			else
				oo=1
			end
		end
	end	
% 	nt
% 	a=dhk(3,:);
% 	dhk(3,:)=dhk(2,:);
% 	dhk(2,:)=a;
% 	pnt=0.94*pnt;
  
%--------------------------------------------------------------------------------%

function [raw_nodes,normals] = rw_stl(filename)

fid = fopen(filename', 'rt');

%raw_nodes = zeros(61128,3);
%normals = zeros(20376,3);
raw_nodes = zeros(5100,3);
normals = zeros(1700,3);

if fid~=-1
  % read header
  fgetl(fid);
  i = 0;
  j = 0;
  while~(feof(fid))
    [fid,x] = rw_normal(fid);
    if (~feof(fid))
			i = i + 1;
			normals(i,1:3) = x;
			j = j + 1; [fid,raw_nodes(j,1:3)] = rw_vertex(fid); 
			j = j + 1; [fid,raw_nodes(j,1:3)] = rw_vertex(fid); 
			j = j + 1; [fid,raw_nodes(j,1:3)] = rw_vertex(fid); 
		end
  end;
  fclose(fid);
%   savenormals('newatrium.norm',normals);
end

