% tri2gra.m 
% computes a matrix based on the graph associated with a triangulated surface
% G=tri2gra(VER,ITRI,mode); 
% The elements of G specify: 
% mode = 1: the adjacency matrix
% mode = 2: the adjacency matrix weighted by the length of the connections
% mode = 3: the (3-D) distances between all nodes
% mode = 4: the (3-D) distances along a line connecting
%             the nodes if the line lies entirely within the surface; else it is zero
%             run distgra.m to make compute the shortest distances between all nodes
%             with the connections running entirely within the surface.
% run gra2dist, using the output of this function, to obtain the shaortest
% path distances between all nodes.
% 20041028; A. van Oosterom; mode 4 streamlined; initial step now in fact
% effective

function  G=tri2gra(VER,ITRI,mode)
nver=size(VER,1);
G=zeros(nver);
ntri=size(ITRI,1);

if mode==1,
	% adjacency matrix
	for i=1:ntri;
		for j=1:3;
			ja=icyc(j+1,3);
			G(ITRI(i,j),ITRI(i,ja))=1;
			G(ITRI(i,ja),ITRI(i,j))=1;
		end
	end
end

if mode==2,
	% adjacency weighted by actual distance between neighbours
	for i=1:ntri;
		for j=1:3;
			ja=icyc(j+1,3);
			r1=VER(ITRI(i,j),:);
			r2=VER(ITRI(i,ja),:);
			G(ITRI(i,j),ITRI(i,ja))=norm(r1-r2);
			G(ITRI(i,ja),ITRI(i,j))=G(ITRI(i,j),ITRI(i,ja));
		end
	end
end

if mode==3,
	% distances between all nodes in 3-D
	for j1=1:nver,
		r1=VER(j1,:);
		for j2=j1:nver;
			r2=VER(j2,:);
			dist=sqrt((r1-r2)*(r1-r2)');
			G(j1,j2)=dist;
			G(j2,j1)=dist;
		end
	end
end

if mode==4,
	% distances between all nodes if they can be connected by a straight
	% line passing entirely through the interior
	% start with connecting the direct neigbours as for mode=2
	for i=1:ntri;
		for j=1:3;
			ja=icyc(j+1,3);
			r1=VER(ITRI(i,j),:);
			r2=VER(ITRI(i,ja),:);
			G(ITRI(i,j),ITRI(i,ja))=sqrt((r1-r2)*(r1-r2)');
            G(ITRI(i,ja),ITRI(i,j))=G(ITRI(i,j),ITRI(i,ja));
		end
	end

	% now fill in any of the other values 
	for j1=1:nver,
		for j2=j1+1:nver,
            if G(j1,j2)==0,
		        view=nodeview(VER,ITRI,j1,j2);
		        if view==1,
			        G(j1,j2)=sqrt((VER(j1,:)-VER(j2,:))*(VER(j1,:)-VER(j2,:))');
			        G(j2,j1)=G(j1,j2);
		        end
            end
        end
	end
end




