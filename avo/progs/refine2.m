% refine2.m
% function [VER,ITRI,btri1,btri2]=refine2(VER,ITRI,edge)
% refine both triangles that share edge
% by inserting a node half-way the edge;
% btri1 and btri2  are the original indices of the refined triangles

% 20140421

function [VER,ITRI,tri1,tri2]=refine2(VER,ITRI,edge)

nver=size(VER,1);
ntris=size(ITRI,1);

%itri1; itri2; % original indices of the refined triangles

% identify neighbouring triangle nbtri (sharing edge)


%edge


btri=find(sum(ismember(ITRI(:,1:3),edge),2)==2);

tri1=btri(1);
tri2=btri(2);

if isempty(btri), ' no triangle matches the edge', return, end

BTRIS=ITRI(btri,:);

nver=nver+1;
VER(nver,:)=mean(VER(edge,:));

odd1=BTRIS(1,:);
odd1(ismember(odd1,edge))=[];

odd2=BTRIS(2,:);
odd2(ismember(odd2,edge))=[];

NEWTRIS=[odd1 edge(1) nver;      % replace 1  
         odd1  nver   edge(2);   % ADD
         nver edge(1) odd2;      % replace 2
         odd2 edge(2) nver;];    % ADD

ITRI(btri(1),:)= NEWTRIS(1,:);
ITRI(btri(2),:)= NEWTRIS(3,:);

ITRI=[ITRI;NEWTRIS([2 4],:)];



     
   
