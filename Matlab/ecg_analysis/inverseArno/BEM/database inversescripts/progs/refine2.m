% refine2.m
% function [VER,ITRI,itri1,itri2]=refine2(VER,ITRI,edge)
% refine both triangles that have edge as an edge edge=[node1 node2]
% by inserting a node half-way the edge;
% itri1 and itri ~0 are the original indices of the refined triangles

% 20111018

function [VER,ITRI,itri1,itri2]=refine2(VER,ITRI,edge)

nver=size(VER,1);
ntris=size(ITRI,1);
itri1=0; itri2=0; % original indices of the refined triangles
ADD1=[]; ADD2=[];

% identify neighbouring triangle nbtri (sharing edge) 
btri=find(sum(ismember(ITRI(:,1:3),edge),2)==2);
btri=sort(btri)
if isempty(btri), ' no trinagle matching the edge found', return, end

nver=nver+1;
VER(nver,:)=mean(VER(edge,:));

TRIS=[ITRI (1:ntris)'];
TRIS=[TRIS(btri,1:4);  TRIS(btri,[2 3 1 4]);  TRIS(btri,[3 1 2 4])]                    



 
    k=find(TRIS(:,1)==edge(1)&TRIS(:,2)==edge(2)&TRIS(:,4));
    
    % specify the two parts of the refined triangle btri(i)
    
    if ~isempty(k)
       k=k(1);
       % specify the two parts of the refined triangle 
       ADD1=[edge(1) nver  TRIS(k,3);
            nver edge(2) TRIS(k,3)]
       itri1=TRIS(k,4);
    end
    
    k=find(TRIS(:,1)==edge(2)&TRIS(:,2)==edge(1));
    if ~isempty(k)
        k=k(1);
        % specify the two parts of the refined triangle 
        ADD2=[nver   edge(1) TRIS(k,3) ;
             nver    TRIS(k,3) edge(2)]; 
        itri2=TRIS(k,4);
    end
    
    if itri1~=0, ITRI=[ITRI(1:itri1-1,:);ADD1;ITRI(itri1+1:end,:)]; end
    if itri1~=0, itri2=itri2+1; end
    if itri2>1, ITRI=[ITRI(1:itri2-1,:);ADD2;ITRI(itri2+1:end,:)]; end
    
    
    
    
    
          
end         
   
   
        



       

   
   
   


   

