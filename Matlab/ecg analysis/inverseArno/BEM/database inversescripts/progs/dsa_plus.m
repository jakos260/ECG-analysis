% dsa_plus.m
% [dsa_s,jsing,obs]=dsa_plus(VER,ITRI,obs,sin_mode);
% dsa_s: the dsa values of all triangles (ITRI) as viewed from obs
% size:(3,ntri)
% treatment of the auto_solid angle by replacing singular obs by a nearby,
% interior position; for a closed surface,
% the interior values are defined positive with sum(sum(dsa_s))=4 pi; 
% for any othe location it is zero
% jsing~0 denotes the index of the singular vertex; 
% in this case obs(output) is the nearby observation point actually used used 
% function calls dsa, a function that does not treat singular observation points 
% if sin_mode <0, the replaced field point is positionned just outside the
% surface
% A. van Oosterom; 2012_02_14

    
    function [DSA,jsing,obs]=dsa_plus(VER,ITRI,obs,sin_mode,sap)
   
    jsing=0;
    
    [DSA,sing]=dsa(VER,ITRI,obs,sap); % distrib  dsa values(3 vertices of all triangles)
    k=find(sing==4);
    
    if ~isempty(k)
        % handling of singular observations such that contributions from the individual nodes of
        % any triangle carrying the singularities is accounted for
        
        TRIS=ITRI(k(1),:);
        trisnodes=unique(TRIS(:));
        ntrn=size(trisnodes,1);
        jsing=find(norm3d(VER(trisnodes,:)-ones(ntrn,1)*obs)<1.e-7);
        jsing=trisnodes(jsing);
        
        
        [bver,BTRI]=buren(ITRI,jsing);
        edge=findedge(VER,BTRI(:,2:4),bver(1)); % edge is a closed string
        edge(end)=[];  % open the string
        ne=size(edge,2);

        VERN=VER(edge,:);
        
%         normalize lengths of edges meeting at jsing (possibly not required)
%         normi=norm3d(VERN-ones(ne,1)*VER(jsing,:));
%         VERN=ones(ne,1)*VER(jsing,:)+(VERN-ones(ne,1)*VER(jsing,:))./(normi*ones(1,3));
      
        % split_up edge into two halves
        lista=1:round(ne/2);
        listb=reverse(round(ne/2)+1: ne);
        % perform a tesselation of the interior of the edge
        TILES=make_strip(VERN,lista,listb);
        sas=solida(VERN,TILES,VER(jsing,:));
        [ma,ima]=max(sas);
        
        lamb=1.e-5;
        if ma>0,
            pivot=mean(VERN(TILES(ima,:),:));
            obs=(1-lamb)*VER(jsing,:)+lamb*pivot;
        else
            if min(sas)==0,
               tris=TILES(ima,:);
               center=mean(VERN(tris,:));
               normal=cross(VERN(tris(1),:)-center,VERN(tris(2),:)-center);
               obs=VER(jsing,:)-1.e-5*normal/norm(normal);
            end
        end
        if ma<0,
            [mi,imi]=min(sas);
            pivot=mean(VERN(TILES(imi,:),:));
            obs=(1+lamb)*VER(jsing,:)-lamb*pivot;
        end
        
        if sin_mode<0, obs=2*VER(jsing,:)-obs;end
        
        
        DSA=dsa(VER,ITRI,obs,sap);  % size(3,ntri)
    end
    
    
    
    
    
    
    
    
    
    
    
    
   
    
    
    
    
    
    
    
  