% dsa_plus.m
% [dsa_s,jsing,obs_used]=dsa_plus(VER,ITRI,obs,sin_mode,sap);
% dsa_s: the dsa values of all three nodes of triangles (ITRI) as viewed from obs
%        size:(3,ntri)
% treatment of the auto_solid angle by replacing singular obs by a nearby,
% interior position;
% for a closed surface, the interior values are defined positive with sum(sum(dsa_s))=4 pi; 
% for any other location it is zero
% jsing~0 denotes the index of the singular vertex; 
% in this case obs(output) is the nearby observation point actually used used 
% function calls dsa, a function that does not treat singular observation points 
% if sin_mode <0, the replaced field point is positionned just outside the surface;
% sap is the solid angle threshold below which dsa uses omega(i)=OMEGA/3 
% A. van Oosterom; 2015_07_12

% TEST VERSION    

    %function [DSA,jsing,obs_used]=dsa_plus(VER,ITRI,obs,sin_mode,sap)
    clear
    
    obs=[0 0 0.0001]
    zz=-0.5
    VER=[0 0 0; 2 0 zz; -1  1 zz; -1 -1 zz ; 2 -1 zz]
    ITRI=[4 3  1; 2 5  4; 2 4  1; 2 1 3];
    ITRI=ITRI(:,[2 1 3])
    sin_mode=0; 
    sap=1.e-3;
    obs_used=obs;
  
    [DSA,sing]=dsa(VER,ITRI,obs,sap); % distrib  dsa values(3 vertices of all triangles ITRI)
    DSA
    sum(sum(DSA(:,[1 3 4])),2)
    pause
    
    
    k=find(sing~=0) % triangles whose nodes comprize obs
    
    if ~isempty(k)  
        % handling of singular observation such that the contributions from the individual nodes of
        % TRIS carrying obs is accounted for:
        TRIS=ITRI(k,:)
        % identify the singular node of VER
        ss=std(sort(TRIS,2));
        jsing=find(ss==0)
        [edge, ntris,nln]=findloop(TRIS,jsing);
        %edge=reverse(edge);
        VERN=VER(edge,:);
        lamb=1.e-5;
        
        if size(edge,2)==3,
           pivot=mean(VERN);
           dist=pivot-obs;
           dist=lamb*dist/norm3d(dist)
           obs;
           obs_sing=obs+dist;
           DSA(:,k)
           DSA_k=dsa(VER,TRIS,obs_sing,sap)
           pause   
           
        else
        % split_up edge into two halves
        lista=1:round(ne/2);
        listb=reverse(round(ne/2)+1: ne);
        
        % perform a tesselation of the interior of the edge
        TILES=make_strip(VERN,lista,listb);
        sas=solida(VERN,TILES,VER(jsing,:));
        [ma,ima]=max(sas);
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
        
        if sin_mode<0,  obs_used=2*VER(jsing,:)-obs;end
        
        DSA=dsa(VER,ITRI,obs,sap)  % size(3,ntri)
        
        end
    end
    
     jsing=k
     obs_used
     sum(sum(DSA(:,[1 2:3])))
     
     
    
    
    
    
    
    
    
    
    
    
    
   
    
    
    
    
    
    
    
  