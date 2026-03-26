% find_edges.m
% ALLEDGES=find_edges(ITRI)
% ALLEDGES(:,1:2): edge indices;  ALLEDGES(:,3): triangle index;

% date: 2014023

function ALLEDGES=find_edges(ITRI)

   ntri=size(ITRI,1);
   ALLEDGES=[ITRI(:,1:2) (1:ntri)'; ITRI(:,2:3) (1:ntri)'; ITRI(:,[3 1]) (1:ntri)'];
   ALLEDGES=sortrows(ALLEDGES,1);
 

   
        