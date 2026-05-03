% refine3.m
% adapt [VER,ITRI] by inserting one extra node at the
% center of gravity of itri and remeshing the triangle
% function [VER,ITRI]=refine3(VER,ITRI,itri)
% VER, ITRI : vertex and triangle indices
% A. van Oosterom; 20090403

function [VER,ITRI]=refine3(VER,ITRI,itri)
         trits=ITRI(itri,:);
         VER=[VER; mean(VER(trits,:))];
         nver=size(VER,1);
         ITRI(itri,:)=[trits(1) trits(2) nver];
         ITRI=[ITRI;
              trits(2) trits(3) nver;
              trits(3) trits(1) nver;];
end
