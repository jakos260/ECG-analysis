% separatetwotris.m
% function [VER1,ITRI1,VER2,ITRI2]=separatetwotris([VER, ITRI, nver2,itri2);
% "inverse" of addtwotris
% 20071016;

function [VER1,ITRI1,VER2,ITRI2]=addtwotris(VER,ITRI,nver2,ntri2);

nver=size(VER, 1);
ntri=size(ITRI,1);

nver1=nver-nver2;
VER1=VER(1:nver1,:);

VER2=VER(nver1+1:nver,:);
ntri1=ntri-ntri2

ITRI1=ITRI(1:ntri1,:);
ITRI2=ITRI(ntri1+1:ntri,:)-nver1;

