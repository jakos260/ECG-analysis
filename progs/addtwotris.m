% addtwotris.m
% function [VER,ITRI]=addtwotris([VER1, ITRI1,VER2,ITR2);
% VER=[VER1;VER2] and matching treatment of ITRI
% 20041101;

function [VER,ITRI]=addtwotris(VER1,ITRI1,VER2,ITRI2);

L=size(VER1,1);
VER=[VER1;VER2];
ITRI=[ITRI1;ITRI2+L];

