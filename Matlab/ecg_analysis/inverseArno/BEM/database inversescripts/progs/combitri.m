% combitri.m
% function [VER,ITRI]=combitri(TRIFILES)
% note: each row, containing filemame triangle file, must have equal length
% 2021107

function [VER,ITRI]=combitri(TRIFILES)
[nfiles nchar]=size(TRIFILES);
VER=[]; ITRI=[];
for i=1:nfiles,
[VERi,ITRIi]=loadtri(TRIFILES(i,:));
l=length(VER);
VER=[VER;VERi];
ITRI=[ITRI;ITRIi+l];
[size(VER) size(ITRI)];
end
