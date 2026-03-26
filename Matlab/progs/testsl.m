% test Surface Laplacian 
% test intripol
clear
file='bol20.tri';
[VER,ITRI]=loadtri(file);
%VER=VER/10.
%SL=surflapl(VER,ITRI,0);
nodes=[12 1 6 5];
T=intripol(VER,ITRI,nodes)
x2=[1; 0; 0; -1;];
x=T*x2