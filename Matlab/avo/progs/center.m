% file center.m
% matlab script supporting triplot
% shift  VER to zero center of gravity of the nodes
[nver idum]=size(VER);
VER=VER-ones(nver,1)*mean(VER);
