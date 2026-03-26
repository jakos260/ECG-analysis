% zeromean.m
% shift potentials to zero mean
% function PSI=zeromean(PHI)
function PSI=zeromean(PHI)
[nlds ntims]=size(PHI);
% shift to zero mean
PSI=PHI-ones(nlds,1)*mean(PHI);
