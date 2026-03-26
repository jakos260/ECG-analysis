% normalize_3d.m
% function R=normalize_3d(R)
% normalize n 3D vectors; size:(n,3)
% A. van Oosterom; 2014_08_17

  function R=normalize_3d(R);
  l=sqrt(sum(R'.*R')');
  R=R./(l*ones(1,3));
  
  


