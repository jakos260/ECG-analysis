% norm3d.m
% function l=norm3d(R)
% norm(s) of 3d vectors; size:(n,3)
  function l=norm3d(R)
  l=sqrt(sum(R'.*R')');


