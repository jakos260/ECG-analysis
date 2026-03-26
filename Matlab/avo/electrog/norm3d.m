% norm3d.m
% function l=norm3d(R)
% norm of 3d vectors
  function l=norm3d(R)
l=sqrt(sum(R'.*R')');


