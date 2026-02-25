% make_cube.m
% function [VER, ITRI]=make_square;
% generates basic triangular mesh for a unit cube:
% centered at origin 
% 2015_01_20; A. van Oosterom

function [VER, ITRI]=make_cube

VER=[ 0.5 -0.5 -0.5; 
      0.5  0.5 -0.5;
     -0.5  0.5 -0.5; 
     -0.5 -0.5 -0.5;
      0.5 -0.5  0.5; 
      0.5  0.5  0.5;
     -0.5  0.5  0.5; 
     -0.5 -0.5  0.5;];
 
 ITRI=[1 5 6; 
       1 6 2;
       2 6 7;
       2 7 3;
       3 7 8;
       3 8 4;
       4 8 5; 
       4 5 1;
       1 2 4;
       2 3 4;
       5 8 6;
       6 8 7;];
   
       
       
     
