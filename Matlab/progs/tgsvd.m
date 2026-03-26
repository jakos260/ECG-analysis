% tgsvd.m
% program for interactive study of 
% the GSVD idea
clear
A=[.2 .3 .5
   .6 .3 .1];
AA=[.2 .3 .5
   .6 .3 .1
   .0 .0 .0
   .0 .0 .0
   .0 .0 .0];
B=[-2  1  1
    1 -2  1
    1  1 -2];
BB=[0  0  0
    0  0  0
   -2  1  1
    1 -2  1
    1  1 -2];
C=[A',B'];
C=C';
[UA,SA,VA]=svd(A)
[UB,SB,VB]=svd(B)
[UC,SC,VC]=svd(C)

[UAA,SAA,VAA]=svd(AA)
[UBB,SBB,VBB]=svd(BB)
