echo off
clear
TAB=loadasci('sterfte.tab');
tvals=TAB(1:77,1);
hvals=TAB(1:77,2);
par=  [1.
         100
	  6];
OPTIONS=foptions;	  
[par,OPTIONS]=curvefit('hazfun',par,tvals,hvals,OPTIONS,'hazgrad')


