% tricompare.m
% function alike=tricompare(ITRI1,ITRI2)
% compare triangles on the basis of their indices; cyclic permutations are 
% taken to reflect identity 
% triangles having the same sense:   alike=1;
% different triangles:               alike=0;
% triangles having reversed sense:   alike=-1;
% 20050421; A. van Oosterom

function alike=tricompare(ITRI1,ITRI2)
ITRI1
ntri1=size(ITRI1,1)
ntri2=size(ITRI2,1)
alike=zeros(min(ntri1,ntri2),1);
ITRI1=[ITRI1 ITRI1(:,1:2)];
[mi imi]=min(ITRI1');

SEL=[(1:ntri1)' imi'];
sel=SEL(:,1)+(SEL(:,2)-1)*ntri1;
ITRI1=[ITRI1(sel) ITRI1(sel+ntri1) ITRI1(sel+2*ntri1)];

ITRI2=[ITRI2 ITRI2(:,1:2)];
[mi imi]=min(ITRI2');

SEL=[(1:ntri2)' imi'];
sel=SEL(:,1)+(SEL(:,2)-1)*ntri2;
ITRI2=[ITRI2(sel) ITRI2(sel+ntri1) ITRI2(sel+2*ntri1)];

alike(sum(ITRI1==ITRI2,2)==3)=1;
alike(sum(ITRI1(:,2:3)==ITRI2(:,[3 2]),2)==2)=-1;

            
