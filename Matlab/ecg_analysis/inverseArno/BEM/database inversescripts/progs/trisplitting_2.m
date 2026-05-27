% file trisplitting_2.m
% [ITRI1,ITRI2]=trisplitting(ITRI,border,hole)
% split triangulated surface ITRI in two sets along a segmentation line:
% border1, which should be (a set of) a closed, non-intersecting contour(s)
% on the surface; each contour spec on border should duplicate the starting point at its end 
% ITRI1 lies to the right, as seen from the outside, while moving along the border
% border2 may be used to specify all nodes demarcating any of internal 'holes',
% if present; such holes should be preferably two node distances away from
% one another as well as from the border.
% A. van Oosterom
% date:051017

function [ITRI1,ITRI2]=trisplitting(ITRI,border,hole)
 
if nargin<3, hole=[];end

nborder=length(border);
if border(1)~=border(nborder),
   border=[border border(1)];
   nborder=nborder-1;
end

ITRI2=ITRI;

% start by identifying triangles that have two nodes on the border 
btris=[];
istart=1;
i=1;
while i<nborder-1,
       %[istart i border(i) border(i+1)]
       labsb=find(sum(ismember(ITRI,border(i:i+1)),2)==2);
       BTRI=ITRI(labsb,:);
       k=0;
       while k<2,
           k=k+1;
           l=find(ismember(BTRI(k,:),border(i))==1);
           if BTRI(k,icyc(l+1,3))==border(i+1),break, end
       end
       btris=[btris labsb(k)];
       i=i+1;
       % treat multiple demarcation lines
       if border(i)==border(istart), i=i+1; istart=i; end
       %pause
end

border(1)=[];

% use only those lying directly to the right of the border 
alltris=find(sum(ismember(ITRI,border),2)==3);
ntr=length(alltris);
keep=[];
for i=1:ntr,
    tris=ITRI(alltris(i),:)
    % check direction
    chain=border;
    chain(~ismember(chain,tris))=[]
    if tricompare(chain(1:3),tris)==1, keep=[keep,alltris(i)]; end
end

btris=unique([keep btris]);
ITRI1=ITRI(btris,:);
set1=unique(ITRI1(:));
ITRI2(btris,:)=[];
set1(ismember(set1,border))=[];
set1(ismember(set1,hole))  =[];

% complete set1 by all nodes to the right of the border
    
    while isempty(ITRI2)==0,
        tris=find(sum(ismember(ITRI2,set1),2)>0);
        
        if isempty(tris)==1, break, end
        FOUND=ITRI2(tris,:);
        ITRI1=[ITRI1;FOUND];
        set1=unique([set1; FOUND(:)]);
        set1(ismember(set1,border))=[];
        set1(ismember(set1,hole))  =[];
        ITRI2(tris,:)=[];
    end
    
    size(ITRI1);
    size(ITRI2);
    
   