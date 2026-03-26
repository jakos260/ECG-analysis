% tri_peeling.m
% 20090428
% peeling  triangulated surface, starting from the edge: "seed"

% under construction
% clear all
%[VER,ITRI]=icosasub(1,0);
% %seed=7;
% % %[VER,ITRI]=loadtri('amtor.tri');
% % [VER,ITRI]=loadtri('../gu/guhar.tri');

% focus=1063;
% ITRIN=ITRI;
% VERIN=VER;
% DIS=edge_dist(ITRIA,focus);
% selnodes=DIS(DIS(:,2)<22,1);
% ITRI=ITRI(ismember(ITRI(:,1),selnodes)&ismember(ITRI(:,2),selnodes)& ...
%           ismember(ITRI(:,3),selnodes),:);
% VER=10*(VERIN-ones(nver,1)*VER(focus,:));
% VALS=DIS(:,2);
%
% figure(1)
% clf
% zebra=-1;
% triplot_contour
% pause

ITRIS=ITRI;

ntri = size(ITRI,1);

EDGES=[[ITRI(:,1:2) (1:ntri)']; [ITRI(:,2:3) (1:ntri)']; [ITRI(:,[3 1]) (1:ntri)']];
EDGES=sortrows(EDGES,3);

%[ADJ,DN]=graphdist(ITRI);        % minimum number of edges to be passed
%[ADJ,DG]=graphdist(ITRI,VER,1);  % minimum distance along edges to be passed

VERIN=VER;
% form edge pairs;
% prepare
PAIRS=zeros(3*ntri,4);  %[triangle_index edge_label  (1  2 or 3)]
for i=1:ntri,
    PAIRS((i-1)*3+1,1:2)  = [i 1];
    PAIRS((i-1)*3+2,1:2)  = [i 2];
    PAIRS((i-1)*3+3,1:2)  = [i 3];
end
% matchmaking
nedg=3*ntri;

for i=1:nedg,
    itri=PAIRS(i,1);
    elab=PAIRS(i,2);
    trits=[elab icyc(elab+1,3) icyc(elab+2,3)];  % cyclic permutations of [1 2 3]
    edg=ITRIS(itri,[trits(1) trits(2) ]);
    ip=find(EDGES(:,1)==ones(nedg,1)*edg(2)& EDGES(:,2)==ones(nedg,1)*edg(1));
    
    if isempty(ip)==0,
        itrib=EDGES(ip,3);
        EDGES(ip,:);
        PAIRS(i,3:4)=[itrib  PAIRS(ip,2)];
        PAIRS(ip,3:4)=[itri elab];
    end
end

% pause
% singles not permitted;
PAIRS(PAIRS(:,4)==0,:)=[];
np=size(PAIRS,1);

% start peeling first one is a "pit"
peel=352;  %nb: a triangle index
npl=1;
ep=find(PAIRS(:,1)==ones(np,1)*peel(npl));
spiral=zeros(ntri,1);
spiral(1)=352;
TREAT=PAIRS(ep,:)
itop=ep(1)


while npl<ntri,
    
    while isempty(TREAT)==0,
        npl=npl+1;
        newtri=TREAT(1,3)
        spiral(npl)=newtri;
        spiral(1:npl)
        %process new triangle
        % find its available edges
        topref=PAIRS(itop,:)
        ep=find(PAIRS(:,1)==newtri*ones(np,1));
        match=ep(PAIRS(ep,3)==topref(1))
        PAIRS(match,:)
        usedge=icyc(PAIRS(match,2)-1,3)
        ep(PAIRS(ep,3)==topref(1))=[];
        PAIRS(ep,:)
        k=find(ep(PAIRS(ep,2)==usedge))
        TREAT(1,:)=[];
        if isempty(k)==0,
            topper=PAIRS(ep(k),:)
            TREAT=[topper; TREAT];
            PAIRS([itop; ep(k)],:)=[];
            itop=ep(k);
            np=np-2;
            ep(k)=[];
        end
        if isempty(ep)==0,
            TREAT=[TREAT; PAIRS(ep,:)];
        end
        TREAT
        pause
        
        figure(1)
        clf
        focus=930;
        VER=20*(VERIN-ones(nver,1)*VERIN(focus,:));
        ITRI=ITRIS(spiral(1:npl),:);
        VALS=VER;
        newboxsize=max(max(abs(VER(unique(ITRIS(:)),:))));
        gridsw=1;
        node=seed(1);
        triplot
        pause
    end
end







%
%        inb=icyc(PAIRS(ep(1),4)+2,3)
%
%        [itop ep(1)]
%        PAIRS([itop ep(1)],:)
%        PAIRS([itop ep(1)],:)=[];
%
%        np=np-2;
%        TREAT(1,:)=[];
%
%
%        k=find(PAIRS(ep,4)==inb);
%        if isempty(k)==0;
%           treat1=PAIRS(PAIRS(:,1)==PAIRS(ep(k),3),:)
%           TREAT=[treat1; PAIRS(ep,:)]
%           pause
%           itop=ep(k);
%           ep(k)=[];
%        end
%        if isempty(ep)==0,
%            TREAT=[TREAT;PAIRS(ep,:)];
%        end
%        ntreat=size(TREAT,1)
%        TREAT(1:ntreat,:)
%        pause
%
%
%
%
%     pause

