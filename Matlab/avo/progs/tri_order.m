
% tri_flipping.m
%[ITRI, flipped]=tri_order(ITRI,seed) 
% check consistent order of trangulation indices; flip indices if required
% seed: index of any triangle having the desired order

% 2009_05_11; A. van Oosterom


function [ITRI, flipped]=tri_flipping(ITRI,seed)

ntri=size(ITRI,1);
seed=seed(1);
LIST1=[ITRI(seed,:) seed 0];
LIST2=[ITRI (1:ntri)' zeros(ntri,1)];
LIST2(seed,:)=[];
flipped=[];
correct=seed;
while isempty(LIST2)==0,
    % search triangles among LIST2 that have an edge in common with LIST1
    n1=size(LIST1,1);
    n2=size(LIST2,1);
    FOUND=[];
    for i=1:n2,
        FOUND=LIST1(sum(ismember(LIST1(:,1:3),ones(n1,1)*LIST2(i,1:3)),2)>1,:);
        if isempty(FOUND)==0, break, end
    end
    if isempty(FOUND),
       'incorrect geometry; single surface, carrying triangle: seed, required'
    end
    a=FOUND(1,:);
    b=LIST2(i,:);
    itri2=LIST2(i,4);
    % flip orientation in ITRI, if required
    for rev=1:3,
        tst=sum(b(1:3)==revolute(a(1:3),rev-1),2);
        if tst==2,
        % flip triangle indices triangle i of LIST2
            flipped=[flipped;itri2];
            ITRI(itri2,1:2)=ITRI(itri2,[2 1]); 
            LIST2(i,1:2)=LIST2(i,[2 1]);
            break;
        end
    end
    LIST1=[LIST1;LIST2(i,:)];
    tris=[FOUND(:,4); LIST2(i,4)];
    k=ismember(LIST1(:,4),tris);
    LIST1(k,5)=LIST1(k,5)+1;
    nf=size(FOUND,1);
    if nf>1,
        nl1=size(LIST1,1);
        LIST1(nl1,5)=LIST1(nl1,5)+nf-1;
    end
    LIST1(LIST1(:,5)>2,:)=[];
    LIST2(i,:)=[];
%     pause
    if isempty(LIST2), break, end
    
    
end










