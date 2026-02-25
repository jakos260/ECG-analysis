% make_patch.m
% function [VER_new,ITRIN]=make_patch(EDGE, nhmax, shrinks, fig)
% triangulate interior of a closed contour: EDGE
% nhmax: maximum harmonic used in the contour lines
% shrinking factors 1:ns specify the subsequent contours
% patch will be closed if final shrinking factor is set to zero
% for filling up an excisting space in a triangulated surface use: add_patch.m
% A. van Oosterom 2015_03_11

function [VER_new,ITRI_new]=make_patch(EDGE, nhmax, shrinks, fig)
nsf=max(size(shrinks));
shrinks=reverse(sort(shrinks));

if fig~=0,
    figure(fig)
    clf
    plot3(EDGE(:,1),EDGE(:,2),EDGE(:,3),'-*k')
    axis equal
    axis square
    hold on
end

CONT1=EDGE;
VER_new=CONT1;
nver=size(VER_new,1);
LIST=[1 nver];
ITRIN=[];

for is=1:nsf,
    if shrinks(is)==0,break, end
    nc1=size(CONT1,1);
    nrecs=round(nc1/2);
    nrecs=max(nrecs,6); % reconstruction nodes
    CONT2=fourier_contours(CONT1,nhmax,0,nrecs,shrinks(is));
    use_nodes=[VER_new;CONT2];
    nver=size(VER_new,1);
    LIST=[LIST; [LIST(is,2)+1 nver]];
    nc2=size(CONT2,1);
    if fig~=0,
        figure(fig)
        plot3(CONT2([1:nc2 1],1),CONT2([1:nc2 1],2),CONT2([1:nc2 1],3),'k*')
        hold on
    end
    ADD=make_peel(VER_new,LIST(is,1):LIST(is,2),LIST(is+1,1):LIST(is+1,2));
    if fig~=0,
        ntri=size(ADD,1);
        for j=1:ntri;
            quatro=[ADD(j,:) ADD(j,1)];
            plot3(VER_new(quatro,1),VER_new(quatro,2),VER_new(quatro,3))
        end
    end
    ITRIN=[ITRIN;ADD];
    if is<nsf,CONT1=CONT2; end
end

if shrinks(is)==0,
    VER_new=[VER_new;mean(CONT2)];
    nver=size(VER_new,1);
    ADD=make_peel(VER_new,LIST(end,1) : LIST(end,2),nver);
    ITRIN=[ITRIN;ADD];
end
    
 


