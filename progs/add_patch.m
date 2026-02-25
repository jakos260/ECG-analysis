% add_patch.m
% function [VER,ITRI]=add_patch(VER,ITRI,edge, nhmax, shrinks, fig)
% triangulate interior of a closed contour: VER(edge,:)
% nhmax: maximum harmonic used in the contour lines
% shrinking factors 1:ns specify the subsequent contours
% new patch will be closed if final shrinking factor is set to zero
% output includes the added patch
% fig~=0 shows the consruction
% 
% A. van Oosterom 2015_03_11

function [VER,ITRI]=add_patch(VER,ITRI,edge, nhmax, shrinks, fig) 

nsf=max(size(shrinks));
shrinks=reverse(sort(shrinks));

if fig~=0,
    figure(fig)
    clf
    plot3(VER(edge,1),VER(edge,2),VER(edge,3),'-*k')
    axis equal
    axis square
    hold on
end

a=edge';
for is=1:nsf,
    if shrinks(is)==0,break, end
    na=size(a,1);
    nrecs=round(na/2);
    nrecons=max(nrecs,6); % reconstruction nodes
    
    %NEWC=fourier_contours(CONT,maximum number hamonics,fig,number reconstructed nodes,contraction_fractions)
    CONT2=fourier_contours(VER(a,:),nhmax,fig,nrecons,shrinks(is));
    nc2=size(CONT2,1)
    
    if fig~=0,
        figure(fig)
        plot3(CONT2([1:nc2 1],1),CONT2([1:nc2 1],2),CONT2([1:nc2 1],3),'k*')
        hold on
    end
   
    nver=size(VER,1)
    VER=[VER; CONT2];
    nver2=nver+nc2;
    b=(nver+1 : nver+nc2)';
    ADD=make_peel(VER,a,b);
    
    if fig~=0,
        ntri=size(ADD,1);
        for j=1:ntri;
            quatro=[ADD(j,:) ADD(j,1)];
            plot3(VER(quatro,1),VER(quatro,2),VER(quatro,3))
        end
    end
    
    ITRI=[ITRI;ADD];
    if is<nsf,a=b; end
end

if shrinks(is)==0,
    VER=[VER;mean(VER(b,:))];
    nver=size(VER,1);
    ADD=make_peel(VER,b,nver);
    ITRI=[ITRI;ADD];
end
    
 





