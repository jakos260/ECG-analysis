% file trituning_3.m
% decrease wall thickness of atri_pp20_3.tri
% date:20050501

clear

%tris_ra:       1   : 1321           (n= 1321) 
%tris_epi:      1322: 2620           (n= 1299) 
%                                    (n= 2620) 

%tris_ra_endo:  1   :  677           (n=  677)
%tris_ra_epi:    678: 1321           (n=  644)
%                                    (n= 1321)

%tris_la_endo:   1322:2004           (n=  683)
%tris_la_epi:    2005:2620           (n=  616)
%                                    (n= 1299)

%markers for edges etc:
mleftright=	792;
mtricusp=	716;
mivc=		156;
msvc=		102;
mmitral=	735;
mlipv=		423;
mlspv=		474;
mripv=		619;
mrspv=		172;
mfossara= 	797;
mfossala=	362;

    [VERIN,ITRIN]=loadtri('atri_pp20_3.tri');
    ITRIRENDO=ITRIN(1:677,:);
    ITRIREPI=ITRIN(678:1321,:);
    ITRILENDO=ITRIN(1322:2004,:);
    ITRILEPI=ITRIN(2005:2620,:);
    
    VER=VERIN; 
    nver=size(VERIN,1)
    
next=1;
if next==1,
    figure(1)
    clf
    grsw=1;
    ITRI=ITRIRENDO;
    triplot
    tricusp=findedge(VERIN,ITRI,mtricusp,1);
    ivc=findedge(VERIN,ITRI,mivc,1);
    svc=findedge(VERIN,ITRI,msvc,1);
    %figure(2)
    %clf
    %crossec
    pause
end

next=1;
if next==1,
    figure(1)
    clf
    grsw=1;
    ITRI=ITRIREPI;
    triplot
    tricusp=findedge(VERIN,ITRI,mtricusp,1);
    ivc=findedge(VERIN,ITRI,mivc,1);
    svc=findedge(VERIN,ITRI,msvc,1);
    leftright=findedge(VERIN,ITRI,mleftright,1);
    %figure(2)
    %clf
    %crossec
    pause
end

next=1;
if next==1,
    figure(1)
    clf
    grsw=1;
    ITRI=ITRILENDO;
    triplot
    mitral=findedge(VERIN,ITRI,mmitral,1);
    lipv=findedge(VERIN,ITRI,mlipv,1);
    lspv=findedge(VERIN,ITRI,mlspv,1);
    ripv=findedge(VERIN,ITRI,mripv,1);
    rspv=findedge(VERIN,ITRI,mrspv,1);
    %figure(2)
    %clf
    %crossec
    pause
end


next=1;
if next==1,
    figure(1)
    clf
    grsw=1;
    ITRI=ITRILEPI;
    triplot
    mitral=findedge(VERIN,ITRI,mmitral,1);
    lipv=findedge(VERIN,ITRI,mlipv,1);
    lspv=findedge(VERIN,ITRI,mlspv,1);
    ripv=findedge(VERIN,ITRI,mripv,1);
    rspv=findedge(VERIN,ITRI,mrspv,1);
    leftright=findedge(VERIN,ITRI,mleftright,1);
    %figure(2)
    %clf
    %crossec
    pause
end



