% crossub
% script preparing plots of cross-sections of n triangulated
% geometries at zlevel;
% A. van Oosterom; 2006-02-19
% kleur(i) in calling program may be used to specify the color of the crossection of surface i
% ngeoms should specify number of subsurfaces, VERA,ITRIA,VERB,ITRIB, etc
% default ngeoms=1 assumes VERA=VER; ITRA=ITRI;

% usefull settings:

linw=1.5;
ad=get(addlev,'value');

if exist('fixedassen') & ad==0, delete(gca); end

%if exist('fixedassen') axis(fixedassen);  end

hold on

if ~exist('ngeoms'), ngeoms=1; end

% if ~exist('VERA'), VERA=VER; ITRIA=ITRI; ngeoms=1; end
% if  isempty(VERA), VERA=VER; ITRIA=ITRI; ngeoms=1; end

for igeom=1:ngeoms,
    if igeom==1 & igeom~=ngeoms,
        VERSECT=VERA; TRI=ITRIA;
    else,
        VERSECT=VER; TRI=ITRI;
    end
    if igeom==2, VERSECT=VERB; TRI=ITRIB; end
    if igeom==3, VERSECT=VERC; TRI=ITRIC; end
    % include more if required %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('kleur'), kleur(igeom)='b'; end
    
    select=[];
    select=(VERSECT(TRI(:,1),3)>=zlevel|VERSECT(TRI(:,2),3)>=zlevel|VERSECT(TRI(:,3),3)>=zlevel)&...
        (VERSECT(TRI(:,1),3)<=zlevel|VERSECT(TRI(:,2),3)<=zlevel|VERSECT(TRI(:,3),3)<=zlevel);
    CONT=[];
    if isempty(select)==0,
        TRI=TRI(select,:);
        TRI=[TRI TRI(:,1)];
        [ntri,jdum]=size(TRI);
        %TRI represents all triangles of geometry(igeom) that intersect the plane z=zlevel
        ntri=size(TRI);
        for i=1:ntri,
            labs=TRI(i,:);
            labsl=[];
            % find edges in the plane of section
            % find nodes in plane z=zlevel
            labsl=labs(abs(VERSECT(labs,3)-zlevel)<=2*eps);
            if isempty(labsl)==0,
                nnline=length(labsl);
                CONT=[CONT;[VERSECT(labsl,1:2) ones(nnline,1)*zlevel ones(nnline,1)*nnline]];
            end
            if isempty(labsl)==1|length(labsl)<=2,
                % treat remaining posibility: line segment connecting two
                % different edges
                ncross=0;
                xcross=[];
                ycross=[];
                for j=1:3,
                    if (VERSECT(labs(j),3)-zlevel)*(VERSECT(labs(j+1),3)-zlevel)<=0,
                        ncross=ncross+1;
                        alpha=(zlevel-VERSECT(labs(j),3))/(VERSECT(labs(j+1),3)-VERSECT(labs(j),3)+eps);
                        xcross(ncross)=VERSECT(labs(j),1)+alpha*(VERSECT(labs(j+1),1)-VERSECT(labs(j),1));
                        ycross(ncross)=VERSECT(labs(j),2)+alpha*(VERSECT(labs(j+1),2)-VERSECT(labs(j),2));
                    end
                end
                if ncross > 0,
                    CONT=[CONT;[xcross' ycross' ones(ncross,1)*zlevel ones(ncross,1)*ncross]];
                end
            end
        end
        
        if isempty(CONT)==0,
            if igeom==1, CONT1=CONT;
            else,
                CONT2=CONT;
            end
            [ncont jdum]=size(CONT);
            nseg=0;
            ibeg=1;
            while ibeg <= ncont,
                nseg=CONT(ibeg,4);
                iend=ibeg+nseg-1;
                XYZ=CONT(ibeg:iend,1:3);
                lijntje=line(XYZ(1:nseg,2),-XYZ(1:nseg,1));
                set(lijntje,'color',kleur(igeom),'linewidth',linw);
                ibeg=iend+1;
            end
            axis equal
            
            figure(1)
            if igeom==1,
                if exist('nivo1'); delete(nivo1);end
                nivo1=[]; % end
            end
            if igeom==2,
                if exist('nivo2'); delete(nivo2);end
                nivo2=[];
            end
            
            nseg=0;
            ibeg=1;
            while ibeg <= ncont,
                nseg=CONT(ibeg,4);
                iend=ibeg+nseg-1;
                XYZ=CONT(ibeg:iend,1:3);
                if igeom==1,
                    nivo1=[nivo1; line(XYZ(1:nseg,1),XYZ(1:nseg,2),XYZ(1:nseg,3),'color','w')];
                else,
                    nivo2=[nivo2; line(XYZ(1:nseg,1),XYZ(1:nseg,2),XYZ(1:nseg,3),'color','w')];
                end
                ibeg=iend+1;
            end
            
            figure(2)
            plotnodes=[];
            
            % adapt if required
            if igeom==1, plotnodes=find(abs(VERSECT(:,3)-zlevel)<=delslab);end
            
            if isempty(plotnodes)==0,
                ie3val=get(ie3,'val');
                hie3=plot(VERSECT(plotnodes,2),-VERSECT(plotnodes,1),'r*');
                if ie3val==0; set(hie3,'vis','off'), end
                ie4val=get(ie4,'val');
                hie4=text(VERSECT(plotnodes,2),-VERSECT(plotnodes,1),num2str(plotnodes),'fontsize',14);
                if ie4val==0, set(hie4,'vis','off'), end
            end
        end
    end
end

% option for plotting points: SELPNTS
if exist('SELPNTS')
    %figure(2)
    if ~exist('delslab'), delslab=sqrt(eps); end
    slab=find(abs(SELPNTS(:,3)-zlevel)<=delslab);
    ie6val=get(ie6,'val');
    hie6=plot(SELPNTS(slab,2),-SELPNTS(slab,1),'k+','linewidth',2);
    if ie6val==0,set(hie6,'vis','off'),end
end

if exist('focus'), plot(VER(focus,2),-VER(focus,1),'m*'),end
set(gca,'buttondownfcn','movenode')
crossaxes
