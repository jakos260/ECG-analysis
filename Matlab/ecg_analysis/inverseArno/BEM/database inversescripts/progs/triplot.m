% file triplot.m
%  specify: VER;ITRI;VALS
%  defaults set: column=1,
%                node=min(unique(ITRI:))); (usually = 1)
%                CMAP=loadmat('avopot.mcm');
%                symcolor=0; symcolor==1 produces colorbar symmetric around zero
%                newboxsize; default: max(max(abs(VER)))
%                default: zebra=0; zebra>nzebra produces nzebra stripes (isofunction
%                lines) zebra<0 produces isofunction lines at abs(zebra)
%                       intervals
%                nobar=0;  nobar=1 suppresses colorbar display
%                if exist('SELPNTS'); these will be shown when using
%                crossec


% A. van Oosterom;
% see also: triplot_contour

% 20080727

%  hs=patch; 'ButtondownFcn','nodesel';
%                                                           defaults:
%  ui1r1=uicontrol('style','radio' ; Callback:  rot     ;
%  ui1r2=uicontrol('style','radio' ; Callback:  gridsw  ;   grsw=1;
%  ui1r3=uicontrol('style','radio' ; Callback:  nodesw  ;   nosw=0;
%  ui1r4=uicontrol('style','radio' ; Callback:  shofunsw;   funsw=0;
%  ui1r5=uicontrol('style','radio' ; Callback:  lightsw;
%  ui1r6=uicontrol('style','radio' ; Callback:  c_thru  ;   cthruval=0;

%  ui7  =uicontrol('style','popup' ; Callback:  setview ;
%  ui100=uicontrol('style','pushb' ;'Callback','close figure 1'
%  ui34 =uicontrol('style','pushb' ;'Callback','node2col');
%  sltri=uicontrol('style','slider';'CalllBack','displcol'

%  ui2a =uicontrol('style','text';string','node selection disabled'
%  ui3  =uicontrol('style','text';'string';'selected node:'
%  ui3a =uicontrol('style','text';'string','value at selected node '
%  ui4  =uicontrol('style','text';'string',sprintf('%3.3f',fun(node)
%  ui5  =uicontrol('style','text';'string',num2str(node)
%  ui13 =uicontrol('style','text';'string','selected triangle: '
%  ui14 =uicontrol('style','text';'string',num2str(itri))
%  ui15 =uicontrol('style','text';'string','function matrix VALS=[]'
%  ui17 =uicontrol('style','text';'string',num2str(column))

dim=size(ITRI);  ntri=dim(1);
dim=size(VER);   nver=dim(1);

% set defaults%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('newboxsize'),
    boxsize=max(norm3d(VER));
else,
    boxsize=newboxsize;;
end

nosw=0;
funsw=0;
fun=zeros(nver,1);
ndat=1;
maxcoll=1;
if exist('VALS'),
    [ndat maxcoll]=size(VALS);
    % show function values
    ui1r4=uicontrol('style','radio');
    uibox1r4=[.0 .85 .09 .05];
    set(ui1r4,'units','norm','position',uibox1r4','string','vals',...
        'val',0,'Callback','shofunsw')
end
if ~exist('nobar') nobar=0; end  % brute force for suppressing colorbar
if ~exist('lsw') lsw=1; end
if ~exist('zebra')  zebra=0; end
if ~exist('node')   node=1 ; end
if ismember(node,unique(ITRI(:)))==0, node=min(unique(ITRI(:))); end
if ~exist('column') column=1;end
if ~exist('symcolor'), symcolor=0; ,end % symcols==1 produces colorbar symmetric around zero
if ~exist('lsw'),lsw=1; end
if  column > maxcoll, column=maxcoll; end
if  column < 1, column=1; end
if exist('VALS') & isempty(VALS)==0,
    if node > nver, nver=1; end
    fun=VALS(:,column);
    set(ui1r4,'val',1); funsw=1;
end
if ~exist('iview'), iview=1; end
if ~exist('azim'), azim=90; elev=0; end


if ~exist('cmap'),
    cmap='avojet.mcm';
    cmapnow=cmap;
    CMAP=loadmat(cmap);
end
if ~exist('cmapnow'),cmapnow='xxx';end
if cmap(1:3)~=cmapnow(1:3),
    CMAP=loadmat(cmap);
    cmapnow=cmap;
end

colormap(CMAP);

if exist('zebra')
    if zebra ~= 0,
        stripes=[];
        ncolors=size(CMAP,1);
        if zebra < 0,
            isostep=abs(zebra);
            mafun=max(fun); mifun=min(fun);
            if symcolor==1,  mafun=max(abs([fun;eps]));mifun=-mafun; end
            spanfun=mafun-mifun+eps;
            beglev=sign(mifun)*floor(abs(mifun)/isostep)*isostep;
            nlev=0;
            test=beglev;
            while test <= mafun,
                nlev=nlev+1;
                stripes(nlev)=min(max(1,round(0.5+(test-mifun)*ncolors/spanfun)),ncolors);
                test=test+isostep;
            end
        else,
            stripes=[1:zebra];
            stripes=round(ncolors*stripes/(zebra+1));
        end
        CMAP(stripes,1:3)=0;
    end
    colormap(CMAP)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ax1a=axes('pos',[.15 .15 .8 .8]);
set(ax1a,'xlim',[-boxsize boxsize],'ylim',[ -boxsize boxsize],'zlim',[-boxsize boxsize]);


%ax1a=axes('pos',[.0 .0 1 1]);
%set(ax1a,'xlim',[-boxsize boxsize],'ylim',[ -boxsize boxsize],'zlim',[-boxsize boxsize]);

axis vis3d
axis off
view(azim,elev)

% colorplot of function values involving linear interpolation
hs=patch('Vertices',VER,'Faces',ITRI,...
    'FaceVertexCData',fun,'FaceColor','inter');
set(hs,'ButtondownFcn','nodesel');
set(hs,'EdgeColor',[.5 .5 .5],'Marker','.','MarkerEdgeColor','none');

view(azim,elev)
lightpos=[cos((azim-90)/180*pi) sin((azim-90)/180*pi) 0];
hlight=light('position',lightpos,'style','infinite');

if symcolor==1,
    % symmetric color bar around zero
    extrmap=max(abs([fun;eps]));
    set(ax1a,'Clim',[-extrmap extrmap]);
end

if nobar==0;
    hcbar=colorbar;
    hcbarpos=[.9 .3 .03 .6];
    set(hcbar,'position',hcbarpos)
end

% highlight selected node
ht=text(1.03*VER(node,1),1.03*VER(node,2),1.03*VER(node,3),'O','color','w');
set(ht,'HorizontalAlignment','center','FontSize',16,'Fontweight','bold')

% light switch
ui1r5=uicontrol('style','radio');
uibox1r5=[.0 .8 .09 .05];
set(ui1r5,'units','norm','position',uibox1r5','string','light',...
    'val',lsw,'Callback','lightsw')
lightsw

% enable/disenable rotating the geometry
ui1r1=uicontrol('style','radiobutton');
uibox1r1=[.0 .55 .12 .05];
set(ui1r1,'units','norm','position',uibox1r1','string','rotate','Callback','rot')
set(ui1r1,'value',0)
rotate3d off;

% warning that rotating the geometry prevents node selection
uibox2a=[.033 .32 .09 .16];
ui2a=uicontrol('style','text','units','norm','position',uibox2a,'string','node selection disabled');
set(ui2a,'visible','off','backgroundcolor','w')

% text: select view
uibox7a=[.0 .65 .12 .05];
ui7a=uicontrol('style','text','units','norm','position',uibox7a,'string',...
    'set view');

% select view
ui7=uicontrol('style','popup');
uibox7=[.0 .6 .12 .05];

set(ui7,'units','norm','position',uibox7,'string',...
    'anterior|left ant|left|left post|posterior|right post|right|right ant|superior|inferior',...
    'val',iview,'Callback','setview');
set(ui7,'foregroundcolor','b');

% text: selected node
ui3=uicontrol('style','text');
uibox3=[.02 .12 .14 .05];
set(ui3,'units','norm','position',uibox3',...
    'string','selected node: ')

% signal: selected node
ui5=uicontrol('style','text');
uibox5=[0.15 .12 .06 .05];
set(ui5,'units','norm','position',uibox5,'string',num2str(node))
set(ui5,'foreground','b','fontsize',12)

% text: value at selected node
ui3a=uicontrol('style','text');
uibox3a=[.02 .05 .35 .05];
set(ui3a,'units','norm','position',uibox3a',...
    'string',' value at selected node ','fontsize',10)

% signal: value at selected node
ui4=uicontrol('style','text');
uibox4=[.4 .05 .12 .05];
set(ui4,'units','norm','position',uibox4,'string',sprintf('%3.3f',fun(node)));
set(ui4,'Foregroundcolor','b','fontsize',12)


if funsw==0, set([ui3a;ui4],'vis','off'); end

viscene=1; % (used in ecgsim)

% show grid
if ~exist('grsw'), grsw=0;,end
ui1r2=uicontrol('style','radio');
uibox1r2=[.0 .95 .09 .05];
set(ui1r2,'units','norm','position',uibox1r2','string','grid',...
    'val',grsw,'Callback','gridsw')
gridsw;

% show nodes
if ~exist('nodesval'), nodeswval=0;, end
ui1r3=uicontrol('style','radio');
uibox1r3=[.0 .9 .09 .05];
set(ui1r3,'units','norm','position',uibox1r3','string','nodes',...
    'val',nodeswval,'Callback','nodesw')
nodesw

% set tranparency
if ~exist('cthruval'), cthruval=0;, end
ui1r6=uicontrol('style','radio');
uibox1r6=[.0 .75 .09 .05];
set(ui1r6,'units','norm','position',uibox1r6','string','c_thru',...
    'val',cthruval,'Callback','cthrusw')
cthrusw

% exit button
ui100=uicontrol('style','pushbutton');
uibox100=[.94 .01 .05 .05];
set(ui100,'units','norm','position',uibox100,'string','exit',...
    'Callback','close figure 1');

% text: selected triangle
ui13=uicontrol('style','text');
uibox13=[.62 .05 .14 .08];
set(ui13,'units','norm','position',uibox13',...
    'string','selected triangle: ')

% signal: selected triangle
ui14=uicontrol('style','text');
uibox14=[0.76 .05 .06 .05];
set(ui14,'units','norm','position',uibox14,'string',num2str(1))
set(ui14,'foreground','b','fontsize',12)
%set([ui13;ui14],'vis','off')

% text: no function specified
ui15=uicontrol('style','text');
uibox15=[.02 .0 .35 .05];
set(ui15,'units','norm','position',uibox15,...
    'string','no function specified')
if funsw==0,set(ui15,'vis','off'); end

% text: 'column'
ui16=uicontrol('style','text');
uibox16=[.7 .95 .08 .05];
set(ui16,'units','norm','position',uibox16,...
    'string','column: ')

% text: column
ui17=uicontrol('style','text');
uibox17=[.78 .95 .1 .05];
set(ui17,'units','norm','position',uibox17,'string',num2str(column))
set(ui17,'foreground','b','fontsize',12)

ui34=uicontrol('style','pushbutton');
uibox34=[.02 .2 .12 .05];
set(ui34,'units','norm','position',uibox34,'string','node2col',...
    'Callback','node2col');
if ndat > maxcoll, set(ui34,'vis','off'); end

sltribox=[.3 .97 .4 .03];
if maxcoll > 1,
    sltri=uicontrol('style','slider','units','norm','position',sltribox,...
        'sliderstep',[1,1]/(maxcoll-1),'min',1,'max',maxcoll,'Value',column,'CallBack','displcol');
end

if exist('VALS') & funsw==1, set(ui15,'vis','off');end

if ~exist('VALS')|isempty(VALS)==1, set([ui16;ui17],'vis','off');
else,              set([ui16;ui17],'vis','on');
end

if exist('SELPNTS')
    hold on
    plot3(SELPNTS(:,1),SELPNTS(:,2),SELPNTS(:,3),'+','linewidth',2)
end



