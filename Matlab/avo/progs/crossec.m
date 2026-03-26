% crossec.m
% script for plotting cross-sections of ngeom triangulated
% geometries at level: zlevel; shown in in figure(1)
% call should specify figure(2)
% default: ngeoms==1; VER and ITRI required
% for ngeoms>1
% kleur(i) in call may specify the color of the crossection of surface i
% default: blue
% ngeoms should specify the number of subsurfaces
% in that case: specify:  VERA,ITRIA,    VERB,ITRIB,      etc.
% currently: 0< n <= 3 supported in script crossub.m
% A. van Oosterom; 200706-23

zoom=1;
clf

% set defaults
if  exist('geom')==0,geom='geom';end
if ~exist('zlevel'),  zlevel=mean(VER(:,3));end
if ~exist('window'),  window=1.1*[min(VER(:,2)) max(VER(:,2)) -max(VER(:,1)) -min(VER(:,1))];end
if ~exist('zincr'),   zincr=(max(VER(:,3))-min(VER(:,3)))/100; end
if ~exist('slmax'),   slmax=ceil(max(VER(:,3))/zincr)*zincr; end
if ~exist('slmin'),   slmin=floor(min(VER(:,3))/zincr)*zincr; end
if ~exist('delslab'), delslab=sqrt(eps); end

% a selection of points (SELPOINTS) to be may be pl

%hold on
% ui's ordered clockwise as viewed in the window 

ie7box=[.92 .9 .08 .05]; fixwin=0;
ie7=uicontrol('style','radio');
set(ie7,'units','norm','position',ie7box',...
'string','fix-win','val',fixwin);

addbox=[.92 .82 .08 .05];
addlev=uicontrol('style','radio','units','norm','position',addbox);
set(addlev,'Value',0,'string','add');

boxt1=[.92 .7 .08 .05];
uit1=uicontrol('style','text','units','norm','position',boxt1);
set(uit1,'string','zlevel');
%set(uit1,'Value',zlevel,'string','zlevel');

boxx=[.94 .29 .03 .41];
sllevel=uicontrol('style','slider','units','norm','position',boxx);
set(sllevel,'min',slmin,'max',slmax,'Value',zlevel);
sstep=zincr/(slmax-slmin); 
slstep=[sstep .1]; 
set(sllevel,'Sliderstep',slstep);
set(sllevel,'callback','crosslev1');

ie1box=[.84 .24 .16 .05];
ie1=uicontrol('style','edit');
set(ie1,'units','norm','position',ie1box',...
'string',num2str(zlevel),'callback','crosslev2');

ie5box=[.84 .0 .16 .05]; 
ie5=uicontrol('style','edit'); % z-increment
set(ie5,'units','norm','position',ie5box',...
'string',num2str(zincr),'callback','crosszincr');

boxt3=[.77 0. .07 .05];
uit3=uicontrol('style','text','units','norm','position',boxt3);
set(uit3,'string','zincr');

boxt2=[.16 0. .1 .05];
uit2=uicontrol('style','text','units','norm','position',boxt2);
set(uit2,'Value',delslab,'string','delslab');

ie1abox=[.00 .0 .16 .05];
ie1a=uicontrol('style','edit');
set(ie1a,'units','norm','position',ie1abox',...
'string',num2str(delslab),'callback','crossdelta');

ie4box=[.0 .06 .09 .05]; 
if ~exist('ie4val'), ie4val=0; end
ie4=uicontrol('style','radio');
% set display node-labels on/off
set(ie4,'units','norm','position',ie4box',...
'string','labels','val',ie4val,'callback','crossub');

ie3box=[.0 .12 .09 .05];
if ~exist('ie3val'), ie3val=0; end
ie3=uicontrol('style','radio');

% set display nodes on/off
set(ie3,'units','norm','position',ie3box',...
'string','nodes','val',ie3val,'callback','crossub');

ie6box=[.0 .18 .09 .05]; ie6val=0;
ie6=uicontrol('style','radio');
% set display extra points on/off
set(ie6,'units','norm','position',ie6box',...
'string','pnts','val',ie6val,'callback','crossub');

plotnodes=find(abs(VER(:,3)-zlevel)<=delslab);
nivo=[]; nivo1=[]; nivo2=[];
crossub;

