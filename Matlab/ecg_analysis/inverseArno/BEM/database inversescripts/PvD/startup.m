%c change the working directory
%cd h:\KUN' 'stuff\

set(0,'DefaultFigurePaperUnits','centimeters');
set(0,'DefaultFigurePaperType','A4');

%set(0,'DefaultFigurePaperOrientation','portrait');si=[21,29.7]; 
% set(0,'DefaultFigurePaperOrientation','landscape');si=[29.7,21];
% set(0,'DefaultFigurePaperSize',si);
% set(0,'DefaultFigurePaperPosition',[1,1,si(1)-1,si(2)-1]);
%set(0,'DefaultFigurePaperPositionMode','auto');
set(0,'DefaultTextFontName','Arial');
set(0,'DefaultTextFontsize',9);
set(0,'DefaultTextErasemode','xor');

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultAxesFontsize',9)

pathdef

% oostep1

% Include path to ps2pdf and Qtriplot for the Mac
path1 = getenv('PATH');
path1 = ['/Applications/MacGhostViewFolder/bin:' '/Applications/qtriplot.app/Contents/MacOS:' path1];
setenv('PATH', path1);
clear path1
format short
