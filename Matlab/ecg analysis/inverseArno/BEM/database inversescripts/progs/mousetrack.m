function z=mousetrack(option)

% MOUSETRACK Demonstrates how to track the mouse position
% key tool: z=get(gca,'currentpoint')


%  Press the left mouse button in the axis and move arround while
%  pressed

% Change History :
% Date		Time		Prog	Note
% 07-Jun-2000	12:04 PM	ATC	Created under
% MATLAB 5.3.1.29215a (R11.1)

% ATC = Ali Taylan Cemgil,
% SNN - University of Nijmegen, Department of Medical Physics and Biophysics 
% e-mail : cemgil@mbfys.kun.nl 

%persistent x y uic ncalls
%persistent x y 

if ~exist('option'), option = 'start'; end;
switch option,
    case 'start',
  %clf
  %uic = uicontrol('style','edit','pos',[10 10 100 20]);
  set(gca,'buttondownfcn',[mfilename '(''down'')']);
  % x = 0; y = 0; ncalls=0;

 case 'down',
  set(gcf,'windowbuttonmotionfcn', [mfilename '(''drag'')']);
  set(gcf,'windowbuttonupfcn', [mfilename '(''release'')']);
  %z = get(gca, 'currentpoint'); 
  %x = z(1,1); y= z(1,2);  
  %set(uic, 'string', sprintf('%3.3f - %3.3f', x, y));
  
 case 'drag',
  %z = get(gca, 'currentpoint'); 
  %x = z(1,1); y= z(1,2);  
  %set(uic, 'string', sprintf('%3.3f - %3.3f', x, y));

 case 'release',
  z = get(gca, 'currentpoint');
  z = get(gca, 'currentpoint'); 
  %x = z(1,1); y= z(1,2);  
  % x = z(1,1); y= z(1,2);     
  % set(gcf, 'userdata', [x y])
 set(gcf,'windowbuttonmotionfcn', ['']);
 set(gcf,'windowbuttonupfcn', ['']);
 %ncalls=ncalls+1;
 
end;

 