%testprogr voor mouse control
clf
x=1:100;
y=sin(x/100*2*pi);
plot(x,y)
h=line([50 50],[0 .1],'color',[1 0 0])
get(h)
%figure
%b=uicontrol('Style','pushbutton', ...
%            'Units','normalized',...
%	    'Position',[.5 .5 .01 .01]);
%s='set(b,''Position'',[.8*rand .9*rand .01 .01])';
%%eval(s)
%set(b,'Callback',s)
%clf
%figure
hold on
i=1;
while i~=0,
[x,y,button]=ginput(1)
plot([x x],[y y],'r*')
if button==2, break;end
end
get(gtext('+'))

	    