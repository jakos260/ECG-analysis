function AxisLegend(varargin)

t=varargin{1};
y=varargin{2};
dT=varargin{3};
dmi=varargin{4};
if length(varargin) > 4 
	doline=varargin{5};
else
	doline=1;
end

delete(findobj(get(gca,'Children'),'Tag','ledendL1'));
delete(findobj(get(gca,'Children'),'Tag','ledendL2'));
delete(findobj(get(gca,'Children'),'Tag','ledendL3'));
delete(findobj(get(gca,'Children'),'Tag','ledendT1'));
delete(findobj(get(gca,'Children'),'Tag','ledendT2'));

maxT=max(t);
mi=min(y);
axis off
if doline
	line([0 maxT], [0 0],'Linewidth',1,'Color','k','Tag','ledendL1');
end
mi=0.85*mi;
% dT=50;
dmi=abs(dmi);
line([ maxT-dT maxT], [mi mi],'Linewidth',2,'Color','k','Tag','ledendL2');
line([ maxT-dT maxT-dT], [mi mi+dmi],'Linewidth',2,'Color','k','Tag','ledendL3');
if dT>1000
	text(maxT-dT*0.65,mi*1.08, [num2str(dT/1000) ' s'],'FontUnits','normalized','Fontname','verdana','FontSize',0.04,'Tag','ledendT1');
elseif dT>0.1
	text(maxT-dT*0.65,mi*1.08, [num2str(dT) ' ms'],'FontUnits','normalized','Fontname','verdana','FontSize',0.04,'Tag','ledendT1');
elseif dT > .1e-3
	text(maxT-dT*0.65,mi*1.08, [num2str(dT) ' \mus'],'FontUnits','normalized','Fontname','verdana','FontSize',0.04,'Tag','ledendT1');
end
if dmi>.1
	text(maxT-dT*1.25,mi+dmi*0.3,[num2str(dmi) ' mV'],'FontUnits','normalized','Fontname','verdana','FontSize',0.04,'Rotation',90,'Tag','ledendT2');
elseif dmi>.1e-3
	text(maxT-dT*1.25,mi+dmi*0.3,[num2str(1000*dmi) ' \muV'],'FontUnits','normalized','Fontname','verdana','FontSize',0.04,'Rotation',90,'Tag','ledendT2');
end
% axis([min(t) maxT min(y) max(y)])