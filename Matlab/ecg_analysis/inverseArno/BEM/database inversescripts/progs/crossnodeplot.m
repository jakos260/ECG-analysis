% crossnodeplot.m
% script of crossec
% plot and or label nodes for which abs(VER(*.3)-zlevel)<delslab
% plot PNTS for which abs(PNTS(*.3)-zlevel)<del
% AvO; 20050212

   ie3val=get(ie3,'val'); % show/hide nodes (if identified)

 if isempty(plotnodes)==0,
   if ie3val==1,  
       set(hie3,'vis','on');
      else,
          set(hie3,'vis','off');
      end
      
   ie4val=get(ie4,'val');
   if ie4val==1,
       set(hie4,'vis','on');
   else,
       set(hie4,'vis','off');
   end
end

if exist('hie6'),
 ie6val=get(ie6,'val');
   if ie6val==1,
       set(hie6,'vis','on');
   else,
       set(hie6,'vis','off');
   end
end