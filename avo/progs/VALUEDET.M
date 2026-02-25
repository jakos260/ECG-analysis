% VAlUEDET shows window, mean + sd, maximum and minimum values.
%        of the window set by lselect
%        These values are written in the variable Values
%        

%function detminmax

global Maximum Minimum

%number of axes children
ax_ch = get(gca,'Children');
si_ax = size(ax_ch);

if exist('Xdat')
   clear Xdat
end

%get XData of Children
for i = 1:si_ax(1)
 Xdat{i} = get(ax_ch(i),'XData');
 Ydat{i} = get(ax_ch(i),'YData');
end

%if size(Xdat{3},2)==2
%   error('to many lines')
%end
tel=0;
%for all lines 1,2 are window lines 
for i = 3 : si_ax(1)
   tel=tel+1;
   %find xdata of detection lines
   L2 = Xdat{1}(1);
   L1 = Xdat{2}(1);

   %if L1 is bigger: L1=l2
   if L1 > L2
      anders = L2;
      L2 = L1;
      L1 = anders;
   end

   %Find out samplefreq and set L1 and L2 to datapoints 
   X_max = max(Xdat{i});
   X_min = min(Xdat{i});
   xreal = [X_max - X_min] / [size(Xdat{i},2) - 1];

   L1 = round([[L1 - X_min] / xreal] +1);
   L2 = round([[L2 - X_min] / xreal] +1);

   %find minima and maxima
   [Ymin, Xmin] = min(Ydat{i}(L1:L2));
   [Ymax, Xmax] = max(Ydat{i}(L1:L2));
   Xmin = [[Xmin - 1] * xreal + [L1 - 1] * xreal] + X_min;
   Xmax = [[Xmax - 1] * xreal + [L1 - 1] * xreal] + X_min;

   Maximum = [Xmax,Ymax];
   Minimum = [Xmin,Ymin];

   %get mean
   Ymean = mean(Ydat{i}(L1:L2));
   Ystd = std(Ydat{i}(L1:L2));

   if Xdat{1}(1) < Xdat{2}(1)
      Maxx=Xdat{2}(1); Minx=Xdat{1}(1);
   else
      Maxx=Xdat{1}(1); Minx=Xdat{2}(1);
   end

   linetype = get(ax_ch(i),'Linestyle');
   lcol= get(ax_ch(i),'Color');
   lwidth=get(ax_ch(i),'LineWidth');

   Lint=sprintf('Linetype : %s',linetype);
   
   disp('Values: ')

   %Display everything
   if si_ax(1) > 3
      tel=tel+1;   
      disp([Lint,'     Color(RGB): ',num2str(lcol),'     Width: ',num2str(lwidth)])
   end

   Win=sprintf('Window  = %8.3f, %8.3f',Minx,Maxx);
   Duur=sprintf('Duration= %8.3f',[Maxx-Minx]);
   Mean=sprintf('Mean    = %8.3f, Sd:%8.3f',Ymean,Ystd);
   Std=sprintf('Std     =  %8.3f',Ystd);
   Max=sprintf('Maximum = %8.3f, %8.3f',Xmax,Ymax);  
   Min=sprintf('Minimum = %8.3f, %8.3f',Xmin,Ymin);
   
   disp(Win)
   disp(Duur)
   disp(Mean)
   disp(Max)
   disp(Min)
   disp(' ');

   if si_ax(1) > 3
      Values{tel}{1}=[linetype,num2str(lcol),num2str(lwidth)];
      Values{tel}{2}=[Minx,Maxx;Ymean,Ystd];
      Values{tel}{3}=[Xmax,Ymax;Xmin,Ymin];   
   else
      Values=[Maximum;Minimum;Ymean,Ystd;Xmax,Ymax;Xmin,Ymin];
   end   
end
