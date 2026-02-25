% show_elecs
% A. van Oosterom; 20130324
% calling script should specify ELPOS and, if needed, ELLABS'

if ~exist('ELPOS'),
    'error: calling script should specify ELPOS and, if needed, ELLABS'
    pause
end

nshow=size(ELPOS,1);
fact=1.01;

for j=1:nshow,
    
    markpos=[ELPOS(j,1) ELPOS(j,2) ELPOS(j,3)];
    line(markpos(1),markpos(2),markpos(3),'LineStyle','none','Marker','o',...
        'MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor',[1 0 0]);
    
%     if ~exist('elnums')
        if ~isempty(elnums),
            lablepos=fact*markpos;
            lable=num2str(j);
            marklabs= text(lablepos(1),lablepos(2),lablepos(3),lable);
            set(marklabs,'FontSize',12);
        end
%     end
    
    if exist('ELLABS')
        if ~isempty(ELLABS),
            lablepos=fact*markpos;
            lable=ELLABS(j,:);
            marklabs= text(lablepos(1),lablepos(2),lablepos(3),lable);
            set(marklabs,'FontSize',12);
        end
    end
    
end


