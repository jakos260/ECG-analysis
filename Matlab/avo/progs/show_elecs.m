% show_elecs
% A. van Oosterom; 2017_01_26
% calling script should specify ELPOS and, if needed, ELLABS'

if ~exist('ELPOS','var'),
    'error: calling script should specify ELPOS and, if desired, ELLABS'
    pause
end

% elec_nodes should be the nodes n_pos to be marked 
% nodes_stdlds should specify the relevant 9 node loctions

n_pos=size(ELPOS,1); 
fact=1.03;

for j=1:n_pos,
    
     elpos=[ELPOS(j,1) ELPOS(j,2) ELPOS(j,3)];
     
     if ismember(elec_nodes(j),nodes_stdlds),
        line(elpos(1),elpos(2),elpos(3),'LineStyle','none','Marker','o',...
        'MarkerSize',8,'MarkerFaceColor','w','MarkerEdgeColor',[1 0 0]);
     else
       line(elpos(1),elpos(2),elpos(3),'LineStyle','none','Marker','o',...
        'MarkerSize',6,'MarkerFaceColor','k','MarkerEdgeColor',[1 0 0]);  
    end
    
end













%     if exist('ELLABS','var'),    
%         
%         ellabs=ELLABS;
%         
%         if size(ellabs)>0 & size(ellabs)<n_elecs,
%             if ismember(j,ellabs),
%                 k=find(ellabs==j);
%                 lable=ELLABS(k,1:3);
%                 lablepos=fact*elpos;
%                 marklabs= text(lablepos(1),lablepos(2),lablepos(3),lable);
%                 set(marklabs,'FontSize',12);
%             end
%             
%         end
%         
%        if size(ellabs)==n_elecs,
%                 lable=num2str(j);
%                 lablepos=fact*elpos;
%                 marklabs= text(lablepos(1),lablepos(2),lablepos(3),lable);
%                 set(marklabs,'FontSize',12);
%         end 
        
    
 




