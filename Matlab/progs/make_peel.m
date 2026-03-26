% make_peel.m
% function [ITRI,lista, listb]=make_peel(VER,lista,listb,mode,fig)
% perform triangulation of a peel bounded on one side by vertices VER(lista,:)
% specified by indices lista, and on the other side by VER(listb,:), both
% lists (node labels; column vectors) must be listed along the peel, while having
% the same sense of rotation. List values should be unique, apart from start and finish!
% if mode==1, lista and listb are used as specified,
%   else, (default) triangulation will start
%         at the nodes in lista and listb that are closest.
% if nargin>4, peel is displayed in figure: fig
% output lista and listb are the ones used, adapted if mode.ne.1
% A. van Oosterom 2013_02_25 ; now forcing orientation of listb (also in output)
% to be consistent with that of lista; now calling make_strip

function [ITRI,lista,listb]=make_peel(VER,lista,listb,mode,fig)
ITRI=[];

if size(lista,1)<size(lista,2), lista=lista'; end % forcing column vectoctors
if size(listb,1)<size(listb,2), listb=listb'; end

na=size(lista,1);
nb=size(listb,1);

% temporarily remove duplicate string endings
if na>1 & norm3d(VER(lista(1),:)-VER(lista(na),:))<=eps, lista(na)=[]; na=na-1; end
if nb>1 & norm3d(VER(listb(1),:)-VER(listb(nb),:))<=eps, listb(nb)=[]; nb=nb-1; end

if size(unique([lista;listb]),1)~=na+nb, 'error: duplicate list entries',pause,  return, end

if na>1 & nb>1;
    na1=size(lista,1);
    nb1=size(listb,1);
    
    if nargin<4, mode=0; end
    
    if mode==0
        % identify the shortest connection between both loops
        ibest=1; jbest=1;
        small=inf;
        for i=1:na1,
            for j=1:nb1,
                d=norm3d(VER(lista(i),:)-VER(listb(j),:) );
                if d<small, ibest=i; jbest=j; small=d; end
            end
        end
        % reshuffle both lists such that lista(1) will be lista(ibest) and listb(1) will be listb(jbest)
        if ibest>1, lista=lista([ibest:na1 1:ibest-1]);end
        if jbest>1, listb=listb([jbest:nb1 1:jbest-1]);end
        
    end
    
    ntri_max=na+nb;
    
    % close the peel edges
    lista=[lista;lista(1)]; na=na+1;
    listb=[listb;listb(1)]; nb=nb+1;

     
% % section below: not full proof   
% %     check consistency of orientation of the loop; reverse listb if needed
%     a1=lista(1); a2=lista(2);
%     b1=listb(1); b2=listb(2);
%     cross1=cross(VER(a1,:)-VER(b1,:),VER(b2,:)-VER(b1,:));
%     cross2=cross(VER(a1,:)-VER(b2,:),VER(a2,:)-VER(b2,:));
%     cross1*cross2';
%     if cross1*cross2'<=0,
%         listb=reverse(listb);
%     end

% string=[lista(1); listb ; reverse(lista)]
% curl=sum(VER(string(2:end),:)-VER(string(1:end-1),:))
% if curl*curl'> 2*eps, 
%    listb=reverse(listb);
%    
%    string=[lista(1); listb ; reverse(lista)]
%    curl=sum(VER(string(2:end),:)-VER(string(1:end-1),:))
   
% end
% 
% pause


     
    if nargin>4,
        figure(fig)
        clf
        plot3(VER(lista,1),VER(lista,2),VER(lista,3),'+-')
        for i=1:na-1,
            text(VER(lista(i),1),VER(lista(i),2),VER(lista(i),3),num2str(lista(i)))
        end
        grid on
        hold on
        plot3(VER(listb,1),VER(listb,2),VER(listb,3),'+-k')
        for i=1:nb-1,
            text(VER(listb(i),1),VER(listb(i),2),VER(listb(i),3),num2str(listb(i)))
        end
        plot3(VER(lista(1),1),VER(lista(1),2),VER(lista(1),3),'*b', 'linewidth',2)
        plot3(VER(listb(1),1),VER(listb(1),2),VER(listb(1),3),'*k', 'linewidth',2) 
    end
    
    
    LINKSa=[]; % link all nodes of lista to the nearest one of listb
    for i=1:max(na-1,1);
        dist=norm3d(VER(listb,:)-ones(nb,1)*VER(lista(i),:));
        [mi,j]=min(dist);
        LINKSa(i,:)=[lista(i) listb(j)];
    end
    
    LINKSa;
    
    LINKSb=[];  % link all nodes of listb to the nearest one of lista
    for j=1:max(nb-1,1),
        dist=norm3d(VER(lista,:)-ones(na,1)*VER(listb(j),:));
        [mi,i]=min(dist);
        LINKSb=[LINKSb;[listb(j) lista(i)]];
    end
    
    % identify common_links
    ALL=[LINKSa;LINKSb(:,[2 1])];
    LINKS=[];
    
    while ~isempty(ALL),
        k=find(ALL(1,1)==ALL(:,1)& ALL(1,2)==ALL(:,2));
        if size(k,1)>1,
            LINKS=[LINKS; ALL(1,:)];
            if nargin>4,
                plot3(VER(ALL(1,:),1),VER(ALL(1,:),2),VER(ALL(1,:),3),'m')
            end
        end
        ALL(k,:)=[];
    end
    
    LINKS=[LINKS;LINKS(1,:)];
    
    ia_beg=1; ib_beg=1;
    
    list_a=lista(1:end-1);
    list_b=listb(1:end-1);
    
    naa=size(list_a,1);
    nbb=size(list_b,1);
    
    
    i=1;
    while size(ITRI,1)<ntri_max, % triangulation of peel in between subsequent segments (LINKS)
        i=i+1;
        ia_end=find(list_a==LINKS(i,1));
        
        if ia_end>ia_beg,
            listaa=list_a(ia_beg:ia_end);
        else
            listaa=list_a([ia_beg:naa 1:ia_end]);
        end
        
        ib_end=find(list_b==LINKS(i,2));
        
        if ib_end>ib_beg,
            listbb=list_b(ib_beg:ib_end);
        else
            listbb=list_b([ib_beg:nbb 1:ib_end]);
        end
        
        ADDTRIS=make_strip(VER,listaa,listbb);
        ntr=size(ADDTRIS,1);
        if nargin>4,
            for j=1:ntr,
                trits=ADDTRIS(j,:);
                plot3(VER(trits,1),VER(trits,2),VER(trits,3),'r')
            end
        end
        ITRI=[ITRI;ADDTRIS];
        size(ITRI);
        ia_beg=ia_end;
        ib_beg=ib_end;
      
    end
    
else
    % disc like geometry
    
    if na==1,
        ntr=size(listb);
        listb=[listb; listb(1)];
        for i=1:ntr,
            ITRI(i,:)=[lista(1) listb(i+1) listb(i)];
        end
    else,
        ntr=size(lista);
        lista=[lista; lista(1)];
        for i=1:ntr,
            ITRI(i,:)=[lista(i) lista(i+1) listb(1)];
        end
    end
   
    'orientation cap triangles may need to be checked'
end


