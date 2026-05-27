% make_peel.m
% function [ITRI,lista, listb]=make_peel(VER,lista,listb,mode,fig)
% perform triangulation of a peel bounded on one side by vertices VER(lista,:)
% specified by indices lista, and on the other side by VER(listb,:), both
% lists (node labels; column vectors) must be listed orderly along the peel, while having
% the same sense of rotation. list values should be unique, apart from start and finish!
% if mode==1, lista and listb are used as specified,
%   else, (default) triangulation will start
%         at the nodes in lista and listb that are closest.
% if nargin>4, peel is displayed in figure; fig
% blue: lista
% output lista and listb are the ones used
% A. van Oosterom 2011_04_12

function [ITRI,lista,listb]=make_peel(VER,lista,listb,mode,fig)



if size(lista,1)<size(lista,2), lista=lista'; end
if size(listb,1)<size(listb,2), listb=listb'; end
na=size(lista,1);
nb=size(listb,1);

if nargin<4, mode=0; end


% temporarily remove duplicate string endings
if na>1 & norm3d(VER(lista(1),:)-VER(lista(na),:))<=eps, lista(na)=[]; na=na-1; end
if nb>1 & norm3d(VER(listb(1),:)-VER(listb(nb),:))<=eps, listb(nb)=[]; nb=nb-1; end

if size(unique([lista;listb]),1)~=na+nb, 'error: duplicate list entries',pause,  return, end

na1=size(lista,1);
nb1=size(listb,1);

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
    
    % check consistency of orientation of the loop; reverse listb if needed
    a1=lista(1); a2=lista(2);
    b1=listb(1); b2=listb(2);
    cross1=cross(VER(a1,:)-VER(b1,:),VER(b2,:)-VER(b1,:));
    cross2=cross(VER(a1,:)-VER(b2,:),VER(a2,:)-VER(b2,:));
    cross1*cross2'
    
    if cross1*cross2'<=0,
        listb=reverse(listb);
        
    end
    
    b1=listb(1); b2=listb(2);
    cross1=cross(VER(a1,:)-VER(b1,:),VER(b2,:)-VER(b1,:));
    cross2=cross(VER(a1,:)-VER(b2,:),VER(a2,:)-VER(b2,:));
    cross1*cross2'
end

% close the peel edges
if na>1, lista=[lista;lista(1)]; na=na+1; end
if nb>1, listb=[listb;listb(1)]; nb=nb+1; end

ntri_max=na+nb-2;


i=1; j=1;
ntri=1;

if nargin>4
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

listb
lista


LINKSa=[]; % link all nodes of lista to the nearest one of listb
for i=1:na-1,
    dist=norm3d(VER(listb,:)-ones(nb,1)*VER(lista(i),:));
    [mi,j]=min(dist);
    plot3(VER([lista(i) listb(j)],1),VER([lista(i) listb(j)],2),VER([lista(i) listb(j)],3),'b')
    LINKSa(i,:)=[lista(i) listb(j)];
end

LINKSb=[];  % link all nodes of listb to the nearbest one of lista
for j=1:nb-1,
    dist=norm3d(VER(lista,:)-ones(na,1)*VER(listb(j),:));
    [mi,i]=min(dist);
    plot3(VER([lista(i) listb(j)],1),VER([lista(i) listb(j)],2),VER([lista(i) listb(j)],3),'k')
    LINKSb=[LINKSb;[listb(j) lista(i)]];
end

LINKSa;
size(LINKSa);

LINKSb;
size(LINKSb);

ALL=[LINKSa; LINKSb]

pause

ALL=sort(ALL,1)
end
% 
% % remove duplicates
% LINKS=[];
% while isempty(ALL)==0,
%     LINKS=[LINKS;ALL(1,:)];
%     nall=size(ALL,1);
%     k=find(sum(ALL==ones(nall,1)*ALL(1,:),2)==2 );
%     ALL(k,:)=[];
% end
% 
% %LINKS
% nlinks=size(LINKS,1);
% 
% for  i=1:nlinks,
%     nods=[lista(LINKS(i,1)) listb(LINKS(i,2))];
%     if nargin>4, plot3(VER(nods,1),VER(nods,2),VER(nods,3),'g'), end
% end
% hold on
% 
% TEST=[LINKS sum(LINKS,2)-1];
% ntest=size(TEST,1);
% 
% %find the missing links
% k=find(TEST(2:ntest,3)-TEST(1:ntest-1,3)>1);
% 
% ADD=zeros(size(k,1),2);
% if isempty(k)==0;
%     for kk=k,
%         d1=norm3d(VER(TEST(kk+1,1),:)-VER(TEST(kk,2),:));
%         d2=norm3d(VER(TEST(kk,1),:)  -VER(TEST(kk+1,2),:));
%         if d1<=d2,
%             ADD=[TEST(kk+1,1) TEST(kk,2)];
%         else
%             ADD=[TEST(kk,1) TEST(kk+1,2)];
%         end
%     end
%     
%     for  i=1:size(ADD,1),
%         nods=[lista(ADD(i,1)) listb(ADD(i,2))];
%         if nargin>4, plot3(VER(nods,1),VER(nods,2),VER(nods,3),'r'), end
%     end
%     
%     LINKS=[LINKS;ADD;];
%     LINKS=sort(LINKS,1);
%     nlinks=size(LINKS,1); 
% end
% 
% pause
% 
% 
% 
% % generate triangle indices
% ntris=nlinks-1;
% ITRI=zeros(ntris,3);
% 
% for itri=1:ntris,
%     quad=[LINKS(itri,1) LINKS(itri+1,1) LINKS(itri,2) LINKS(itri+1,2)];
%     if quad(1)==quad(2),
%         ITRI(itri,:)=[lista(quad(1)) listb(quad(4)) listb(quad(3))];
%     else
%         ITRI(itri,:)=[lista(quad(1)) lista(quad(2)) listb(quad(3))];
%     end
%     
% end





