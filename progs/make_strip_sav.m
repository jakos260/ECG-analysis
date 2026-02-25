% make_strip.m
% function ITRI=make_strip(VER,lista,listb,fig,sub)
% perform triangulation of a strip bounded on one side by vertices
% VER(lista,:), orderly listed along one of the long side
% and VER(listb,:) along the other side;
% lista(1) should be directly opposite listb(1)
% to make a (close) peel set lista(end)=lista(1) and listb(end)=listn(1)
%  lista and listb should be disjunct
% see also make_peel
% if nargin>3,figure(fig) will monitor the process
% sub: (sub)loop label

% A. van Oosterom 2015_11_27

function ITRI=make_strip(VER,lista,listb,fig,sub)
% for a (closed) strip (=peel) check: lista(end)=lista(1) and listb(end)=listb(1)

lista
listb
pause

if nargin>3,
    if fig>0,
        figure(fig)
        if sub==1,
            clf,
        end
        if nargin>4,
            if sub==1,
                ccol='b';
            else
                ccol='m';
            end
        end
        plot3(VER(lista,1),VER(lista,2),VER(lista,3),['-+' ccol])
        hold on
        plot3(VER(listb,1),VER(listb,2),VER(listb,3),['m-+' ])
    end
end

ntri=0;
ITRI=[];

if size(lista,1)<size(lista,2), lista=lista'; end
if size(listb,1)<size(listb,2), listb=listb'; end

na=size(lista,1)
nb=size(listb,1)

if na<2|nb<2,
    beep
    'incorrect input data; na<2|nb<2'
    pause
end

ntri_max=na+nb-2

ntri=0;

if nargin>3,
    if fig>0,
        plot3(VER(lista,1),VER(lista,2),VER(lista,3),'+-')
        hold on
        plot3(VER(lista(1),1),VER(lista(1),2),VER(lista(1),3),'ro')
        plot3(VER(lista(end-1),1),VER(lista(end-1),2),VER(lista(end-1),3),'ko')
        plot3(VER(listb,1),VER(listb,2),VER(listb,3),['+-' ccol])
        plot3(VER(listb(1),1),VER(listb(1),2),VER(listb(1),3),'ro')
        plot3(VER(listb(end-1),1),VER(listb(end-1),2),VER(listb(end-1),3),'ko')
    end
end

% pause

a1=1; b1=1; a2=2; b2=2;

while ntri<ntri_max,
    
    if ntri<ntri_max-1,
        
        
        da=norm3d(VER(lista(a1),:)-VER(listb(b2),:) );
        
        
%         if a2==na,
%             db=inf;
%         else
            db=norm3d(VER(listb(b1),:) -VER(lista(a2),:));
       % end
        
        if da<=db,
            addtri=[ listb(b1) lista(a1) listb(b2)];
            b1=b2;
            b2=min(b2+1,nb);
        else,
            addtri=[ listb(b1) lista(a1) lista(a2)];
            a1=a2;
            a2=min(a2+1,na);
        end
        
    else,
        
        % specfy final triangle
        if a1==a2,
            addtri=[lista(a1) listb(b2) listb(b1)]
        else
            addtri=[listb(b1) lista(a1) lista(a2)]
        end
    end
    
    if nargin>3,
        if fig>0,
            quart=[addtri addtri(1)]
            plot3(VER(quart,1),VER(quart,2),VER(quart,3),'g')
            for j=1:3,
                text(VER(addtri(j),1),VER(addtri(j),2),VER(addtri(j),3),num2str(addtri(j)))
            end
            view(0,30)
        end
    end
    ITRI=[ITRI;addtri]
    ntri=ntri+1
    [a1 a2 na; b1 b2 nb]
    ntri_max
    %pause
end
%pause





