% make_strip.m
% function ITRI=make_strip(VER,lista,listb,fig,sub)
% perform triangulation of a strip bounded on one side by vertices
% VER(lista,:), orderly listed along one of the long side
% and VER(listb,:) along the other side;
% lista(1) should be directly opposite listb(1)
% to make a (close) peel ensure lista(end)=lista(1) and listb(end)=listn(1)
%  lista and listb should be disjunct
% see also make_peel
% if nargin>3,figure(fig) will monitor the process
% sub: (sub)loop label

% A. van Oosterom 2015_12_02

function ITRI=make_strip(VER,lista,listb,fig,sub)
% for a (closed) strip (=peel) check: lista(end)=lista(1) and listb(end)=listb(1)
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

na=size(lista,1);
nb=size(listb,1);

if na<2|nb<2,
    beep
    'incorrect input data; na<2|nb<2'
    %pause
end

%[na nb]
ntri_max=na+nb-2;
%pause

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
        for j=1:na,
            text(VER(lista(j),1),VER(lista(j),2),VER(lista(j),3),num2str(lista(j)))
        end
        for j=1:nb,
            text(VER(listb(j),1),VER(listb(j),2),VER(listb(j),3),num2str(listb(j)))
        end
    end
end
view(100, 40)
%pause

a1=1; b1=1; a2=2; b2=2;

while ntri<ntri_max-1,
    
    
    %'next will be based on step'
    status=[a1 a2 na b1 b2 nb];
    
    
    da=norm(VER(lista(a1),:) -VER(listb(b2),:) );
    if a1==na, da=0; end
    
    db=norm(VER(listb(b1),:) -VER(lista(a2),:));
    if b1==nb  , db=0; end
   
    if da<=db,
        addtri=[listb(b1) lista(a1) listb(b2)];
        %'is addtri based on [b1 a1 b2] :'
        [b1 a1 b2];
        %'next tri :'
        b1=b2;
        b2=min(b2+1,nb);
    else,
        addtri=[ listb(b1) lista(a1) lista(a2)];
        %'is addtri based on [b1 a1 a2]  : '
        [b1 a1 a2];
        %'next tri'
        a1=a2;
        a2=min(a2+1,na);
    end
    
    ITRI=[ITRI;addtri];
    ntri=ntri+1;
    %[ntri ntri_max];
    
    if nargin>3,
        if fig>0,
            quart=[addtri addtri(1)];
            plot3(VER(quart,1),VER(quart,2),VER(quart,3),'g')
            for j=1:3,
                text(VER(addtri(j),1),VER(addtri(j),2),VER(addtri(j),3),num2str(addtri(j)))
            end
            view(160,40)
        end
    end
    
    %beep
    %pause
end 

% 'final step was'
% status

if a1==na & a2==na,
    addtri=[lista(end) listb(end) listb(end-1)];
else
   addtri=[listb(end) lista(end-1) lista(end)]; 
end
  
ITRI=[ITRI;addtri];
ntri=ntri+1;

if nargin>3,
    if fig>0,
        quart=[addtri addtri(1)];
        plot3(VER(quart,1),VER(quart,2),VER(quart,3),'g')
        for j=1:3,
            text(VER(addtri(j),1),VER(addtri(j),2),VER(addtri(j),3),num2str(addtri(j)))
        end
        view(100,30)
    end
    %pause
end

% 'ntri ntri_max'
% [ntri ntri_max]

%pause

if size(ITRI,1)~=ntri_max, 
    beep
    'make_strip ERROR'
    pause
end


