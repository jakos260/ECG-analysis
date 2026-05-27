% make_strip.m
% function ITRI=make_strip(VER,lista,listb)
% perform triangulation of a strip bounded on one side by vertices
% VER(lista,:), orderly listed along one of the long side
% and VER(listb,:) along the other side; 
% lista(1) should be directly opposite listb(1)
% lista should be disjunct with listb
% see also make_peel
% A. van Oosterom 2010_12_06



function ITRI=make_strip(VER,lista,listb)
ntri=0;
ITRI=[];

if size(lista,1)<size(lista,2), lista=lista'; end
if size(listb,1)<size(listb,2), listb=listb'; end
na=size(lista,1);
nb=size(listb,1);

if na<2|nb<2,
    beep
    'incorrect input data; na<2|nb<2'
    pause
end
    
   
ntri_max=na+nb-2;
ibest=1; jbest=1;
ntri=1;
 
while ntri<=ntri_max,
    d1=inf;
    if ibest+1<=na,
       d1=norm3d(VER(lista(ibest+1),:)-VER(listb(jbest),:) );
    end
    d2=inf;
     if jbest+1<=nb,
         d2=norm3d(VER(lista(ibest),:)  -VER(listb(jbest+1),:));
     end
    if d1<=d2,
          ITRI(ntri,:)=[lista(ibest) lista(ibest+1) listb(jbest)];
          if ibest<na,
              ibest=ibest+1;
          end
     else,
          ITRI(ntri,:)=[lista(ibest) listb(jbest+1) listb(jbest)];
          if jbest<nb, 
             jbest=jbest+1;
          end
    end
    ntri=ntri+1;  
end

 



