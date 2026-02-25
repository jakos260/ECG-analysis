% refine4.m
% generates refined version of any triangulation [VER,ITRI] by inserting 
% extra nodes halfway ALL edges
% function [VER,ITRI]=refine4(VER,ITRI)
% VER, ITRI : vertex and triangle indices
% surface need NOT to be closed

% A. van Oosterom; 20090727

function [VER,ITRI]=refine4(VER,ITRI)
        ntri=size(ITRI,1);
        nver=size(VER ,1);
        TRI=[];
        for i=1:ntri,
            for j=1:3,
                jp=icyc(j+1,3);
                i1=ITRI(i,j); i2=ITRI(i,jp);
                h(1:3)=(VER(i1,:)+VER(i2,:))/2;
                %h=h/norm(h);
                n=[];
                n=find(sum(VER==ones(nver,1)*h,2)==3);
                if isempty(n)==1,
                   VER=[VER;h];
                   nver=nver+1;
                   index(j)=nver;
                else,
                   index(j)=n;
                end
            end
            TRI=[TRI; 
            ITRI(i,1)  index(1)   index(3);
            index(1)   ITRI(i,2)  index(2);
            index(3)   index(1)   index(2);
            index(3)   index(2)   ITRI(i,3);];
       end
       ITRI=TRI;
  