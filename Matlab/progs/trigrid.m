% file trigrid.m
% ITRI: triangle indices
% VER: vertices
% pp: plotparameters as in "old" triplo: er, alpha, p, izin

function trigrid(VER,ITRI,pp)
ntri=size(ITRI);
nver=size(VER);
er=pp(1)*pi;
alpha=pp(2)*pi;
p=pp(3);
izin=pp(4);
cer=cos(er);
ser=sin(er);
cal=cos(alpha);
sal=sin(alpha);
for i=1:ntri(1),
    for j=1:3,
        k=ITRI(i,j);
        for l=1:3,
            h(l)=VER(k,l);
        end
        hpk1=h(1)*cer-h(2)*ser;
        hpk2=h(1)*ser+h(2)*cer;
        hpk3=h(3);
        yt(j)=hpk3+p*hpk2*sal;
        xt(j)=hpk1+p*hpk2*cal;
     end
     zin=(xt(2)-xt(1))*(yt(3)-yt(1))-(xt(3)-xt(1))*(yt(2)-yt(1));
     if zin ~= 0,
        xt(4)=xt(1);
        yt(4)=yt(1);
        if zin > 0,
           line([xt(1) xt(2) xt(3) xt(4)],[yt(1) yt(2) yt(3) yt(4)]);
        elseif zin < 0,
           if izin == 2,
              set(line([xt(1) xt(2) xt(3) xt(4)],[yt(1) yt(2) yt(3) yt(4)]),'color',[.5 .5 0]);
           end
        end
     end
end


