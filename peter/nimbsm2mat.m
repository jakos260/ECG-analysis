function nimbsm2mat

doexchange2829=0;
doexchange2122=0;

if doexchange2829
    disp('exchanging lead 28 and 29 ');
end
if doexchange2122
    disp('exchanging lead 21 and 22 ');
end

direc= dir('*.bsm');
[VER,ITRI]= loadtri('C:\Users\damp2\Documents\Data\geometries\ppd\ppd2\ppd2_thorax65_Nijmegen.tri');
T = intripol(VER,ITRI,1:65);

for i=1:length(direc)
    A=loadbsmnim_VRL(direc(i).name);
    if size(A,2) ==9999
%         A=[A; zeros(1,size(A,2))];  


        A=[A; -sum(A([1 2],:))];
        if  doexchange2829
            B= A(29,:);
            A(29,:) = A(28,:);
            A(28,:) = B;
        end
        if  doexchange2122
            B= A(21,:);
            A(21,:) = A(22,:);
            A(22,:) = B;
        end

%         savemat([direc(i).name(1:end-4) '.mat'], A);
        hf=figure(1);
        pos = get(hf,'Position');
%         set(hf,'position',[pos(1) pos(2) 1150 700]);
        A = A - ones(size(A,1),1)*mean(A([1 2 65],:));
        A = A - mean(A')'*ones(1,size(A,2));
        figure(1);
        leadv16(A(:,5000:7000),'leadsys','nim','max',2.5,'paperspeed',25)
        saveas(hf,[direc(i).name(1:end-4) '.png']);
%         if i==1
            figure(2);
            Aext= T*A;
            tt=find(abs(Aext(19,:)) == max(abs(Aext(19,:)) ));
            showpatch(VER,ITRI,Aext(:,tt(1)))
            saveas(figure(2),[direc(i).name(1:end-4) 'bsm.png']);
%         end
%         pause
    end
end
