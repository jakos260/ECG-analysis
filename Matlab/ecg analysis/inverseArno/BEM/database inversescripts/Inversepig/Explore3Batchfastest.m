clearvars
close all
% load('BatchFastestPathOutput20130924165409');
% load('BatchFastestPathOutput20130929114044'); % pace start
% load('BatchFastestPathOutput20130930111304'); % QRS start
% load('BatchFastestPathOutput20131001094235'); % Only initvelocity==0, qrsonset
load('BatchFastestPathOutput20131001105836');%Only initvelocity==0, pace
overlayopp=2; % if 1 plot opposing foci as overlay in different color, 2 plots ratio

for i=1:length(FileName)
    indendo(i)=~isempty(findstr(FileName{i},'Endo'));
end
indepi=~indendo;

if qrsstart
    startstr='Start at QRSonset';
else
    startstr='Start at stimulus';
end

% RDfull(ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,igeom)
varstr={'RDfull','RDhalf','RDWfull','RDWhalf','CORfull','CORhalf'};
s=size(RDfull);
for i=1:prod(s)
    [ianisinit,iinitialVelocity,iinitialActTime,ianis,iopposite,igeom]=...
        ind2sub(size(RDfull),i);
    %     descriptavg{i_CORRECT FOR IGEOM}=sprintf('AnisI %0.2f,VeloI %.01f,DurI %d.Anis:%0.2f. Opposite%d',...
    %         aniss(ianisinit),velocities(iinitialVelocity),inittimes(iinitialActTime), aniss(ianis),iopposite-1);
    %     siganisinit(i)=aniss(ianisinit);
    %     siginitialVelocity(i)=velocities(iinitialVelocity);
    %     siginitialActTime(i)=inittimes(iinitialActTime);
    %     siganis(i)=aniss(ianis);
    %     sigopposite(i)=iopposite-1;
    sigt(1,i)=aniss(ianisinit);
    sigt(2,i)=velocities(iinitialVelocity);
    sigt(3,i)=inittimes(iinitialActTime);
    sigt(4,i)=aniss(ianis);
    sigt(5,i)=iopposite-1;
    sigt(6,i)=igeom;
end
sig=sigt(1:5,sigt(6,:)==1); % only per geom.
sigmin=repmat(min(sig,[],2),1,size(sig,2));
sigmax=repmat(max(sig,[],2),1,size(sig,2));

for i=1:length(varstr);
    V=eval(varstr{i});
    %     figure('name',[varstr{i} ' ratio individual plots']);
    %     plot(reshape(V,[],size(V,6)));
    avgV=mean(V,6);
    avgendoV=mean(V(:,:,:,:,:,indendo),6);
    avgepiV=mean(V(:,:,:,:,:,indepi),6);
    sem=std(V,1,6)/sqrt(size(V,6));
    plotsig=((sig-sigmin)./(sigmax-sigmin)*0.2+repmat([0:0.2:0.8]'-1+min(avgV(:)),1,size(sig,2)));
    
    fh(1)=figure('Name',[varstr{i} ' ' startstr]);
    if overlayopp==1
        plot(avgV(1:end/2),'b');
        hold on
        plot(avgV(1:end/2)+sem(1:end/2),':b');
        plot(avgV(end/2+1:end),'r');
        plot(avgV(end/2+1:end)+sem(end/2+1:end),':r');
        legend(varstr{i},[varstr{i} '+SEM'],[varstr{i} '-Opp'],[varstr{i} '-Opp+SEM'],['anisinit[' num2str(aniss,'%0.2f/') ']'],['InitVelocity[' num2str(velocities,'%0.2f/') ']'],['initActivationTime['  num2str(inittimes,'%0.2f/') ']'],['anisotropoy[' num2str(aniss,'%0.2f/') ']'],'opposite[0/1]');

        
    elseif overlayopp==0
        plot(avgV(1:end/2),'k');
        hold on
        plot(avgV(1:end/2)+sem(1:end/2),':k');
        plot(plotsig(:,1:end/2)');
        legend(varstr{i},[varstr{i} '+SEM'],['anisinit[' num2str(aniss,'%0.2f/') ']'],['InitVelocity[' num2str(velocities,'%0.2f/') ']'],['initActivationTime['  num2str(inittimes,'%0.2f/') ']'],['anisotropoy[' num2str(aniss,'%0.2f/') ']'],'opposite[0/1]');

        
        
        fh(2)=figure('Name',[varstr{i} ' Epi ' startstr]);
        plot(avgepiV(1:end/2),'b');
        hold on
        plot(avgepiV(end/2+1:end),'r');
        plot(plotsig(:,1:end/2)');
        legend(varstr{i},[varstr{i} '_Opp'],['anisinit[' num2str(aniss,'%0.2f/') ']'],['InitVelocity[' num2str(velocities,'%0.2f/') ']'],['initActivationTime['  num2str(inittimes,'%0.2f/') ']'],['anisotropoy[' num2str(aniss,'%0.2f/') ']'],'opposite[0/1]');

        
        fh(3)=figure('Name',[varstr{i} ' Endo ' startstr]);
        plot(avgendoV(1:end/2),'b');
        hold on
        plot(avgepiV(end/2+1:end),'r');
        plot(plotsig(:,1:end/2)');
        legend(varstr{i},[varstr{i} '_Opp'],['anisinit[' num2str(aniss,'%0.2f/') ']'],['InitVelocity[' num2str(velocities,'%0.2f/') ']'],['initActivationTime['  num2str(inittimes,'%0.2f/') ']'],['anisotropoy[' num2str(aniss,'%0.2f/') ']'],'opposite[0/1]');
          
        fh(4)=figure('name',[varstr{i} '  individual plots']);
        plot(reshape(V,[],size(V,6)));
        hold on
        set(gca,'XGrid','on');
        set(gca,'XTick',1:54:length(V(:))/2);
        set(gca,'XMinorTick','off');
        plot(plotsig(:,1:end/2)');
        legend([FileName {['anisinit[' num2str(aniss,'%0.2f/') ']'],['InitVelocity[' num2str(velocities,'%0.2f/') ']'],['initActivationTime['  num2str(inittimes,'%0.2f/') ']'],['anisotropoy[' num2str(aniss,'%0.2f/') ']'],'opposite[0/1]'}]);


        
        set(0,'CurrentFigure',fh(1));
        
        
        
    elseif overlayopp ==2
        ratio=V(:,:,:,:,1,:)./V(:,:,:,:,2,:);
        fh(2)=figure('name',[varstr{i} ' ratio individual plots']);
        plot(reshape(ratio,[],size(V,6)));
        hold on
        set(gca,'XGrid','on');
        set(gca,'XTick',1:54:length(V(:))/2);
        set(gca,'XMinorTick','off');
        plot(plotsig(:,1:end/2)');
        legend([FileName {['anisinit[' num2str(aniss,'%0.2f/') ']'],['InitVelocity[' num2str(velocities,'%0.2f/') ']'],['initActivationTime['  num2str(inittimes,'%0.2f/') ']'],['anisotropoy[' num2str(aniss,'%0.2f/') ']'],'opposite[0/1]'}]);
             
        
        %         set(0,'CurrentFigure',fh(1));
        fh(3)=figure('name',[varstr{i} ' Difference individual plots']);
        plot(reshape(V(:,:,:,:,1,:)-V(:,:,:,:,2,:),[],size(V,6)));
        hold on
        set(gca,'XGrid','on');
        set(gca,'XTick',1:54:length(V(:))/2);
        set(gca,'XMinorTick','off');
        plot(plotsig(:,1:end/2)');
        legend([FileName {['anisinit[' num2str(aniss,'%0.2f/') ']'],['InitVelocity[' num2str(velocities,'%0.2f/') ']'],['initActivationTime['  num2str(inittimes,'%0.2f/') ']'],['anisotropoy[' num2str(aniss,'%0.2f/') ']'],'opposite[0/1]'}]);

        
        set(0,'CurrentFigure',fh(1));
        set(fh(1),'Name',[varstr{i} ' ratio ' startstr]);
        avgratio=mean(ratio,6);
        plot(avgratio(:),'k');
        hold on
        semratio=std(ratio,1,6)/sqrt(size(V,6));
        plot(avgratio(:)+semratio(:),':k');
        plot(plotsig(:,1:end/2)');
        legend([varstr{i} '(focus/Opp)'],[varstr{i} '(focus/Opp)+SEM'],['anisinit[' num2str(aniss,'%0.2f/') ']'],['InitVelocity[' num2str(velocities,'%0.2f/') ']'],['initActivationTime['  num2str(inittimes,'%0.2f/') ']'],['anisotropoy[' num2str(aniss,'%0.2f/') ']'],'opposite[0/1]');

        
    end
    %     errorbar(avgV(:),1.96*sem(:),'k');
    set(gca,'XGrid','on');
    set(gca,'XTick',1:54:length(V(:))/2);
    set(gca,'XMinorTick','off');
    
    
    %     plotscaled(siganisinit);
    %     plotscaled(siginitialVelocity);
    %     plot(siginitialActTime);
    %     plot(siganis);
    %     sigopposite;
   
end

