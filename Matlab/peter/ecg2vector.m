function ecg2vector(ECGin,varargin)


KI_9=[ -0.5267    0.1634   -0.2867   -0.1300    0.0500   -0.0100    0.1400    0.0600    0.5400;   % Vx
	  -0.8633   -0.0733    0.9266    0.0600   -0.0200   -0.0500    0.0600   -0.1700    0.1300;   % Vy
	   0.3300    0.3200   -0.0200   -0.4300   -0.0600   -0.1400   -0.2000   -0.1100    0.3100;];  % VzR,VL VF V1-V6]
cols=['b  ';'r  ';'k  ';'g  ';'m  ';'y  ';'r--';'k--';'g--';'c--';'m--';'y--'];
styles=['k- ';'r- ';'k--' ; 'r--' ;'k: ' ; 'k-.'];        
figure(105); clf
if ~iscell(ECGin)
    ECG{1} = ECGin;
else
    ECG = ECGin;
end

maxV =0;
for i=1:length(ECG)
    ecg = ECG{i};

    
    if size(ecg,1)==12
        ecg=ecg(4:end,:);
    end

    V= KI_9 * ecg;
    maxV=max([maxV; abs(V(:))])
end

for i=1:length(ECG)
    ecg = ECG{i};

    
    if size(ecg,1)==12
        ecg=ecg(4:end,:);
    end

    V= KI_9 * ecg;
    V=V';

    subplot(1,3,1)
    plot(V(:,1),V(:,2),cols(i,:))
    hold on
    axis([-maxV maxV -maxV maxV])
    
    subplot(1,3,2)
    plot(V(:,1),V(:,3),cols(i,:))
    hold on
    axis([-maxV maxV -maxV maxV])

    
    subplot(1,3,3)
    plot(V(:,2),V(:,3),cols(i,:))
    hold on
    axis([-maxV maxV -maxV maxV])

end

subplot(1,3,1)
title('horizontal plane')
subplot(1,3,2)
title('frontal plane')
subplot(1,3,3)
title('sagital plane')
