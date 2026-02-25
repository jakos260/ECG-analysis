

d = dir('*1061-1068*.selecg');
LAY= loadmat('nim65.mla');
close all
figure(1);clf;
for i=1:length(d)
    basename = d(i).name(1:end-7);
    ECG = loadmat([ basename '.selecg']);
    specs =loadmat([ basename '.spe']);
    ECG = baselinecor(ECG(:,specs(2):specs(5)));
    rrms = rms(ECG);
    figure(1)
    plot(rrms);
    hold on;
    eval(['ECG' num2str(i) '= ECG;']);
    figure(3)
    sigplot_p(ECG,'lay',LAY);
    hold on
end
