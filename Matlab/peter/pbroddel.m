
global DATA
ECG=readECG(fn ,'erasmus',0, 300);
selectBSMBeats(ECG,'lay',LAY,'filename',['.\beats\' fn]);
leadv16(DATA.BSMOUT(:,DATA.BEATS(min(find(DATA.SELBEATS==1))):DATA.BEATS(max(find(DATA.SELBEATS==1))+1)),'leadsys','9lds','max',[-2.5 2.5],'paperspeed',25);
saveas(gcf,[fn '.png']);