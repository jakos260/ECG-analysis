function plotnimBSM(PHI)

clf
maxphi= max(max(max(max(PHI),0)),abs(min(min((PHI)))));
LAY       =loadmat('nim64.mla');
dd=0;
for i=2:length(LAY)
% 	subplot(LAY(1,2),	LAY(1,1),(floor(LAY(i,2))-1)*LAY(1,1)+floor(LAY(i,1)));
	axes('Position',[(LAY(i,1)/LAY(1,1)), 1-(LAY(i,2)/LAY(1,2)),0.9/LAY(1,1),1/LAY(1,2)]);
	plot( [1:length(PHI)],PHI(LAY(i,3),:),'b')
 	text(0.8*length(PHI(LAY(i,3),:)),0.8*maxphi,num2str(LAY(i,3)),'Fontname','verdana','FontSize',6);
	axis([0,length(PHI),-maxphi,maxphi]);
	axis off
	
% 	if LAY(i,3)==63
% 		hl=legend('measured','simulated','location', 'SE');
% 		legend('boxoff')
% 		set(hl,'Fontname','verdana','FontSize',6);
% 	end
end
% subplot(LAY(1,2), LAY(1,1),LAY(1,2)*LAY(1,1))
axes('Position',[1-1.5*0.9/LAY(1,1), 0.5*1/LAY(1,2),0.9/LAY(1,1),1/LAY(1,2)]);
line([ length(PHI)-100 length(PHI)], [-maxphi*0.9 -maxphi*0.9],'Linewidth',2,'Color','k');
line([ length(PHI) length(PHI)], [-maxphi*0.9 -maxphi*0.9+0.05],'Linewidth',2,'Color','k');
text(length(PHI)*0.4,-maxphi*0.75,'100 ms','Fontname','verdana','FontSize',7);
text(length(PHI)+2,-maxphi*0.7,'50 \muV','Fontname','verdana','FontSize',7);
axis([0,length(PHI),-maxphi,maxphi]);axis off

