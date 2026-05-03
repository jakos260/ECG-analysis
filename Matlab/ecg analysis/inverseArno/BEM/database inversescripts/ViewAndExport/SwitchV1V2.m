clear all
% used to correct V1 and V2 switch in original amsterdam65.mla

folder='/Users/peteroosterhoff/Documents/Werk/Brugada/Ajm009_5346982/BSPM/hagen/export';
d=dir(fullfile(folder,'*.mat'));
for i=1:length(d)
    S=load(fullfile(folder,d(i).name));
    if S.mlas(2,3)==18
        S.mlas([2,3],:)=S.mlas([3,2],:);
        save(fullfile(folder,d(i).name),'-struct','S');
    else
        display('V2,V1 switch not needed for %f',d(i).name);
    end
end