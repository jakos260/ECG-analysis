function ECG = nim2015_2_12lead(BSM)
ECG=[];
% index=[65 66 67 13 25 33 41 48 53];
index=[65 66 67 14 26 34 43 49 54];

if size(BSM,1)==67
    ECG = BSM(index,:);  
    ECG = bsxfun(@minus,ECG,mean(ECG')');

    wct = mean(ECG(1:3,:));
    ECG = bsxfun(@minus,ECG,wct);
    ECG= [ ECG(2,:) - ECG(1,:);...
           ECG(3,:) - ECG(1,:);...
           ECG(3,:) - ECG(2,:);...
           ECG(1:3,:) *1.5;...
           ECG(4:end,:)];
end