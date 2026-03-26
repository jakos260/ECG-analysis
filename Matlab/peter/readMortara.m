function ECG = readMortara(varargin)

fn = varargin{1};
tbegin = 1;
if nargin >= 2
    tbegin = max(1,round(varargin{2}*1000.0));
end
tend = inf;
if nargin == 3
    tend = round(varargin{3}*1000);
end

fid=fopen(fn,'r');
ECG=fread(fid,'int16')*15/16000;
fclose(fid);
ECG=reshape(ECG,8,length(ECG)/8);
tend = min(tend,size(ECG,2));
ECG=ECG(:,tbegin:tend);

Vr = -1/3 * ( ECG(1,:) + ECG(2,:) );
Vl = ECG(1,:) + Vr; %1/3 * ( ECG(2,:) - 2 * ECG(1,:));
Vf = ECG(2,:) + Vr; 

% II = Vl - Vr
% Vf = Vr - ECG(2,:);  II = Vr-Vf
% III= ECG(1,:) - ECG(2,:);
% ECG = [ECG(1:2,:); III; Vr; Vl; Vf ;ECG(3:end,:); ];
ECG = [Vr; Vl; Vf ;ECG(3:end,:); ];
