function [H,SIGS]= readHagaBDF(filename,t0,t1)



if ~isempty(strfind(filename,'.bdf'))
	filename=filename(1:end-4);
end

H=ImportBDFHeader(filename);
realname = [filename, '.bdf'];
sizeHeader = 256 + 256 * H.channels;
fin = fopen(realname,'r');
fseek(fin,0,'eof');
endOfFile = ftell(fin);
fclose(fin);


Exg1 = 66;H.channels - 8;
Exg2 = 65;H.channels - 7;
Exg3 = 67;H.channels - 6;
channels= [ 1:48 60 61 62 49 : 59 Exg1 Exg2 Exg3];
chnim65=[10:64 3:9 1 2 65 ];
if H.channels < 65 
    error('not enough channels')
end


channel=1;
sig = H.sensor.gain(channels(channel)) * ChannelReaderBDF(realname,H.channels,H.nSamples,H.nTrials,channels(channel),H.sampleRate,endOfFile,sizeHeader)';
sig = sig(max(1,t0*H.sampleRate):t1*H.sampleRate);
sig = resample(sig,1,length(sig),length(sig)*1000/H.sampleRate);
SIGS=zeros(length(chnim65),length(sig));
ichan = 1;
SIGS(chnim65(ichan),:)=sig;

% sensors zijn  7 kanalen
% 64 kanalen waarvan kanaal 63 en 64 niet gebruikt
% van de EXECG de eerste 3 kanalen
hw=waitbar(0,'reading data');


for channel=2:size(SIGS,1)%H.channels
	sig = H.sensor.gain(channels(channel)) * ChannelReaderBDF(realname,H.channels,H.nSamples,H.nTrials,channels(channel),H.sampleRate,endOfFile,sizeHeader)';
    sig = sig(max(1,t0*H.sampleRate):t1*H.sampleRate);
    SIGS(chnim65(channel),:) = resample(sig,1,length(sig),length(sig)*1000/H.sampleRate);
	waitbar(channel/size(SIGS,1),hw);drawnow
end
close(hw)

SIGS = SIGS/1000;



% msigs=mean(SIGS);
% for i=1:size(SIGS,1)
% 	SIGS(i,:)=SIGS(i,:)-msigs;
% end
