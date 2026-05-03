function [H,SIGS]= ReadBDF(filename)

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
channel=1;
hw=waitbar(0,'reading data');
sig = H.sensor.gain(channel)*ChannelReaderBDF(realname,H.channels,H.nSamples,H.nTrials,channel,H.sampleRate,endOfFile,sizeHeader)';
SIGS=zeros(H.channels,length(sig));
SIGS(1,:)=sig;

for channel=2:H.channels
	SIGS(channel,:)=H.sensor.gain(channel)*ChannelReaderBDF(realname,H.channels,H.nSamples,H.nTrials,channel,H.sampleRate,endOfFile,sizeHeader)';
% 	SIGS(i,:)=ReadBDFchannel(filename,i);
	waitbar(channel/H.channels,hw);drawnow
end
close(hw)

msigs=mean(SIGS);
for i=1:size(SIGS,1)
	SIGS(i,:)=SIGS(i,:)-msigs;
end
