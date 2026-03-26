copdir='reconstruct678\';
t0=678;

a=dir;
for i=1:length(a)
	if ~a(i).isdir && isdicom(a(i).name)
		info=dicominfo(a(i).name);
		if isfield(info,'TriggerTime') && ~strcmp(info.SeriesDescription,'CSPAMM') &&...
			info.TriggerTime>t0-5 && info.TriggerTime<t0+5
			copyfile(a(i).name,[copdir a(i).name]);
			disp(num2str(i))
		end
	end
end
