function slote2mat(name);
% SLOTEMAT	Convert an ATTL-file to an MFBF matrix file.
%
%		Usage: slote2mat('file');
%
%		SLOTE2MAT('file') converts the ATTL-data file 'file'
%		to the MFBF matrix file 'file.mat'.
%		The sample frequency is also stored in that file (using
%		the standard MFBF tail structure).
%
%		Thom Oostendorp, MF&BF University of Nijmegen, the Netherlands

[rec, info, sampint, rejected] = rload(name, 1);
outname=[name '.mat'];
extra=char('Sample frequency', num2str(1000/sampint));
savemat(outname, rec', extra);
