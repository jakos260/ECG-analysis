function [subject, bsmdir, layfile] = subjectdir(sub, bp, issinus)

dirbeats_nice   = '/Biosemi/export/beats/';
subject         = ['NICE00' num2str(sub)];
bsmdir          = [bp subject dirbeats_nice];
if sub == 1, layfile = 'ams65_NICE001.mla'; else layfile = 'ams65.mla'; end

%% PIG DATA FROM PETER O.
% dirbeatspig     = '/Biosemi/export/AVG/beats/';
% dirbeatsbucket  = '/Biosemi/export/beats/';
% 
% if sub{1} == 7
%     subject     = 'Pig07';
%     bsmdir      = [bp 'Data/' subject dirbeatspig];
%     layfile     = 'pig_adam.mla';
% elseif sub{1} == 8
%     subject     ='Pig08';
%     bsmdir      = [bp 'Data/' subject dirbeatspig];
%     layfile     = 'pig_adam.mla';
% elseif sub{1} == 9
%     subject     ='Pig09';
%     bsmdir      = [bp 'Data/' subject dirbeatspig];
%     layfile     = 'pig_adam.mla';
% elseif sub{1} == 10
%     subject     ='Pig10';
%     bsmdir      = [bp 'Data/' subject dirbeatspig];
%     layfile     = 'pig_adam.mla';
% elseif strcmp(sub{1},'b01')
%     subject     ='Bucket01';
%     bsmdir      = [bp 'Data/' subject dirbeatsbucket];
%     layfile     = 'bucket.lay';
% elseif strcmp(sub{1},'b01s')
%     subject     ='Bucketsmall01';
%     bsmdir      = [bp 'Data/Bucket01/Biosemi/export/beats/'];
%     layfile     = 'bucket.lay';
% elseif strcmp(sub{1},'b03')
%     subject     = 'Bucket03';
%     bsmdir      = [bp 'Data/' subject dirbeatsbucket];
%     layfile     = 'bucket.lay';
% else
%     subject     = '';
%     bsmdir      = '';
% end
% 
% % if sinusrythm, change bsmdir:
% if issinus, bsmdir = regexprep(bsmdir,'/export/','/export/SR/'); end
% 
% %'pig64.mla'; 'Bucket01Removed.lay';