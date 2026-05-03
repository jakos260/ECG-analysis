function [blocks estr Aout]=index2blocks(Ain,preproc)
% Convert 1D-array with increasing numbers (presumably indices to a 1D array)
% into blocks of consecutive indices/numbers.
% if Y=X(Aout), then blocks.org* are indices to X and block.dest* to Y.
% If preproc is true (default false) indices order is sorted. Duplex indices are deleted. Changing Ain into Aout.
% estr is a string that can be used by with eval to reconstruct the
% indices: Aout=eval(['[' estr ']' ]
if ~exist('preproc','var')
    preproc=false;
end

A=squeeze(Ain);
if size(A,1)>1
    error('Only 1D index arrays allowed');
end
if preproc
    A=sort(A);
    izero=find(diff(A)==0);
    if ~isempty(izero)
        warning('Duplex indices detected');
        A(izero+1)=[];
    end
end
Aout=A;

deststartblock=find(diff(A)>1 | diff(A)<=0);
destendblock=[deststartblock length(A)]; % last is allways end of block
deststartblock=[1 (deststartblock+1)]; % first sample always start of block, diff index is shifted by 1

orgstartblock=A(deststartblock);
orgendblock=A(destendblock);

blocks.deststart=deststartblock;
blocks.destend=destendblock;
blocks.orgstart=orgstartblock;
blocks.orgend=orgendblock;
blocks.length=orgendblock-orgstartblock+1;
if length(orgstartblock)==1
    blocks.gap=0;
else
    blocks.gap=orgstartblock(2:end)-orgendblock(1:end-1)-1;
    blocks.gap=[blocks.gap 0];
end

estr='';
for i=1:length(orgstartblock)
    if orgstartblock(i)==orgendblock(i) % descending series handled as individual numbers
        estr=[estr sprintf(' %d',orgstartblock(i))];
    else
        estr=[estr sprintf(' %d:%d',orgstartblock(i),orgendblock(i))];
    end
end



