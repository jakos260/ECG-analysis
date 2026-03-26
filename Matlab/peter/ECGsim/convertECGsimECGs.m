function ECG = convertECGsimECGs(fn)

ECGsim= loadmat(fn);
ECG = ECGsim([ 1:3 10:12 4:9],:);
savemat([fn '.ecg'],ECG);


%%
function M = loadmat(name)

f=fopen(name);
if (f==-1)
    fprintf('\nCannot open %s\n\n', name);
    M=0;
else
    [N,nr]=fscanf(f,'%d',2);
    if (nr~=2)
        fclose(f);
        f=fopen(name);
        N=fread(f,2,'long');
        M=fread(f,[N(2),N(1)],'float');       
    else
        M=fscanf(f,'%f',[N(2) N(1)]);
    end
    fclose(f);
    M=M';
end

%%
function savemat(name, M)

M=M';
f=fopen(name, 'wb');
n=size(M);
fwrite(f, n(2), 'long');
fwrite(f, n(1), 'long');
fwrite(f, M, 'float');
fclose(f);
