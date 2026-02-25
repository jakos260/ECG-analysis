% load_mixed_mat
% reads asci file containing numbers and characters 
% script needs to be tuned to specific input formats
% 20070206; Mathieu Lemay


% name_I =  ['PP08.ELS'];
% name_II = ['QQ08.ELS'];

f=fopen(name_I, 'r');

[N,nr]=fscanf(f,'%d',1);

for I = 1:N(1),
    for J = 1:4,
        tmp = str2num(fscanf(f,'%s',1));
        if J < 4,
            Matrix(I,J) = tmp;
        end
    end
end

f=fopen(name_II,'wt');
index = [N(1) 3];
fprintf(f, '%4.0f %6.0f\n',index)';
fprintf(f, '%4.0f %6.4f %6.4f\n', Matrix');
fclose(f)