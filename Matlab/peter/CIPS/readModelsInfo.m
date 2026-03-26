dirname ='C:\Users\peter\Documents\Data\measurements\UCLA\prosCIPS';
[INFO, info] = read_info(dirname);

for i=1 : length(info)
%     disp(info{i}.filename);
    name = info{i}.filename(length(dirname)+2:end);
    disp(name)
end


xlswrite(fullfile(dirname,'modelinfo.xlsx'),INFO,'info','B2')

figure(1)
clf
hold on
plot(INFO(:,7),INFO(:,2)-mean(INFO(:,2)),'x')
plot(INFO(:,7),INFO(:,3)-mean(INFO(:,3)),'o')
plot(INFO(:,7),INFO(:,4)-mean(INFO(:,4)),'d')

p = polyfit(x,y,1)




figure(2)
clf
hold on
plot(INFO(:,8),INFO(:,2)-mean(INFO(:,2)),'x')
plot(INFO(:,8),INFO(:,3)-mean(INFO(:,3)),'o')
plot(INFO(:,8),INFO(:,4)-mean(INFO(:,4)),'d')

figure(3)
clf
hold on
plot(INFO(:,9),INFO(:,2)-mean(INFO(:,2)),'x')
plot(INFO(:,9),INFO(:,3)-mean(INFO(:,3)),'o')
plot(INFO(:,9),INFO(:,4)-mean(INFO(:,4)),'d')

figure(4)
clf
hold on
plot(INFO(:,7).*INFO(:,9),INFO(:,2)-mean(INFO(:,2)),'x')
plot(INFO(:,7).*INFO(:,9),INFO(:,3)-mean(INFO(:,3)),'o')
plot(INFO(:,7).*INFO(:,9),INFO(:,4)-mean(INFO(:,4)),'d')