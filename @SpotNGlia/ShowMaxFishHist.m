function ShowMaxFishHist(obj)

load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo')

b = zeros(numel(RegistrationInfo),3);
for k2=1:numel(RegistrationInfo)
    b(k2, 1:3) = RegistrationInfo{k2}(find(strcmp({RegistrationInfo{k2}.name}, 'MaxFishColor'))).value;
end
bm = mean(b(:));

figure; 
bar(mean(b, 2));
%bar(b(:,3));
set(gca,'YLim',[0,265]);
hold on
line(get(gca, 'Xlim'),[bm, bm], 'color', [0, 0, 0]);  

end
