function ShowFishHeadHist(obj, fishnumbers)

load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo')

if ~exist('fishnumbers', 'var')
    fishnumbers = 1:numel(obj.StackInfo);
elseif max(fishnumbers) > numel(obj.StackInfo)
    error('at least one fish does not exist in RegistrationInfo')
end
nfishes = numel(fishnumbers);


for k1 = 1:nfishes
    fn = fishnumbers(k1);
    AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
    threshold = RegistrationInfo{k1}(find(strcmp({RegistrationInfo{k1}.name}, 'BackgroundThreshold'))).value;
    
    figure
    bn=1;
    
    for k3 = 1:3
        Img = AlignedFish(:, 700:end, k3);
        h = hist(Img(:), 0:bn:255);    
               
        subplot(3, 1, k3);
        %set(gca,'Xlim',[0,])
        bar(h, 'Facecolor', [0.4, 0.4, 0.8], 'Edgecolor', 'none', 'BarWidth', 1);
        hold on
        line([threshold(k3), threshold(k3)], [get(gca, 'Ylim')], 'color', [0, 0, 0]);        
        set(gca, 'XLim', [0,((255 + 1)/bn)]);
        set(gca, 'YLim', [0, bn * 8000]);
    end

end

end
