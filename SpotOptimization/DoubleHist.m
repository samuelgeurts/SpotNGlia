    function DoubleHist(f1, f2, titl)
    %subfunction for DefineSpotThreshold
%{        
f1 = SpotMeanhsvC(:, 2)
        f2 = SpotMeanhsvI(:, 2)
        
        f1 = SpotMeanC(:, 3)
        f2 = SpotMeanI(:, 3)
        
        f1 = AzimC
        f2 = AzimI
        
        f1 = ElevC
        f2 = ElevI
        
        f1 = ColorC %<- overfitted
        f2 = ColorI
        
        f1 = RC
        f2 = RI
        
        f1 = AreaC
        f2 = AreaI
        %}
        
        x = linspace(min([f1(:, 1); f2(:, 1)]), max([f1(:, 1); f2(:, 1)]), 100); %histogram bins
        hc = histcounts(f1(:, 1), x); %correct spot area
        hi = histcounts(f2(:, 1), x); %incorrect spot area
        sng_DistDifHistPlot((x(1:end-1) + ((x(2) - x(1)) / 2)), hc, hi);
        title(titl)
        
        %hold on;line([MinSpotSize MinSpotSize],get(gca,'Ylim'),'color',[0 0 0]);
        %hold on;line([MaxSpotSize MaxSpotSize],get(gca,'Ylim'),'color',[0 0 0]);
    end
