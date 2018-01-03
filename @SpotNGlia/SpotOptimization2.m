function objTemp = SpotOptimization2(obj,...
    fishnumbers,...
    zfinputlist,...
    brain_sw,...
    SpotsDetected,...
    name)
%% spotoptimization based on minimizing RMSD
% first for a single object
% than apply in SNGbatch
% than apply
%obj            SpotNGlia object
%fishnumber     fishnumbers to proces, if input is nothing than all are processed
%zfinputlist    structure with cell with parameters to test: see SpotOptimizationScript.m
%brain          ['Annotation','Computation','Correction']
%               Annotation: when spot and brain annotations are known, optimized on fscore
%               Computation: use the computed brain en optimize on RMSD
%               Correction: use corrected braindata in checkout.mat and optimize on RMSD
%Example
%   objTemp = obj.SpotOptimizationRMSD([],zfinputlist,'Correction')
%minspotsize:   0-20-60
%maxspotsize    200-300-470
%MinProbability 0.03-0.06-0.12
%nfishes is the number of fishes to process
%obj.nfishes is the total number of available fishes (nfishes <= obj.nfishes) = true

%{
brain_sw = 'Correction'
%}


%{
        fishnumbers = 1:5
    name = 'SpotOptList2'
%}


if exist([obj.SavePath, '/', name, '.mat'], 'file')
    load([obj.SavePath, '/', name, '.mat'], 'SpotOptList2')
else
    SpotOptList2 = [];
    save([obj.SavePath, '/', 'SpotOptList2', '.mat'], 'SpotOptList2')
end


if ~exist('fishnumbers', 'var') || isempty(fishnumbers)
    fishnumbers = 1:numel(obj.StackInfo);
elseif max(fishnumbers) > numel(obj.StackInfo)
    error('at least one fish does not exist in RegistrationInfo')
end
nfishes = numel(fishnumbers);


if exist('zfinputlist', 'var') && ~isempty(zfinputlist)
    ColorToGrayVectorL = zfinputlist.ColorToGrayVectorL;
    ScaleBaseL = zfinputlist.ScaleBaseL;
    KthresholdL = zfinputlist.KthresholdL;
    MPlevelsL = zfinputlist.MPlevelsL;
    MPthresholdL = zfinputlist.MPthresholdL;
    MinSpotSizeL = zfinputlist.MinSpotSizeL;
    MaxSpotSizeL = zfinputlist.MaxSpotSizeL;
    MinProbabilityL = zfinputlist.MinProbabilityL;
else
    ColorToGrayVectorL = {[0; 1; 0]};
    ScaleBaseL = {0.5};
    KthresholdL = {0};
    MPlevelsL = {5:7};
    MPthresholdL = {200};
    MinSpotSizeL = 8.4;
    MaxSpotSizeL = 458;
    MinProbabilityL = 0.066;
end


%stores every combination of indices given by the size of the variables above
[a, b, c, d, e] = ndgrid(1:numel(ColorToGrayVectorL), ...
    1:numel(ScaleBaseL), ...
    1:numel(KthresholdL), ...
    1:numel(MPlevelsL), ...
    1:numel(MPthresholdL));

%create different variable-sets to test select only spotinfo parameters
ZFParametersTemp = cell(1, numel(a));
for k2 = 1:numel(a)
    ZFParametersTemp{k2} = obj.ZFParameters(strcmp({obj.ZFParameters.stage}, 'SpotDetection'));
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'RgbToGray', 'ColorToGrayVector', ColorToGrayVectorL{a(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'Wavelet', 'ScaleBase', ScaleBaseL{b(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'Wavelet', 'Kthreshold', KthresholdL{c(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'MultiProduct', 'MPlevels', MPlevelsL{d(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'MultiProduct', 'MPthreshold', MPthresholdL{e(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'SpotSelection', 'MinSpotSize', 1, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'SpotSelection', 'MaxSpotSize', 1000, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'SpotSelection', 'MinProbability', 0, ''); %color selection
end
if ~exist('SpotsDetected', 'var')
    CompleteTemplate = LoadTemplateSNG([obj.SourcePath, '/', 'Template 3 dpf']);
    
    %% Choose which brain is used
    MidBrain = cell(1, obj.nfishes);
    switch brain_sw
        case 'Annotation'
            for k1 = 1:nfishes
                fn = fishnumbers(k1);
                MidBrain{fn} = fliplr([obj.Annotations(fn).MidBrain]);
            end
        case 'Computation'
            load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo')
            for k1 = 1:nfishes
                fn = fishnumbers(k1);
                MidBrain{fn} = fliplr(BrainSegmentationInfo(fn).BrainEdge);
            end
        case 'Correction'
            load([obj.SavePath, '/', obj.InfoName, '.mat'], 'checkup');
            include = boolean([checkup.Include]);
            fishnumbers = fishnumbers(include(fishnumbers)); %removes excluded number according checkup.include, slim he :)
            nfishes = numel(fishnumbers);
            for k1 = 1:nfishes
                fn = fishnumbers(k1);
                MidBrain{fn} = fliplr(checkup(fn).Midbrain);
            end
        otherwise
            error('choose between "Annotation", "Computation", "Correction"')
    end
    
    %store AlignedFish in memory for speed
    AlignedFish = cell(1, obj.nfishes);
    for k1 = 1:nfishes
        fn = fishnumbers(k1);
        AlignedFish{fn} = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
    end
    
    
else
    ZFParametersTemp = 1;
end

%% compute spot parameters for certain ZFParameter input
for k2 = 1:numel(ZFParametersTemp)
    if ~exist('SpotsDetected', 'var')
        
        SpotsDetected = cell(1, obj.nfishes);
        SpotN = NaN(1, obj.nfishes);
        
        for k1 = 1:nfishes
            fn = fishnumbers(k1);
            %fprintf([num2str(a(k2)), ',', num2str(b(k2)), ',', num2str(c(k2)), ',', num2str(d(k2)), ',', num2str(e(k2)), ',', num2str(k1), '\n'])
            fprintf('%d ', fn)
            [SpotsDetected{fn}, ~, ~, SpotN(fn)] = SpotDetectionSNG(AlignedFish{fn}, CompleteTemplate, MidBrain{fn}, ZFParametersTemp{k2});
        end
        %{
         fn = fishnumbers;
         parfor k1 = 1:nfishes
             %fprintf([num2str(a(k2)), ',', num2str(b(k2)), ',', num2str(c(k2)), ',', num2str(d(k2)), ',', num2str(e(k2)), ',', num2str(k1), '\n'])
             fprintf('%d ', fn(k1))
             [SpotsDetected{fn(k1)},~,~,SpotN(fn(k1))] = SpotDetectionSNG(AlignedFish{fn(k1)}, CompleteTemplate, MidBrain{fn(k1)}, ZFParametersTemp{k2});
         end
        %}
    end
    
    %{
 
 
       for k5 = 1:obj.nfishes
           temp = num2cell([SpotParameters{k5}.Area] >= MinSpotSize);
           [SpotParameters{k5}.LargerThan] = temp{:};
           temp = num2cell([SpotParameters{k5}.Area] <= MaxSpotSize);
           [SpotParameters{k5}.SmallerThan] = temp{:};
           temp = num2cell([SpotParameters{k5}.ColorProbability] >= MinProbability);
           [SpotParameters{k5}.MinProbability] = temp{:};
       end
 
 
         %% selects spot sizes and color probabilities, correct & incorrect
         AreaC = [];
         AreaI = [];
         ColorC = [];
         ColorI = [];
 
         for k5 = 1:obj.nfishes
             %Spotpar = SpotParameters{k5}([SpotParameters{k5}.Insite]); %selects only insite brain region
             SpotCom = reshape([SpotsDetected{k5}.Centroid], 2, numel(SpotsDetected{k5}))';
             [atemp] = ismember(SpotCom, obj.SpotInfo(k5).CorrectSpots);
 
             AreaC = [AreaC, [Spotpar(atemp(:, 1)).Area]];
             AreaI = [AreaI, [Spotpar(~atemp(:, 1)).Area]];
 
             ColorC = [ColorC, [Spotpar(atemp(:, 1)).ColorProbability]];
             ColorI = [ColorI, [Spotpar(~atemp(:, 1)).ColorProbability]];
         end
    %{
                                   if fig_tf
                                       figure;imagesc(Spotpar(3).Image)
                                   end
    %}
 
         %% Compute Threshold for area, larger than and smaller than!
         x = linspace(0, 1000, 1000); %histogram bins
         AreaHistCorrect = histcounts(AreaC, x); %correct spot area
         AreaHistIncorrect = histcounts(AreaI, x); %incorrect spot area
         %finds threshold, index is part of the right site of the histogram
         for j = 1:size(AreaHistCorrect, 2)
             FR(j) = sum(AreaHistIncorrect(1:j-1)) + sum(AreaHistCorrect(j:end));
             FL(j) = sum(AreaHistCorrect(1:j-1)) + sum(AreaHistIncorrect(j:end));
         end
         [mxR, indexR] = max(FR);
         MinSpotSize = x(indexR); %every area equal or larger than "LargerThan" is large enough
         [mxL, indexL] = max(FL(indexR:end));
         MaxSpotSize = x(indexL+indexR); %every area equal or larger than "SmallerThan" is small enough
    %{
                                   if fig_tf2
                                       sng_DistDifHistPlot((x(1:end-1)+((x(2)-x(1))/2)),AreaHistCorrect,AreaHistIncorrect);
                                       hold on;line([MinSpotSize MinSpotSize],get(gca,'Ylim'),'color',[0 0 0]);
                                       hold on;line([MaxSpotSize MaxSpotSize],get(gca,'Ylim'),'color',[0 0 0]);
                                   end
    %}
 
         %% Compute Threshold for color probability
         x = linspace(0, 1.2, 1000);
         ColorHistCorrect = histcounts(ColorC, x); %create bins
         ColorHistIncorrect = histcounts(ColorI, x); %
         for j = 1:size(ColorHistCorrect, 2)
             FR(j) = sum(ColorHistIncorrect(1:j-1)) + sum(ColorHistCorrect(j:end));
         end
         [mxR, indexR] = max(FR);
         MinProbability = x(indexR); %every area equal or larger than "LargerThan" is large enough
    %{
                                   if fig_tf2
                                       sng_DistDifHistPlot((x(1:end-1)+((x(2)-x(1))/2)),ColorHistCorrect,ColorHistIncorrect);
                                       hold on;line([MinProbability MinProbability],get(gca,'Ylim'),'color',[0 0 0]);
                                   end
    %}
 
         %using a machine learning technique would be an improvement
         %DS = [[AreaC';AreaI'],[ColorC';ColorI']]
         %lab = [ones(numel(AreaC),1);2*ones(numel(AreaI),1)]
         %a = prdataset(DS,lab)
         %a = a*scalem(a,'variance') %scaling
         %figure;
         %scatterd(a)
         %wp = parzenc(a,0.25)
         %plotc(wp)
         %err{l1} = a*wp*testc %error of testset c on trained classifier wb
 
    %}
    
    %% Optimize color and size thresholds
    
    
    [f, g, h] = ndgrid(1:numel(MinSpotSizeL), ...
        1:numel(MaxSpotSizeL), ...
        1:numel(MinProbabilityL));
    
    objTemp = obj;
    objTemp.BrainInfo = [];
    objTemp.BrainStats = [];
    objTemp.SpotInfo = [];
    objTemp.SpotStats = [];
    objTemp.SpotBrainInfo = [];
    objTemp.SpotBrainStats = [];
    
    minRMSD = 1000;
    maxttest = 0;
    for k3 = 1:numel(f)
        if round(k3/100) == k3 / 100; fprintf('%d ', k3); end
        SpotsDetectedTemp = SpotsDetected;
        %%% assign new threshold logical array to obj.SpotParameters
        for k1 = 1:nfishes
            fn = fishnumbers(k1);
            SpotsDetectedTemp{fn} = SpotsDetectedTemp{fn}([SpotsDetectedTemp{fn}.Area] >= MinSpotSizeL(f(k3)));
            SpotsDetectedTemp{fn} = SpotsDetectedTemp{fn}([SpotsDetectedTemp{fn}.Area] <= MaxSpotSizeL(g(k3)));
            SpotsDetectedTemp{fn} = SpotsDetectedTemp{fn}([SpotsDetectedTemp{fn}.ColorProbability] >= MinProbabilityL(h(k3)));
            SpotN(fn) = numel(SpotsDetectedTemp{fn});
        end        
        
        objTemp = objTemp.TtestVal(SpotN);
        
        %for selecting on RMSD
        meanRMSDlist = objTemp.SpotBrainStats.RMSD;
        if meanRMSDlist < minRMSD
            minRMSD = meanRMSDlist;
            %indx = k3;
            minobjTemp = objTemp;
            StoreSpotOpt(objTemp, k3, 'RMSD')
        end
        
        
        %for selecting on ttest
        ttestlist = objTemp.SpotBrainStats.ttestval;
        if ttestlist >= maxttest
            maxttest = ttestlist;
            %indx = k3;
            %minobjTemp = objTemp;
            
            if maxttest >= 1
                %optimode = 'ttest';
                StoreSpotOpt(objTemp, k3, 'ttest')
            end
        end
        %}
        % for selecting on f1score
        % objTemp = objTemp.SpotVal(SpotsDetectedTemp);
        % meanf1list(k3) = objTemp.SpotStats.MeanF1score;
    end
    clear SpotsDetected
    %optimode = 'RMSD';
    %StoreSpotOpt(minobjTemp,indx)
    disp(num2str(minobjTemp.SpotBrainStats.RMSD));
    
end

    function StoreSpotOpt(Object2Store, index, optimode)
        
        %if nargin == 0
        %    objTemp = minobjTemp;
        %end
        
        MinSpotSize = MinSpotSizeL(f(index));
        MaxSpotSize = MaxSpotSizeL(g(index));
        MinProbability = MinProbabilityL(h(index));
        
        SpotOpt.ttest = Object2Store.SpotBrainStats.ttest;
        SpotOpt.ttestval = Object2Store.SpotBrainStats.ttestval;
        SpotOpt.MeanAbsDifference = Object2Store.SpotBrainStats.MeanAbsDifference;
        SpotOpt.StdAbsDifference = Object2Store.SpotBrainStats.StdAbsDifference;
        SpotOpt.RMSD = Object2Store.SpotBrainStats.RMSD;
        
        %test if one of the thresholds found is on the edge of the seeking range
        if (MinSpotSize == min(MinSpotSizeL)) && (MinSpotSize ~= 0) || (MinSpotSize == max(MinSpotSizeL)) || ...
                (MaxSpotSize == min(MinSpotSizeL)) && (MinSpotSize ~= 0) || (MinSpotSize == max(MinSpotSizeL)) || ...
                (MinProbability == min(MinProbabilityL)) && (MinProbability ~= 0) || (MinProbability == max(MinProbabilityL))
            SpotOpt.WarningRange = true;
        else
            SpotOpt.WarningRange = false;
        end
        
        %{
    SpotOpt.MeanF1score = mean([objTemp.SpotInfo.F1score]);
    SpotOpt.MeanPrecision = mean([objTemp.SpotInfo.Precision]);
    SpotOpt.MeanRecall = mean([objTemp.SpotInfo.Recall]);
        %}
        SpotOpt.optimode = optimode;
        
        SpotOpt.ColorToGrayVectorL = ColorToGrayVectorL{a(k2)};
        SpotOpt.ScaleBaseL = ScaleBaseL{b(k2)};
        SpotOpt.KthresholdL = KthresholdL{c(k2)};
        SpotOpt.MPlevelsL = MPlevelsL{d(k2)};
        SpotOpt.MPthresholdL = MPthresholdL{e(k2)};
        SpotOpt.MinSpotSize = MinSpotSize;
        SpotOpt.MaxSpotSize = MaxSpotSize;
        SpotOpt.MinProbability = MinProbability;
        
        SpotOpt.datetime = string(datetime);
        
        %{
    SpotOpt.MeanAbsDifference = mean(abs([objTemp.SpotInfo.AbsDifference]));
    SpotOpt.MeanRelDifference = mean(abs([objTemp.SpotInfo.RelDifference]));
 
    SpotOpt.StdPrecision = std([objTemp.SpotInfo.Precision]);
    SpotOpt.StdRecall = std([objTemp.SpotInfo.Recall]);
    SpotOpt.StdF1score = std([objTemp.SpotInfo.F1score]);
    SpotOpt.StdAbsDifference = std(abs([objTemp.SpotInfo.AbsDifference]));
    SpotOpt.StdRelDifference = std(abs([objTemp.SpotInfo.RelDifference]));
 
    SpotOpt.Precision = [objTemp.SpotInfo.Precision];
    SpotOpt.Recall = [objTemp.SpotInfo.Recall];
    SpotOpt.F1score = [objTemp.SpotInfo.F1score];
    SpotOpt.AbsDifference = [objTemp.SpotInfo.AbsDifference];
    SpotOpt.RelDifference = [objTemp.SpotInfo.RelDifference];
        %}
        
        if isempty(SpotOptList2)
            SpotOptList2 = SpotOpt;
            %elseif exist('linen','var') & linen ~= 0
            %    SpotOptList2(linen) = SpotOpt;
        else
            SpotOptList2(numel(SpotOptList2)+1) = SpotOpt;
        end
        
        save([obj.SavePath, '/', 'SpotOptList2', '.mat'], 'SpotOptList2');
        
    end


end

