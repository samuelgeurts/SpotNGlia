
%% spotoptimization based on minimizing RMSD
% first for a single object
% than apply in SNGbatch
% than apply


function objTemp = SpotOptimization(obj, fishnumbers, zfinputlist, brain_sw)
%obj            SpotNGlia object
%fishnumber     fishnumbers to proces, if input is nothing than all are processed
%zfinputlist    structure with cell with parameters to test: see SpotOptimizationScript.m
%brain          ['Annotation','Computation','Correction']
%               Annotation: when spot and brain annotations are known, optimized on fscore
%               Computation: use the computed brain en optimize on RMSD
%               Correction: use corrected braindata in checkout.mat and optimize on RMSD

%minspotsize:   0-20-60
%maxspotsize    200-300-470
%MinProbability 0.03-0.06-0.12


%{
 fishnumbers = 1:5
%}
CompleteTemplate = LoadTemplateSNG([obj.SourcePath, '/', 'Template 3 dpf']);


if exist([obj.SavePath, '/', 'SpotOptList', '.mat'], 'file')
    load([obj.SavePath, '/', 'SpotOptList', '.mat'], 'SpotOptList')
else
    SpotOptList = [];
    save([obj.SavePath, '/', 'SpotOptList', '.mat'], 'SpotOptList')
end


if ~exist('fishnumbers', 'var') || isempty(fishnumbers)
    fishnumbers = 1:numel(obj.StackInfo);
elseif max(fishnumbers) > numel(obj.StackInfo)
    error('at least one fish does not exist in RegistrationInfo')
end

nfishes = numel(fishnumbers);

if exist('zfinputlist', 'var')
    ColorToGrayVectorL = zfinputlist.ColorToGrayVectorL;
    ScaleBaseL = zfinputlist.ScaleBaseL;
    KthresholdL = zfinputlist.KthresholdL;
    MPlevelsL = zfinputlist.MPlevelsL;
    MPthresholdL = zfinputlist.MPthresholdL;
else
    ColorToGrayVectorL = {[0; 1; 0]};
    ScaleBaseL = {0.5};
    KthresholdL = {0};
    MPlevelsL = {5:7};
    MPthresholdL = {200};
    MinSpotSizeL = {8.4084}
    MaxSpotSizeL = {458}
    MinProbabilityL = {0.066}
end


%stores every combination of indices given by the size of the variables above
[a, b, c, d, e] = ndgrid(1:numel(ColorToGrayVectorL), ...
    1:numel(ScaleBaseL), ...
    1:numel(KthresholdL), ...
    1:numel(MPlevelsL), ...
    1:numel(MPthresholdL));

%create different variable-sets to test select only spotinfo parameters
for k2 = 1:numel(a)
    ZFParametersTemp{k2} = obj.ZFParameters(strcmp({obj.ZFParameters.stage}, 'SpotDetection'));
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'RgbToGray', 'ColorToGrayVector', ColorToGrayVectorL{a(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'Wavelet', 'ScaleBase', ScaleBaseL{b(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'Wavelet', 'Kthreshold', KthresholdL{c(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'MultiProduct', 'MPlevels', MPlevelsL{d(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'MultiProduct', 'MPthreshold', MPthresholdL{e(k2)}, ''); %color selection
    %ZFParametersTemp = sng_zfinput(ZFParametersTemp,0,'SpotDetection','SpotSelection','MinSpotSize',MinSpotSizeL{g},''); %color selection
    %ZFParametersTemp = sng_zfinput(ZFParametersTemp,0,'SpotDetection','SpotSelection','MaxSpotSize',MaxSpotSizeL{h},''); %color selection
    %ZFParametersTemp = sng_zfinput(ZFParametersTemp,0,'SpotDetection','SpotSelection','MinProbability',MinProbabilityL{i},''); %color selection
end

%compute Spots


switch brain_sw
    case 'Annotion'

    case 'Computation'
        load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo')     
    case 'Correction'
        load([obj.SavePath, '/', obj.InfoName, '.mat'], 'checkup')
    otherwise
        error('choose between "Annotation", "Computation", "Correction"')
end

AlignedFish = cell(1,50);
for k1 = 1:nfishes
    fn = fishnumbers(k1);
    AlignedFish{k1} = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
end


for k2 = 1:numel(ZFParametersTemp)
    %compute spot parameters for certain ZFParameter input
    
    SpotParameters = cell(1,numel(obj.StackInfo));
    for k1 = 1:nfishes
        fn = fishnumbers(k1);
        %fprintf([num2str(a(k2)), ',', num2str(b(k2)), ',', num2str(c(k2)), ',', num2str(d(k2)), ',', num2str(e(k2)), ',', num2str(k1), '\n'])
        fprintf('%d ', k1)
                
        %choose which midbrain to use
        switch brain_sw
            case 'Annotion'
                MidBrain = fliplr([obj.Annotations(fn).MidBrain]);
            case 'Computation'
                MidBrain = fliplr(BrainSegmentationInfo(fn).BrainEdge);
            case 'Correction'
                MidBrain = checkup(fn).Midbrain;       
        end        
        [~, SpotParameters{k1}] = SpotDetectionSNG(AlignedFish{k1}, CompleteTemplate, fliplr(MidBrain), ZFParametersTemp{k2});        
    end
    
    %% assign new threshold logical array to obj.SpotParameters
        
    if isfield(zfinputlist, 'MinSpotSizeL') && ...
            isfield(zfinputlist, 'MaxSpotSizeL') && ...
            isfield(zfinputlist, 'MinProbabilityL')
        MinSpotSize = zfinputlist.MinSpotSizeL{1};
        MaxSpotSize = zfinputlist.MaxSpotSizeL{1};
        MinProbability = zfinputlist.MinProbabilityL{1};
        
        for k5 = 1:nfishes
            temp = num2cell([SpotParameters{k5}.Area] >= MinSpotSize);
            [SpotParameters{k5}.LargerThan] = temp{:};
            temp = num2cell([SpotParameters{k5}.Area] <= MaxSpotSize);
            [SpotParameters{k5}.SmallerThan] = temp{:};
            temp = num2cell([SpotParameters{k5}.ColorProbability] >= MinProbability);
            [SpotParameters{k5}.MinProbability] = temp{:};
        end
        objTemp = obj;
        
        f = 1; h = 1; g = 1; indx = 1;
        
    else
        
        %% selects spot sizes and color probabilities, correct & incorrect
        AreaC = [];
        AreaI = [];
        ColorC = [];
        ColorI = [];
        for k5 = 1:nfishes
            Spotpar = SpotParameters{k5}([SpotParameters{k5}.Insite]); %selects only insite brain region
            SpotCom = reshape([Spotpar.Centroid], 2, numel(Spotpar))';
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
        
        %% Optimize color and size thresholds
        
        
        MinSpotSizeT = MinSpotSize * zfinputlist.MinSpotRange;
        MaxSpotSizeT = MaxSpotSize * zfinputlist.MaxSpotRange;
        MinProbabilityT = MinProbability * zfinputlist.MinProbabilityRange;
        
        [f, g, h] = ndgrid(1:numel(MinSpotSizeT), ...
            1:numel(MaxSpotSizeT), ...
            1:numel(MinProbabilityT));
        
        meanf1list = zeros(1, numel(f));
        
        objTemp = obj;
        objTemp.BrainInfo = [];
        objTemp.BrainStats = [];
        objTemp.SpotInfo = [];
        objTemp.SpotStats = [];
        objTemp.SpotBrainInfo = [];
        objTemp.SpotBrainStats = [];
        
        for k3 = 1:numel(f)
            %%% assign new threshold logical array to obj.SpotParameters
            for k5 = 1:nfishes
                temp = num2cell([SpotParameters{k5}.Area] >= MinSpotSizeT(f(k3)));
                [SpotParameters{k5}.LargerThan] = temp{:};
                temp = num2cell([SpotParameters{k5}.Area] <= MaxSpotSizeT(g(k3)));
                [SpotParameters{k5}.SmallerThan] = temp{:};
                temp = num2cell([SpotParameters{k5}.ColorProbability] >= MinProbabilityT(h(k3)));
                [SpotParameters{k5}.MinProbability] = temp{:};
            end
            objTemp = objTemp.SpotVal(SpotParameters);
            meanf1list(k3) = objTemp.SpotStats.MeanF1score;
        end
        [mx, indx] = max(meanf1list);
        
        for k5 = 1:nfishes
            temp = num2cell([SpotParameters{k5}.Area] >= MinSpotSizeT(f(indx)));
            [SpotParameters{k5}.LargerThan] = temp{:};
            temp = num2cell([SpotParameters{k5}.Area] <= MaxSpotSizeT(g(indx)));
            [SpotParameters{k5}.SmallerThan] = temp{:};
            temp = num2cell([SpotParameters{k5}.ColorProbability] >= MinProbabilityT(h(indx)));
            [SpotParameters{k5}.MinProbability] = temp{:};
            
            MinSpotSize = MinSpotSizeT(f(indx));
            MaxSpotSize = MaxSpotSizeT(g(indx));
            MinProbability = MinProbabilityT(h(indx));
        end
        
        
    end
    
    objTemp = objTemp.SpotVal(SpotParameters);
    
    
    disp(num2str(mean([objTemp.SpotInfo.F1score])));
    
    %disp([num2str(l2),' ',num2str(l1),' ',num2str(mean([objTemp.SpotInfo.F1score]))]);
    
    SpotOpt.MeanF1score = mean([objTemp.SpotInfo.F1score]);
    SpotOpt.MeanPrecision = mean([objTemp.SpotInfo.Precision]);
    SpotOpt.MeanRecall = mean([objTemp.SpotInfo.Recall]);
    
    if (f(indx) == 1 && zfinputlist.MinSpotRange(1) ~= 0) || f(indx) == numel(MinSpotSizeT) || ...
            (g(indx) == 1 && zfinputlist.MaxSpotRange(1) ~= 0) || g(indx) == numel(MaxSpotSizeT) || ...
            (h(indx) == 1 && zfinputlist.MinProbabilityRange(1) ~= 0) || h(indx) == numel(MinProbabilityT)
        SpotOpt.WarningRange = true;
    else
        SpotOpt.WarningRange = false;
    end
    
    SpotOpt.ColorToGrayVectorL = ColorToGrayVectorL{a(k2)};
    SpotOpt.ScaleBaseL = ScaleBaseL{b(k2)};
    SpotOpt.KthresholdL = KthresholdL{c(k2)};
    SpotOpt.MPlevelsL = MPlevelsL{d(k2)};
    SpotOpt.MPthresholdL = MPthresholdL{e(k2)};
    SpotOpt.MinSpotSize = MinSpotSize;
    SpotOpt.MaxSpotSize = MaxSpotSize;
    SpotOpt.MinProbability = MinProbability;
    SpotOpt.MinSpotRange = zfinputlist.MinSpotRange;
    SpotOpt.MaxSpotRange = zfinputlist.MaxSpotRange;
    SpotOpt.MinProbabilityRange = zfinputlist.MinProbabilityRange;
    
    SpotOpt.MinSpotFactor = f(indx);
    SpotOpt.MaxSpotFactor = g(indx);
    SpotOpt.MinProbabilityFactor = h(indx);
    
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
    
    SpotOpt.date = date;
    
    if isempty(SpotOptList)
        SpotOptList = SpotOpt;
        %elseif exist('linen','var') & linen ~= 0
        %    SpotOptList(linen) = SpotOpt;
    else
        SpotOptList(numel(SpotOptList)+1) = SpotOpt;
    end
    
    save([obj.SavePath, '/', 'SpotOptList', '.mat'], 'SpotOptList');
end

end
