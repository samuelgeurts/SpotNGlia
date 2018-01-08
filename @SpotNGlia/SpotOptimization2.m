function objTemp = SpotOptimization2(obj, ...
    fishnumbers, ...
    inp, ...
    SpotsDetected, ...
    name)
%%dropboc spotoptimization based on minimizing RMSD
% first for a single object
% than apply in SNGbatch
% than apply
%obj            SpotNGlia object
%fishnumber     fishnumbers to proces, if input is nothing than all are processed
%inp            structure with cell with parameters to test: see SpotOptimizationScript.m
%inp.brain_sw       ['Annotation','Computation','Correction']
%spot_sw        ['Annotation','Counts']
%               Annotation: when spot and brain annotations are known, optimized on fscore
%               Computation: use the computed brain en optimize on RMSD
%               Correction: use corrected braindata in checkout.mat and optimize on RMSD
%Example
%   objTemp = obj.SpotOptimizationRMSD([],inp,'Correction')
%minspotsize:   0-20-60
%maxspotsize    200-300-470
%MinProbability 0.03-0.06-0.12
%nfishes is the number of fishes to process
%obj.nfishes is the total number of available fishes (nfishes <= obj.nfishes) = true


%{
           fishnumbers = 1:5
           name = 'SpotOptList2'
%}

%load necessary information about template and fish computations
INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo', 'BrainSegmentationInfo', 'checkup');
TEMPLATE = load([obj.SourcePath, '/', 'Template3dpf.mat'], 'ref_temp', 'SVAP_index', 'SpotVectorArrayProbability');

%create new optimization list and save if not exist already
if exist([obj.SavePath, '/', name, '.mat'], 'file')
    load([obj.SavePath, '/', name, '.mat'], 'SpotOptList2')
else
    SpotOptList2 = [];
    save([obj.SavePath, '/', name, '.mat'], 'SpotOptList2')
end

%parameters to determine on which fishes optimization takes place
if ~exist('fishnumbers', 'var') || isempty(fishnumbers)
    fishnumbers = 1:numel(obj.StackInfo);
elseif max(fishnumbers) > numel(obj.StackInfo)
    error('at least one fish does not exist in StackInfo')
end
nfishes = numel(fishnumbers);


%stores every combination of indices given by the size of the variables above
[a, b, c, d, e] = ndgrid(1:numel(inp.ColorToGrayVectorL), ...
    1:numel(inp.ScaleBaseL), ...
    1:numel(inp.KthresholdL), ...
    1:numel(inp.MPlevelsL), ...
    1:numel(inp.MPthresholdL));

%create different variable-sets to test select only spotinfo parameters
ZFParametersTemp = cell(1, numel(a));
for k2 = 1:numel(a)
    ZFParametersTemp{k2} = obj.ZFParameters(strcmp({obj.ZFParameters.stage}, 'SpotDetection'));
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'SpotSelection', 'ComputeOnSlice', inp.ComputeOnSlice, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'SpotSelection', 'SpotDistLimit', inp.SpotDistLimit, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'RgbToGray', 'ColorToGrayVector', inp.ColorToGrayVectorL{a(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'Wavelet', 'ScaleBase', inp.ScaleBaseL{b(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'Wavelet', 'Kthreshold', inp.KthresholdL{c(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'MultiProduct', 'MPlevels', inp.MPlevelsL{d(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'MultiProduct', 'MPthreshold', inp.MPthresholdL{e(k2)}, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'SpotSelection', 'MinSpotSize', 1, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'SpotSelection', 'MaxSpotSize', 1000, ''); %color selection
    ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'SpotSelection', 'MinProbability', 0, ''); %color selection
end


if ~exist('SpotsDetected', 'var') || isempty(SpotsDetected)
    %load Midbrain and Images
    %[MidBrain, AlignedFish] = loadFishAndBrain(inp);
    [MidBrain, AlignedFish] = loadFishAndBrain(obj, fishnumbers, INFO,TEMPLATE, inp.brain_sw, inp.ComputeOnSlice);

    %% SpotDetection
    SpotN = NaN(numel(ZFParametersTemp), obj.nfishes);
    SpotsDetected = cell(numel(ZFParametersTemp), obj.nfishes);
    
    %{
    for k2 = 1:numel(ZFParametersTemp)
        for k1 = 1:nfishes
            fn = fishnumbers(k1);
            fprintf('%d ', fn)
            if inp.ComputeOnSlice == 1
                [SpotsDetected{k2, fn}, ~, ~, SpotN(k2, fn)] = SpotDetectionSliceSNG(AlignedFish(fn, 1:obj.StackInfo(fn).stacksize), TEMPLATE, MidBrain{fn}, ZFParametersTemp{k2});
            else
                [SpotsDetected{k2, fn}, ~, ~, SpotN(k2, fn)] = SpotDetectionSNG(AlignedFish{fn}, TEMPLATE, MidBrain{fn}, ZFParametersTemp{k2});
            end
        end
    end
    %}
    %% parfor alternative
    alf = cell(1, nfishes);
    mdb = cell(1, nfishes);
    for k0 = 1:nfishes
        fn = fishnumbers(k0);
        if inp.ComputeOnSlice == 1
            alf{k0} = AlignedFish(fn, 1:obj.StackInfo(fn).stacksize);
        else
            alf{k0} = AlignedFish{fn};          
        end
        mdb{k0} = MidBrain{fn};
        zfp = ZFParametersTemp{k2};
        cos = inp.ComputeOnSlice;
    end
    parfor k0 = 1:nfishes
        disp(num2str(k0))
        if cos
            [SpotsDetected{k0}, ~, ~, SpotN(k0)] = SpotDetectionSliceSNG(alf{k0}, TEMPLATE, mdb{k0}, zfp);
        else
            [SpotsDetected{k0}, ~, ~, SpotN(k0)] = SpotDetectionSNG(alf{k0}, TEMPLATE, mdb{k0}, zfp);
        end
    end
    for k0 = 1:nfishes
        fn = fishnumbers(k0);
        SpotsDetected{k2, fn} = SpotsDetected{k0};
        SpotN(k2, fn) = SpotN(k0);
    end    
end

%% precompute spot threshold parameters





















%% compute spot parameters for certain ZFParameter input


[f, g, h] = ndgrid(1:numel(inp.MinSpotSizeL), ...
    1:numel(inp.MaxSpotSizeL), ...
    1:numel(inp.MinProbabilityL));

objTemp = obj;
objTemp.BrainInfo = [];
objTemp.BrainStats = [];
objTemp.SpotInfo = [];
objTemp.SpotStats = [];
objTemp.SpotBrainInfo = [];
objTemp.SpotBrainStats = [];


%% Selection of SpotsDetectedTemp based on thresholds
for k2 = 1:numel(ZFParametersTemp)
    %% Optimize color and size thresholds on RMSD and ttest
    disp('k2')
    
    meanRMSDlist = zeros(1, numel(f));
    ttestlist = zeros(1, numel(f));
    meanf1list = zeros(1, numel(f));
    for k3 = 1:numel(f)
        %%% assign new threshold logical array to obj.SpotParameters
        SpotsDetectedTemp = SpotsDetected;
        for k1 = 1:nfishes
            fn = fishnumbers(k1);
            SpotsDetectedTemp{k2, fn} = SpotsDetectedTemp{k2, fn}([SpotsDetectedTemp{k2, fn}.Area] >= inp.MinSpotSizeL(f(k3)));
            SpotsDetectedTemp{k2, fn} = SpotsDetectedTemp{k2, fn}([SpotsDetectedTemp{k2, fn}.Area] <= inp.MaxSpotSizeL(g(k3)));
            SpotsDetectedTemp{k2, fn} = SpotsDetectedTemp{k2, fn}([SpotsDetectedTemp{k2, fn}.ColorProbability] >= inp.MinProbabilityL(h(k3)));
            SpotN(k2, fn) = numel(SpotsDetectedTemp{k2, fn});
        end
        
        if strcmp(inp.spot_sw, 'Counts')
            objTemp = objTemp.TtestVal(SpotN(k2, :));
            meanRMSDlist(k3) = objTemp.SpotBrainStats.RMSD;
            ttestlist(k3) = objTemp.SpotBrainStats.ttestval;
        end
        if strcmp(inp.spot_sw, 'Annotation')
            objTemp = objTemp.SpotVal(SpotsDetectedTemp(k2, :));
            meanf1list(k3) = mean([objTemp.SpotInfo(fishnumbers).F1score]);
        end
    end
    
    if strcmp(inp.spot_sw, 'Counts')
        [mx, indx] = min(meanRMSDlist);
        %[mx, indx] = max(ttestlist);
    end
    if strcmp(inp.spot_sw, 'Annotation')
        [mx, indx] = max(meanf1list);
    end
    
    %compute the best values again
    SpotsDetectedTemp = SpotsDetected;
    for k1 = 1:nfishes
        fn = fishnumbers(k1);
        SpotsDetectedTemp{k2, fn} = SpotsDetectedTemp{k2, fn}([SpotsDetectedTemp{k2, fn}.Area] >= inp.MinSpotSizeL(f(indx)));
        SpotsDetectedTemp{k2, fn} = SpotsDetectedTemp{k2, fn}([SpotsDetectedTemp{k2, fn}.Area] <= inp.MaxSpotSizeL(g(indx)));
        SpotsDetectedTemp{k2, fn} = SpotsDetectedTemp{k2, fn}([SpotsDetectedTemp{k2, fn}.ColorProbability] >= inp.MinProbabilityL(h(indx)));
        SpotN(k2, fn) = numel(SpotsDetectedTemp{k2, fn});
    end
    
    objTemp = objTemp.TtestVal(SpotN(k2, :));
    
    if strcmp(inp.spot_sw, 'Counts')
        StoreSpotOpt(objTemp, indx, 'RMSD');
        %StoreSpotOpt(objTemp, mx, 'ttest');
        disp(['RMSD: ', num2str(mx)]);
    end
    
    if strcmp(inp.spot_sw, 'Annotation')
        objTemp = objTemp.SpotVal(SpotsDetectedTemp(k2, :));
        StoreSpotOpt(objTemp, indx, 'f1score');
        disp(['f1score: ', num2str(mx)]);
    end
end


%%


%     function [MidBrain, AlignedFish] = loadFishAndBrain(obj,fishnumbers,INFO)
%         % Choose which brain is used
%         MidBrain = cell(1, obj.nfishes);
%         switch inp.brain_sw
%             case 'Annotation'
%                 for k1 = 1:nfishes
%                     fn = fishnumbers(k1);
%                     MidBrain{fn} = fliplr([obj.Annotations(fn).MidBrain]);
%                 end
%             case 'Computation'
%                 for k1 = 1:nfishes
%                     fn = fishnumbers(k1);
%                     MidBrain{fn} = fliplr(INFO.BrainSegmentationInfo(fn).BrainEdge);
%                 end
%             case 'Correction'
%                 include = boolean([INFO.checkup.Include]);
%                 fishnumbers = fishnumbers(include(fishnumbers)); %removes excluded number according checkup.include, slim he :)
%                 nfishes = numel(fishnumbers);
%                 for k1 = 1:nfishes
%                     fn = fishnumbers(k1);
%                     MidBrain{fn} = fliplr(INFO.checkup(fn).Midbrain);
%                 end
%             otherwise
%                 error('choose between "Annotation", "Computation", "Correction"')
%         end
%         
%         % Choose if optimization is on Slice or on Combination and preload images
%         if inp.ComputeOnSlice == 1
%             load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo')
%             for k1 = 1:nfishes
%                 fn = fishnumbers(k1);
%                 
%                 tform_complete = INFO.RegistrationInfo{fn}(strcmp({INFO.RegistrationInfo{fn}.name}, 'tform_complete')).value;
%                 CorrectedFish = sng_openimstack2([obj.SavePath, '/', 'CorrectedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
%                 AlignedSlice = cell(1, numel(CorrectedFish));
%                 for k3 = 1:numel(CorrectedFish)
%                     AlignedFish{fn, k3} = imwarp(CorrectedFish{k3}, tform_complete, 'FillValues', 255, 'OutputView', TEMPLATE.ref_temp);
%                 end
%             end
%         else
%             store AlignedFish in memory for speed
%             AlignedFish = cell(1, obj.nfishes);
%             for k1 = 1:nfishes
%                 fn = fishnumbers(k1);
%                 AlignedFish{fn} = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
%             end
%         end
%     end

    function StoreSpotOpt(Object2Store, index, optimode)
        
        %if nargin == 0
        %    objTemp = minobjTemp;
        %end
        
        MinSpotSize = inp.MinSpotSizeL(f(index));
        MaxSpotSize = inp.MaxSpotSizeL(g(index));
        MinProbability = inp.MinProbabilityL(h(index));
        
        if strcmp(inp.spot_sw, 'Annotation')
            %SpotVal does not see difference between empty en non empty cells in SpotsDetectedTemp
            %therefore, f1score etc has to ben computed on a selection by fishnumbers. SpotVal does fill
            %the SpotInfo field while TTestVal does not fill BrainSpotInfo.
            SpotOpt.MeanF1score = mean([objTemp.SpotInfo(fishnumbers).F1score]);
            SpotOpt.MeanPrecision = mean([objTemp.SpotInfo(fishnumbers).Precision]);
            SpotOpt.MeanRecall = mean([objTemp.SpotInfo(fishnumbers).Recall]);
            
            SpotOpt.StdF1score = std([objTemp.SpotInfo(fishnumbers).F1score]);
            SpotOpt.StdPrecision = std([objTemp.SpotInfo(fishnumbers).Precision]);
            SpotOpt.StdRecall = std([objTemp.SpotInfo(fishnumbers).Recall]);
        end
        
        %SpotOpt.ttest = Object2Store.SpotBrainStats.ttest;
        SpotOpt.ttestval = Object2Store.SpotBrainStats.ttestval;
        SpotOpt.MeanAbsDifference = Object2Store.SpotBrainStats.MeanAbsDifference;
        SpotOpt.StdAbsDifference = Object2Store.SpotBrainStats.StdAbsDifference;
        SpotOpt.RMSD = Object2Store.SpotBrainStats.RMSD;
        
        %test if one of the thresholds found is on the edge of the seeking range
        if (MinSpotSize == min(inp.MinSpotSizeL)) && (MinSpotSize ~= 0) || (MinSpotSize == max(inp.MinSpotSizeL)) || ...
                (MaxSpotSize == min(inp.MinSpotSizeL)) && (MinSpotSize ~= 0) || (MinSpotSize == max(inp.MinSpotSizeL)) || ...
                (MinProbability == min(inp.MinProbabilityL)) && (MinProbability ~= 0) || (MinProbability == max(inp.MinProbabilityL))
            SpotOpt.WarningRange = true;
        else
            SpotOpt.WarningRange = false;
        end
        
        SpotOpt.optimode = optimode;
        
        SpotOpt.ColorToGrayVectorL = inp.ColorToGrayVectorL{a(k2)};
        SpotOpt.ScaleBaseL = inp.ScaleBaseL{b(k2)};
        SpotOpt.KthresholdL = inp.KthresholdL{c(k2)};
        SpotOpt.MPlevelsL = inp.MPlevelsL{d(k2)};
        SpotOpt.MPthresholdL = inp.MPthresholdL{e(k2)};
        SpotOpt.MinSpotSize = MinSpotSize;
        SpotOpt.MaxSpotSize = MaxSpotSize;
        SpotOpt.MinProbability = MinProbability;
        
        SpotOpt.ComputeOnSlice = inp.ComputeOnSlice;
        SpotOpt.SpotDistLimit = inp.SpotDistLimit;
        SpotOpt.brain_sw = inp.brain_sw;
        SpotOpt.spot_sw = inp.spot_sw;
        
        
        datet = uint16(clock);
        SpotOpt.date = [num2str(datet(2)), '-', num2str(datet(3)), '-', num2str(datet(1))];
        SpotOpt.time = [num2str(datet(4)), ':', num2str(datet(5)), ':', num2str(datet(6))];
        
        %{
           SpotOpt.Precision = [objTemp.SpotInfo.Precision];
           SpotOpt.Recall = [objTemp.SpotInfo.Recall];
           SpotOpt.F1score = [objTemp.SpotInfo.F1score];
           SpotOpt.AbsDifference = [objTemp.SpotInfo.AbsDifference];
        %}
        
        if isempty(SpotOptList2)
            SpotOptList2 = SpotOpt;
            %elseif exist('linen','var') & linen ~= 0
            %    SpotOptList2(linen) = SpotOpt;
        else
            SpotOptList2(numel(SpotOptList2)+1) = SpotOpt;
        end
        
        save([obj.SavePath, '/', name, '.mat'], 'SpotOptList2');
    end


end

