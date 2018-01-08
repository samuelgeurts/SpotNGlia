function [MidBrain, AlignedFish] = loadFishAndBrain(obj, fishnumbers, INFO, TEMPLATE, brain_sw, ComputeOnSlice)
%load midbrains and alignedfishes
%obj:   SpotNGlia object
%INFO:  INFO.(BrainSegmentation,checkup,RegistrationInfo)
%TEMPLATE:  TEMPLATE.ref_temp
%brain_sw:  'Annotation','Computation','Correction'
%ComputeonSlice: [1,0]
%EXAMPLE: [MidBrain, AlignedFish] = loadFishAndBrain(obj, fishnumbers, INFO,TEMPLATE, 'Annotation', 1)

if ~exist('fishnumbers', 'var') || isempty(fishnumbers)
    fishnumbers = 1:numel(obj.StackInfo);
elseif max(fishnumbers) > numel(obj.StackInfo)
    error('at least one fish does not exist in StackInfo')
end
nfishes = numel(fishnumbers);


%% Choose which brain is used
MidBrain = cell(1, obj.nfishes);
switch brain_sw
    case 'Annotation'
        for k1 = 1:nfishes
            fn = fishnumbers(k1);
            MidBrain{fn} = fliplr([obj.Annotations(fn).MidBrain]);
        end
    case 'Computation'
        if ~isfield(INFO, 'BrainSegmentationInfo')
            INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo');
        end
        for k1 = 1:nfishes
            fn = fishnumbers(k1);
            MidBrain{fn} = fliplr(INFO.BrainSegmentationInfo(fn).BrainEdge);
        end
    case 'Correction'
        if ~isfield(INFO, 'checkup')
            INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'checkup');
        end
        include = logical([INFO.checkup.Include]);
        fishnumbers = fishnumbers(include(fishnumbers)); %removes excluded number according checkup.include, slim he :)
        nfishes = numel(fishnumbers);
        for k1 = 1:nfishes
            fn = fishnumbers(k1);
            MidBrain{fn} = fliplr(INFO.checkup(fn).Midbrain);
        end
    otherwise
        error('choose between "Annotation", "Computation", "Correction"')
end

%% Choose if optimization is on Slice or on Combination and preload images
if ComputeOnSlice
    if ~isfield(INFO, 'checkup')
        INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo');
    end
    for k1 = 1:nfishes
        fn = fishnumbers(k1);
        
        tform_complete = INFO.RegistrationInfo{fn}(strcmp({INFO.RegistrationInfo{fn}.name}, 'tform_complete')).value;
        CorrectedFish = sng_openimstack2([obj.SavePath, '/', 'CorrectedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
        %AlignedSlice = cell(1, numel(CorrectedFish));
        for k3 = 1:numel(CorrectedFish)
            AlignedFish{fn, k3} = imwarp(CorrectedFish{k3}, tform_complete, 'FillValues', 255, 'OutputView', TEMPLATE.ref_temp);
        end
    end
else
    %store AlignedFish in memory for speed
    AlignedFish = cell(1, obj.nfishes);
    for k1 = 1:nfishes
        fn = fishnumbers(k1);
        AlignedFish{fn} = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
    end
end
end