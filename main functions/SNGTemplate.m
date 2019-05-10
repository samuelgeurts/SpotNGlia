classdef SNGTemplate < handle
    % class for generating template object. Based on the function MakeTemplateMat
    % which it replaces. Running the function load the needed template
    % variables and creates the template on a chosen location. In future the load function has to be replaced by the functions
    % really generating the template instead of loading.
    
    properties
        fileName = 'Template3dpf'
        sourcePath
        TemplatePath3dpf
        RoiBrainPath
        RoiSpotPath
        savePath
        
        %Object
        SpotNGliaObject
        
        %Template
        Template
        Size
        ref_temp
        
        %Midbrain
        CenterMidBrain
        %MidBrainDistanceMap
        BandMidBrain
        polarMidbrainBandWithGaussianWithBlurr
        %Edge Template
        %EdgeTemplate1
        EdgeFilter
        
        %Spot
        SpotVectorArrayProbability
        SVAP_index
        
        spotcolorsTT;
        backrcolorsTT;
        
    end %basic
    properties(Constant = true) %inputparameters
        %these input parameters maybe should be added to the main input
        %structure. But because it is a template variable only maybe to a
        %separete variable. -maybe, i have to think about it
        
        %the sigma based variable of the gaussian added to the polar transform of the
        %probability band of the midbrains. Used to increase the probability
        %area without changing the basic probability area based on the 50 fish
        %batch
        sigmabasedVariable = 6;
        extraSmoothingKernel = [5, 5];
        filterwidth = 21; %width of the average edge template
        
        
        %belonging to spot parameters
        CubeSizeSpotErode = 4; %spot structure element for erode
        CubeSizeSpotDilate = 4; %background structure element for dilate
        winsiz = 20; %size of spot box is 2*winsiz + 1 = 51x51
    end
    properties(Transient = true, Hidden = true) %additional output
        % Eyes
        EyeRegion1
        EyeRegion2
        EyeMask1
        EyeMask2
        EyeCenter1
        EyeCenter2
        
        % MidBrain
        MidbrainInfo = struct('Poly', [], 'Area', [], 'Centroid', [], 'BoundingBox', []);
        MeanMidBrain
        MidBrainDistanceMap
        
        distancemapblurm
        skeletonm
        midbrainXCoordList
        midbrainYCoordList
        midbrainInnerContour
        midbrainOuterContour
        
        squareMidbrainBand
        polarMidbrainBand
        polarMidbrainBandWithGaussian
        listGaussians
        
        EdgeTemplate1;
        
        
        % ForeBrain
        ForbrainInfo = struct('Poly', [], 'Area', [], 'Centroid', [], 'BoundingBox', []);
        MeanForeBrain
        ForeBrainDistanceMap
        BandForeBrain
        distancemapblurf
        skeletonf
        forebrainXCoordList
        forebrainYCoordList
        
        % Spot
        SpotContrastVector
        Classifier
        Spot1
        SpotF
        BackF
        Mask2
    end
    properties(Transient = true) %image output
        % Brain Images
    end
    properties
        %input properties from object
        ColorToGrayVector = [];
        ScaleBase
        ScaleLevels
        Kthreshold
        MPlevels
        MPthreshold
        MinSpotSize
        MaxSpotSize
        MinProbability
        SpotDistLimit
        ComputeOnSlice
    end
    
    
    methods
        function objt = SNGTemplate(SpotNGliaObject)
            disp('select template source path');
            objt.sourcePath = uigetdir(pwd, 'select template source path');
            %folder with all 3dpf template information
            objt.TemplatePath3dpf = [objt.sourcePath, '/', 'Template 3 dpf'];
            %folder with all brain rois
            objt.RoiBrainPath = [objt.sourcePath, filesep, 'Roi brain'];
            %folder with all spot rois
            objt.RoiSpotPath = [objt.sourcePath, filesep, 'Roi microglia'];
            %folder to save
            objt.savePath = uigetdir(pwd, 'select template save path');
            
            loadSpotNGliaObject;
            loadInputParameters;
        end
        function loadSpotNGliaObject(objt, SpotNGliaObject)
            %added for the TemplateSpot function
            objt.SpotNGliaObject = SpotNGliaObject;
        end
        function loadInputParameters(objt)
            inputFields = fields(objt.SpotNGliaObject.SngInputParameters.SpotDetection);
            for iInputField = 1:numel(inputFields)
                objt.(inputFields{iInputField}) = objt.SpotNGliaObject.SngInputParameters.SpotDetection.(inputFields{iInputField});
            end
        end
        
        function value = get.midbrainXCoordList(objt)
            %caller function to run BrainTemplate
            if isempty(objt.midbrainXCoordList)
                disp('compute Brain parameters');
                objt.BrainTemplate; %compute all Brain paramaters including midbrainXCoordList
            end
            value = objt.midbrainXCoordList;
        end
        function value = get.squareMidbrainBand(objt)
            %caller function to run polarProbabilityMidbrain
            %asumed is that checking for squareMidbrainBand is enough
            if isempty(objt.squareMidbrainBand)
                disp('compute Brain parameters');
                objt.polarProbabilityMidbrain; %compute polar band
            end
            value = objt.squareMidbrainBand;
        end
        function value = get.polarMidbrainBandWithGaussian(objt)
            %caller function to run polarProbabilityMidbrain
            %asumed is that checking for squareMidbrainBand is enough
            if isempty(objt.polarMidbrainBandWithGaussian)
                disp('compute Brain parameters');
                objt.polarProbabilityMidbrain; %compute polar band
            end
            value = objt.polarMidbrainBandWithGaussian;
        end
        function value = get.Spot1(objt)
            if isempty(objt.Spot1)
                disp('compute additional Spot parameters');
                objt.TemplateSpot; %compute polar band
            end
            value = objt.Spot1;
        end
        
        function loadTemplate(objt)
            % Mean Fish Registration Template
            %based on older template, can be updated with 50 fish batch but lead
            
            objt.Template = imread([objt.TemplatePath3dpf, '/template.tif']);
            SIZE = size(objt.Template);
            objt.Size = SIZE;
            objt.ref_temp = imref2d(SIZE);
        end
        function loadEyeParameters(objt)
            %Eye region, annotated by Anouk H
            %annotated on the older template
            temp = sng_html2roi2([objt.TemplatePath3dpf, '/edof_eyes.html']);
            objt.EyeRegion1 = temp{1}{2};
            objt.EyeRegion2 = temp{1}{1};
            
            eye1 = objt.EyeRegion2;
            eye2 = objt.EyeRegion1;
            
            objt.EyeMask1 = poly2mask(eye1(:, 1), eye1(:, 2), SIZE(1), SIZE(2));
            objt.EyeMask2 = poly2mask(eye2(:, 1), eye2(:, 2), SIZE(1), SIZE(2));
            clear eye1 eye2
            
            [x1c, y1c] = sng_CenterOfMass(objt.EyeMask1);
            [x2c, y2c] = sng_CenterOfMass(objt.EyeMask2);
            
            objt.EyeCenter1 = [x1c, y1c];
            objt.EyeCenter2 = [x2c, y2c];
            
            %{
              figure;imagesc(eyem1+eyem2)
            %}
        end
        function loadBrainEdge(objt)
            % Brain Edge Template
            % has te be removed later on
            %a more normal function would be better
            temp = load([objt.TemplatePath3dpf, '/EdgeTemplate.mat']); %polar transforms and edge template
            objt.EdgeTemplate1 = temp.EdgeTemplate;
            
            
        end
        function loadSpotParameters(objt)
            % Spot Parameters
            % generated with SpotTemplate2
            temp = load([objt.TemplatePath3dpf, '/SpotTemplate.mat'], 'SpotTemplateVar');
            objt.SpotContrastVector = temp.SpotTemplateVar.SpotContrastVector;
            objt.SpotVectorArrayProbability = temp.SpotTemplateVar.SpotVectorArrayProbability;
            objt.SVAP_index = temp.SpotTemplateVar.SVAP_index;
            %clear SpotTemplateVar
            
        end
        function loadMachineLearningParameters(objt)
            %Based on the function MakeTemplateMat. Running the function load the needed template
            %variables and creates the template on a chosen location. In future the load function has to be replaced by the functions
            %really generating the template instead of loading.
            
            
            % Machine learning SpotParameters
            ClassifierInfo = load([objt.sourcePath, filesep, 'Microglia classifiers', filesep, 'pksvc5050.mat'], 'wp');
            objt.Classifier = ClassifierInfo.wp;
        end
        function loadBrainParameters(objt)
            % Brain
            %generated with the function TemplateBrainParameters
            
            %Midbrain info (poly,area,centroid,box) of individual fishes
            temp = load([objt.TemplatePath3dpf, '/MidbrainInfo.mat']);
            objt.MidbrainInfo = temp.MidbrainInfo;
            
            %Forbrain info (poly,area,centroid,box) of individual fishes
            temp = load([objt.TemplatePath3dpf, '/ForbrainInfo.mat']);
            objt.ForbrainInfo = temp.ForbrainInfo;
            
            %MeanMidbrain
            temp = load([objt.TemplatePath3dpf, '/MeanMidbrain.mat']);
            objt.MeanMidBrain = temp.MeanMidBrain;
            objt.CenterMidBrain = temp.Centroidm;
            objt.MidBrainDistanceMap = temp.distancemapm;
            objt.BandMidBrain = temp.bandm;
            %Meanforbrain makes not much sense as it variates largely
            %it must be dependt on the midbrean and eyes
        end
        function loadAll(objt)
            %Based on the function MakeTemplateMat. Running the function load the needed template
            %variables and creates the template on a chosen location. In future the load function has to be replaced by the functions
            %really generating the template instead of loading.
            
            objt.loadBrainEdge;
            objt.loadBrainParameters;
            objt.loadEyeParameters;
            objt.loadMachineLearningParameters;
            objt.loadSpotParameters;
            objt.loadTemplate;
        end
        
        function save(objt)
            %save(strcat(objt.sourcePath, filesep, objt.fileName, '.mat'), 'objt');
            save(strcat(objt.savePath, filesep, objt.fileName, '.mat'), 'objt');
        end
        
        function BrainTemplate(objt)
            % function based on TemplateBrainParameters
            % this function needs aan SpotNGlia object to retrieve information from, stackinfo information
            % transform raw midbrain coords to horizontal fish coords and create skeleton
            % and distancemap. Also other values are computed.
            % based on GenerateSkeletonDistancemapFromRawCoordinates3 specially applied
            % on the 20170327_3dpf batch of 0 fishes
            
            %subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.03], [0.04 0.03], [0.03 0.03]);
            
            
            %obj.LoadParameters('RegistrationInfo');
            %objt.SpotNGliaObject.LoadParameters('RegObject');
            
            tformList = [objt.SpotNGliaObject.RegObject.tform_1234];
            
            SIZE = objt.Size; %the size of the Template
            
            
            %% unzip and open roi files
            roiFiles = dir(fullfile(objt.RoiBrainPath, '*.zip'));
            nRoiFiles = numel(roiFiles);
            
            ROI = cell(nRoiFiles, 2);
            for iRoiFile = 1:nRoiFiles
                [~, name, ~] = fileparts(roiFiles(iRoiFile).name);
                
                
                if ~exist([objt.RoiBrainPath, '/', name], 'dir')
                    %creates unzipped roi in a folder with the name of the zipfile
                    a = unzip([objt.RoiBrainPath, filesep, roiFiles(iRoiFile).name], [objt.RoiBrainPath, '/', name]);
                else
                    a = sng_FilesFromMap3([objt.RoiBrainPath, '/', name], 'roi');
                end
                ROI{iRoiFile, 1} = ReadImageJROI(a{1}); %midbrain? no variates
                ROI{iRoiFile, 2} = ReadImageJROI(a{2}); %forbrain? no variates
            end
            
            
            %% corrects for htmls which for and midbrains are in reversed order i.e. place largest area at first
            for iRoiFile = 1:nRoiFiles
                %unsorted
                %figure;mapshow(ROI{k,1}.mnCoordinates(:,1),ROI{k,1}.mnCoordinates(:,2))
                
                BW1 = poly2mask(ROI{iRoiFile, 1}.mnCoordinates(:, 1), ROI{iRoiFile, 1}.mnCoordinates(:, 2), SIZE(2), SIZE(1));
                BW2 = poly2mask(ROI{iRoiFile, 2}.mnCoordinates(:, 1), ROI{iRoiFile, 2}.mnCoordinates(:, 2), SIZE(2), SIZE(1));
                if sum(BW2(:)) > sum(BW1(:))
                    temp = ROI{iRoiFile, 1};
                    ROI{iRoiFile, 1} = ROI{iRoiFile, 2};
                    ROI{iRoiFile, 2} = temp;
                end
                %sorted s.t. largest first (midbrain)
                %figure;mapshow(ROI{iRoiFile,1}.mnCoordinates(:,1),ROI{iRoiFile,1}.mnCoordinates(:,2))
            end
            
            %% Computes mean skeleton and band of batch + stores individual brain information
            
            %preallocation
            
            rpm(50) = struct('Area', [], 'Centroid', [], 'BoundingBox', []);
            rpf(50) = struct('Area', [], 'Centroid', [], 'BoundingBox', []);
            
            midbrainInnerArea = true(SIZE(1), SIZE(2));
            midbrainOuterArea = false(SIZE(1), SIZE(2));
            forebrainInnerArea = true(SIZE(1), SIZE(2));
            forebrainOuterArea = false(SIZE(1), SIZE(2));
            
            for iRoiFile = 1:nRoiFiles
                %affine transformation tform object
                
                iTform = tformList(iRoiFile);
                %MIDBRAIN
                %transform to registrated fish
                [midbrainXCoordList{iRoiFile}, midbrainYCoordList{iRoiFile}] = transformPointsForward(iTform, ROI{iRoiFile, 1}.mnCoordinates(:, 1), ROI{iRoiFile, 1}.mnCoordinates(:, 2));
                %{
                 Icombined = imread([ExtendedDeptOfFieldPath,'/',stackinfo(k).stackname,'-ExtendedDeptOfField.tif']);
                 Ialigned = imread([RegistratedPath,'/',stackinfo(k).stackname,'-Registration.tif']);
 
                 figure
                 hold on
                 imagesc(Icombined)
                 plot(ROI{k,1}.mnCoordinates(:,1),ROI{k,1}.mnCoordinates(:,2));
 
                 figure
                 hold on
                 imagesc(Ialigned)
                 plot(xc,yc)
                %}
                %generates outer mask and inner mask
                midbrainAreaMask = poly2mask(double(midbrainXCoordList{iRoiFile}), double(midbrainYCoordList{iRoiFile}), SIZE(1), SIZE(2));
                midbrainInnerArea = midbrainAreaMask & midbrainInnerArea;
                midbrainOuterArea = midbrainAreaMask | midbrainOuterArea;
                %figure;imagesc(BWm)
                
                midbrainInnerContour = bwboundaries(midbrainInnerArea);
                midbrainOuterContour = bwboundaries(midbrainOuterArea);
                
                %stores individual information
                MidbrainInfo(iRoiFile).Poly = ROI{iRoiFile, 1}.mnCoordinates;
                rpm(iRoiFile) = regionprops(midbrainAreaMask, 'Area', 'Centroid', 'BoundingBox');
                
                %FORBRAIN
                [forebrainXCoordList{iRoiFile}, forebrainYCoordList{iRoiFile}] = transformPointsForward(iTform, ROI{iRoiFile, 2}.mnCoordinates(:, 1), ROI{iRoiFile, 2}.mnCoordinates(:, 2));
                %generates outer mask and inner mask
                BWf = poly2mask(double(forebrainXCoordList{iRoiFile}), double(forebrainYCoordList{iRoiFile}), SIZE(1), SIZE(2));
                forebrainInnerArea = BWf & forebrainInnerArea;
                forebrainOuterArea = BWf | forebrainOuterArea;
                %stores individual information
                ForbrainInfo(iRoiFile).Poly = ROI{iRoiFile, 2}.mnCoordinates;
                rpf(iRoiFile) = regionprops(BWf, 'Area', 'Centroid', 'BoundingBox');
            end
            
            %{
              figure;
              imagesc(objt.Template)
              for iRoiFile = 1:nRoiFiles
                  hold on
                  plot(midbrainXCoordList{iRoiFile},midbrainYCoordList{iRoiFile})
              end
 
              figure;
              imagesc(objt.Template)
              for iRoiFile = 1:nRoiFiles
                  hold on
                  plot(forebrainXCoordList{iRoiFile},forebrainYCoordList{iRoiFile})
              end
            %}
            
            
            temp = {rpm.Area}; [MidbrainInfo.Area] = temp{:};
            temp = {rpm.Centroid}; [MidbrainInfo.Centroid] = temp{:};
            temp = {rpm.BoundingBox}; [MidbrainInfo.BoundingBox] = temp{:};
            
            temp = {rpf.Area}; [ForbrainInfo.Area] = temp{:};
            temp = {rpf.Centroid}; [ForbrainInfo.Centroid] = temp{:};
            temp = {rpf.BoundingBox}; [ForbrainInfo.BoundingBox] = temp{:};
            
            
            %%
            %Midbrain
            bandm = boolean(midbrainOuterArea-midbrainInnerArea);
            distancemapm = bwdist(~bandm);
            distancemapblurm = imgaussfilt(distancemapm, 20);
            skeletonm = bwmorph(bandm, 'shrink', Inf);
            [X, Y] = find(skeletonm);
            MeanMidBrain = sng_OrderContourCoordinates([X, Y]);
            
            
            % centrum based on the region in the center of brain with the largest
            % distance from all brain edges
            dd = bwdist(~midbrainInnerArea);
            [~, cx] = max(max(dd));
            [~, cy] = max(max(dd'));
            Centroidm = [cx, cy];
            
            %Forbrain
            bandf = boolean(forebrainOuterArea-forebrainInnerArea);
            distancemapf = bwdist(~bandf);
            distancemapblurf = imgaussfilt(distancemapf, 20);
            skeletonf = bwmorph(bandf, 'shrink', Inf);
            [X, Y] = find(skeletonf);
            MeanForBrain = sng_OrderContourCoordinates([X, Y]);
            
            %it makes not sence to forbrain center based on all fishes as there is no
            %center which is enclosed by all forbrains
            
            
            %{
         figure;imagesc(C1m)
         figure;imagesc(C2m)
         figure;imagesc(bandm);axis off equal tight;
         figure;imagesc(distancemapm)
         figure;imagesc(distancemapblurm)
         figure;imagesc(skeletonm)
 
         %to show skeleton more clear
         skeletonblur = imgaussfilt(double(skeletonm),2);
         figure;imagesc(skeletonblur);axis off equal tight;
 
         figure;imagesc(C1f)
         figure;imagesc(C2f)
         figure;imagesc(bandf);axis off equal tight;
         figure;imagesc(distancemapf)
         figure;imagesc(distancemapblurf)
         figure;imagesc(skeletonf)
            %}
            
            %             if saveTF
            %                 save([TemplatePath3dpf, '/MeanMidbrain', '.mat'], 'bandm', 'distancemapm', ...
            %                     'distancemapblurm', 'skeletonm', 'MeanMidBrain', 'Centroidm')
            %                 save([TemplatePath3dpf, '/MeanForbrain', '.mat'], 'bandf', 'distancemapf', ...
            %                     'distancemapblurf', 'skeletonf', 'MeanForBrain')
            %                 save([TemplatePath3dpf, '/MidbrainInfo', '.mat'], 'MidbrainInfo')
            %                 save([TemplatePath3dpf, '/ForbrainInfo', '.mat'], 'ForbrainInfo')
            %             end
            
            objt.MeanMidBrain = MeanMidBrain;
            objt.CenterMidBrain = Centroidm;
            objt.MidBrainDistanceMap = distancemapm;
            objt.BandMidBrain = bandm;
            
            objt.midbrainXCoordList = midbrainXCoordList;
            objt.midbrainYCoordList = midbrainYCoordList;
            objt.forebrainXCoordList = forebrainXCoordList;
            objt.forebrainYCoordList = forebrainYCoordList;
            
            objt.distancemapblurm = distancemapblurm;
            objt.skeletonm = skeletonm;
            objt.MeanForeBrain = MeanForBrain;
            objt.ForeBrainDistanceMap = distancemapf;
            objt.BandForeBrain = bandf;
            objt.distancemapblurf = distancemapblurf;
            objt.skeletonf = skeletonf;
            
            objt.MidbrainInfo = MidbrainInfo;
            objt.ForbrainInfo = ForbrainInfo;
            
            objt.midbrainInnerContour = midbrainInnerContour{1};
            objt.midbrainOuterContour = midbrainOuterContour{1};
            
        end
        function polarProbabilityMidbrain(objt)
            %% Polar transform distancemap
            
            %DONE: use more suffisticated method to exclude area in shortest path
            %finding, multiply shortest path with based on brains blurred distancemap
            %to exclude some area from pathfinding
            
            %TODO: add band to generateskeleton function output instead of computing from distancemap
            
            %using the distancemap wont help as the probability goes to much down
            %for areas where the real probability is high.
            %IsquareMB = sng_boxaroundcenter(objt.MidBrainDistanceMap,objt.CenterMidBrain,[],0);
            
            %1000x1000 box around midbrain band
            squareMidbrainBand = sng_boxaroundcenter(objt.BandMidBrain, objt.CenterMidBrain, [], 0);
            
            
            %{
              figure;imagesc(objt.BandMidBrain)
              figure;imagesc(IsquareMB)
              figure;imagesc(boolean(IsquareMB));sng_imfix;colormap gray
            %}
            
            polarMidbrainBand = sng_Im2Polar3(squareMidbrainBand); %polar transform
            polarMidbrainBand = permute(polarMidbrainBand, [2, 1, 3]); %horizontal orientation
            %{
              figure;imagesc(Polarmidbrainband)
            %}
            
            %IDistmap3 contains added gausian filter at te boudary i.e. the values of
            %IDistmap which are 1 stays 1, i.e. the plane is preserved
            %sigmabasedVariable = 6; added to propeties
            
            [nRows, nColumns] = size(polarMidbrainBand);
            polarMidbrainBandWithGaussian = polarMidbrainBand;
            
            center = 501;
            listGaussians = cell(1, nColumns);
            
            for iColumn = 1:nColumns
                a = find(polarMidbrainBand(:, iColumn) == 1, 1, 'first');
                b = find(polarMidbrainBand(:, iColumn) == 1, 1, 'last');
                
                %a gaussian formula
                %y = gaussmf([-500:1:500],[(b-a)/blurvar 0])';
                %for versions without fuzzy logic toolbox
                listGaussians{iColumn} = exp(-((-500:1:500) - 0).^2/(2 * ((b - a) / objt.sigmabasedVariable).^2))';
                polarMidbrainBandWithGaussian(1:a, iColumn) = listGaussians{iColumn}(center - a + 1:center);
                polarMidbrainBandWithGaussian(b:nRows, iColumn) = listGaussians{iColumn}(center:center + nRows - b);
            end
            %a bit of gaussian blur added
            
            %{
                  figure;imagesc(polarMidbrainBandWithGaussian)
            %}
            
            %extraSmoothingKernel = [5,5] %added to properties
            polarMidbrainBandWithGaussianWithBlurr = imgaussfilt(polarMidbrainBandWithGaussian, objt.extraSmoothingKernel);
            %figure;plot(y) % a gaussian
            
            %{
                  figure;imagesc(IDistmap4)
            %}
            
            objt.squareMidbrainBand = squareMidbrainBand;
            objt.polarMidbrainBand = polarMidbrainBand;
            objt.polarMidbrainBandWithGaussian = polarMidbrainBandWithGaussian;
            objt.polarMidbrainBandWithGaussianWithBlurr = polarMidbrainBandWithGaussianWithBlurr;
            objt.listGaussians = listGaussians;
            
        end
        function BrainEdgeTemplate(objt)
            %% generate useable edge filter from edge template
            %add more templates later on
            %filterwidth = 21; %has to be an odd numberl
            
            EdgeTemp = objt.EdgeTemplate1;
            EdgeTempMean = mean(EdgeTemp, 2);
            objt.EdgeFilter = repmat(EdgeTempMean, 1, objt.filterwidth);
            
            
            %{
               figure;
               imagesc(uint8(EdgeTemp));axis off tight equal
               set(gca,'position',[0 0 1 1],'units','normalized')
               truesize(gcf,12*sng_size(get(get(gca,'Children'),'Cdata'),[1,2]))
 
               figure;imagesc(uint8(EdgeTemp));axis off equal tight
               figure;imagesc(uint8(EdgeTempMean));sng_imfix
 
               figure;
               imagesc(uint8(EdgeFilter));axis off tight equal
               set(gca,'position',[0 0 1 1],'units','normalized')
               truesize(gcf,12*sng_size(get(get(gca,'Children'),'Cdata'),[1,2]))
            %}
            
        end
        
        function TemplateSpot(objt)
            %this function computes the template variables of the test batch with annotations.
            
            
            %TODO split deze functie verder op
            %save alleen het masker per spot want dat duurd het langst
            %bereken eroded spotmasker en background masker apart
            %bereken colors later
            
            
            nFishes = 1:numel(objt.SpotNGliaObject.StackInfo);
            %nFishes = 1;
            
            CubeSizeSpotErode = objt.CubeSizeSpotErode; %spot structure element for erode
            CubeSizeSpotDilate = objt.CubeSizeSpotDilate; %background structure element for dilate
            winsiz = objt.winsiz; %size of spot box is 2*winsiz + 1 = 51x51
            
            f = int16(-winsiz:winsiz);
            
            subplot = @(m, n, p) subtightplot(m, n, p, [0, 0], [0, 0], [0, 0]);
            image_TF = false;
            Save_TF = false;
            HistEq_TF = false; %Wavelet spotdetection doesn work very well when used histogram equalization
            
            %ColorToGrayVector = [-0.378629920630393;0.865619921584407;-0.327630179561694]
            
            spotcolorsTT = [];
            backrcolorsTT = [];
            
            
            for iFish = nFishes
                disp(num2str(iFish))
                %transformation matrix to aligned fish
                %tform_1234 = RegistrationInfo{k1}(strcmp({RegistrationInfo{k1}.name},'tform_complete')).value;
                tform_1234 = objt.SpotNGliaObject.RegObject(iFish).tform_1234;
                
                %open correctedfish
                %CorrectedSlice = sng_openimstack2([obj.SavePath, '/', 'CorrectedFish', '/', obj.StackInfo(k1).stackname, '.tif']);
                CorrectedSlice = objt.SpotNGliaObject.PreprocessingObject(iFish).correctedSlice;
                
                %annotated spots
                RoiMicroglia = ReadImageJROI([objt.SpotNGliaObject.AnnotatedSpotPath, '/', objt.SpotNGliaObject.StackInfo(iFish).stackname, '.roi']);
                spota = RoiMicroglia.mfCoordinates;
                spots = RoiMicroglia.vnSlices;
                
                %%it seems that this does nothing
                %[SpotAnn(:,1),SpotAnn(:,2)] = transformPointsForward(tform_1234,spota(:,1),spota(:,2));
                
                
                %for k2 = 1:4
                % mean(CorrectedSlice{k2}(CorrectedSlice{k2}(:) < bkt(1)))
                %end
                
                
                %histogram equalization
                %Wavelet spotdetection doesn work very well when used histogram equalization
                if HistEq_TF
                    bkt = RegistrationInfo{iFish}(strcmp({RegistrationInfo{iFish}.name}, 'BackgroundThreshold')).value;
                    %sl = RegistrationInfo{k1}(strcmp({RegistrationInfo{k1}.name},'stretchlim_0procent')).value;
                    %sl(1,1:3) = 0;
                    for k2 = 1:numel(CorrectedSlice)
                        %histogram stretching
                        %IStrech{k2} = imadjust(CorrectedSlice{1},sl,[]);
                        
                        %histogram equalization excluding background
                        bp = sum(CorrectedSlice{k2}(:) > bkt(1)) / numel(CorrectedSlice{k2}); %fraction of background pixels
                        hgram = [(1 - bp) / 254 * ones(1, 254), bp]; %flat histogram unless large peak containing background
                        IEqual{k2} = histeq(CorrectedSlice{1}, hgram);
                        %{
                  figure;imagesc(CorrectedSlice{k2}(:,:,2))
                  figure;imagesc(IStrech{k2}(:,:,2))
                  figure;imagesc(IEqual{k2}(:,:,2))
                        %}
                    end
                else
                    IEqual = CorrectedSlice;
                end
                
                
                % % frame around annotated spots
                
                backrcolorsT = [];
                spotcolorsT = [];
                
                
                if image_TF
                    subplotvar = ceil(sqrt(2*numel(spots)));
                    figure;
                end
                
                
                for l = 1:numel(spots)
                    Spot1 = IEqual{spots(l)}(f + spota(l, 2), f + spota(l, 1), 1:3);
                    Spot2 = imgaussfilt(Spot1, 0.7);
                    GSpot2 = sng_RGB2Gray(Spot2, objt.ColorToGrayVector, false); %gray transformed spot
                    Mask2 = sng_SpotWavelet(GSpot2, objt.ScaleLevels, objt.ScaleBase, objt.MPlevels, objt.MPthreshold, false);
                    SE1 = strel('cube', objt.CubeSizeSpotErode); %spot structure element for erode
                    SE2 = strel('cube', objt.CubeSizeSpotDilate); %background structure element for dilate
                    SpotMask = imerode(Mask2, SE1);
                    BackMask = ~imdilate(Mask2, SE2);
                    
                    SpotF = Spot1; SpotF(repmat(SpotMask, 1, 1, 3) == 0) = 0;
                    BackF = Spot1; BackF(repmat(BackMask, 1, 1, 3) == 0) = 0;
                    
                    %{
          figure;imagesc(Spot1)
          figure;imagesc(Spot2)
          figure;imagesc(GSpot2)
          figure;imagesc(Mask2)
          figure;imagesc(SpotMask)
          figure;imagesc(BackMask)
          figure;imagesc(SpotF)
          figure;imagesc(BackF)
                    %}
                    
                    spotcolorsarray = Spot1(repmat(SpotMask, 1, 1, 3));
                    spotcolors = reshape(spotcolorsarray, numel(spotcolorsarray)/3, 3);
                    
                    backrcolorsarray = Spot1(repmat(BackMask, 1, 1, 3));
                    backrcolors = reshape(backrcolorsarray, numel(backrcolorsarray)/3, 3);
                    
                    %used for previous color measure probability
                    spotcolorsT = [spotcolorsT; spotcolors];
                    backrcolorsT = [backrcolorsT; backrcolors];
                    
                    
                    if image_TF
                        subplot(subplotvar, subplotvar, (2 * l)-1); imagesc(SpotF); axis off tight equal
                        subplot(subplotvar, subplotvar, (2 * l)); imagesc(BackF); axis off tight equal
                    end
                end
                
                drawnow
                spotcolorsTT = [spotcolorsTT; spotcolorsT];
                backrcolorsTT = [backrcolorsTT; backrcolorsT];
                
                objt.spotcolorsTT = spotcolorsTT;
                objt.backrcolorsTT = backrcolorsTT;
                
                objt.Spot1 = Spot1;
                objt.SpotF = SpotF;
                objt.BackF = BackF;
                objt.Mask2 = Mask2;
                
                
                %     %
                %{
                 spotcolorsTT2 = unique(objt.spotcolorsTT, 'rows');
                 backrcolorsTT2 = unique(objt.backrcolorsTT, 'rows');
 
                 figure; scatter3(spotcolors(:, 1), spotcolors(:, 2), spotcolors(:, 3));
                 hold on;scatter3(backrcolors(:,1),backrcolors(:,2),backrcolors(:,3));
                 axis equal;xlim([1,255]);ylim([1,255]);zlim([1,255]);xlabel('Red');ylabel('Green');zlabel('Blue')
 
                 figure; scatter3(spotcolorsTT2(:, 1), spotcolorsTT2(:, 2), spotcolorsTT2(:, 3));
                 hold on;scatter3(backrcolorsTT2(:,1),backrcolorsTT2(:,2),backrcolorsTT2(:,3));
                 axis equal;xlim([1,255]);ylim([1,255]);zlim([1,255]);xlabel('Red');ylabel('Green');zlabel('Blue')
                 drawnow
                %}
                
            end
            
            
        end
        function TemplateSpotpart2(objt)
            S = double(objt.spotcolorsTT);
            B = double(objt.backrcolorsTT);
            
            %
            
            % % feature segmentation.
            %{
  data = double([backrcolorsTT2;spotcolorsTT2]);
  label = [(ones(size(backrcolorsTT2,1),1));2*(ones(size(spotcolorsTT2,1),1))];
  A = prdataset(data(:,1:3),label)
 
  W2 =  perlc(A)
 
  figure;scatterd(A(:,1:2))
  plotc({W2})
 
  mat = eye(3)
  mat(:,1:2) = struct(W2).data.rot
 
  D = A(:,1:2)*W2
 
  figure;scatterd(D)
            %}
            % % histogram spot/background contrast.
            %{
  bin = 255;
 
 
  Shist=sng_chistcount(S',linspace(0,255,bin+1));      %create bins
  Bhist=sng_chistcount(B',linspace(0,255,bin+1));      %
 
  Snorm = Shist/size(S,1);
  Bnorm = Bhist/size(B,1);
 
  [minC1,C1,Q1] = sng_threshold2(Snorm,Bnorm);                %find threshold
 
  sng_chistplot(Snorm,Bnorm,minC1);set(gcf,'numbertitle','off','name','RGB')
 
  %% a new way to transform to gray instead of taking the green values
  bm = mean(B)
  sm = mean(S)
 
 
  b = (bm-sm)
  bnorm = b/norm(b)
 
  mx = dot([255,255,255],bnorm)
  mn = dot([0,0,0],bnorm)
 
  st = double(S)*bnorm'; %transformed spot
  bt = double(B)*bnorm'; %transformed background
 
  st = st* 255/mx;
  bt = bt* 255/mx;
 
  SThist = histcounts(st,linspace(0,255,bin+1));
  BThist = histcounts(bt,linspace(0,255,bin+1));
 
  STnorm = SThist/size(S,1);
  BTnorm = BThist/size(B,1);
 
  [minC2,C2,Q2] = sng_threshold2(STnorm,BTnorm);                %find threshold
 
  sng_chistplot(STnorm,BTnorm,minC2);set(gcf,'numbertitle','off','name','RGB')
 
 
  Overlap = min(STnorm,BTnorm);
  Combination = max(STnorm,BTnorm);
  Jaccard = 1-sum(Overlap)/sum(Combination)
            %}
            
            
            [cartvec, bestJac, polarvec, index] = objt.FindBestContrastVector(S, B, 7, 11, false, true);
            bestvec = cartvec(:, end);
            
            if image_TF
                Jc1 = sng_Rgb2BwContrast(S, B, bestvec, true); set(gca, 'XLim', [75, 160])
                Jc2 = sng_Rgb2BwContrast(S, B, [1; 0; 0], true);
                Jc3 = sng_Rgb2BwContrast(S, B, [0; 1; 0], true);
                Jc4 = sng_Rgb2BwContrast(S, B, [0; 0; 1], true);
            end
            
            
            % % Compute BW image based on transform vector
            k10 = 10
            CorrectedSlice = sng_openimstack2([PreprocessionPath, '/', stackinfo(k10).stackname, '.tif']);
            Icombined = sng_SliceCombine(CorrectedSlice, stackinfo(k10).ExtendedDeptOfField.IndexMatrix);
            
            Img2 = sng_RGB2Gray(Icombined, bestvec, true);
            %Img2 = sng_RGB2Gray(Icombined,[0;1;0],true);
            
            
            % % create volume that contains which a color vector can compared with to determine if it is in the range of spots
            
            %spot colorvectors scatter plot
            %{
  figure;scatter3(spotcolorsTT2(:,1),spotcolorsTT2(:,2),spotcolorsTT2(:,3),'filled',1);
  axis equal;xlim([1,255]);ylim([1,255]);zlim([1,255]);xlabel('Red');ylabel('Green');zlabel('Blue')
 
  %mask spot colorvectors
  SpotMask = false(255,255,255);
  for n = 1:size(spotcolorsTT2,1)
      SpotMask(spotcolorsTT2(n,1),spotcolorsTT2(n,2),spotcolorsTT2(n,3)) = true;
  end
  figure;hi = imagesc(SpotMask(:,:,1))
  for m=1:255
      set(hi,'CData',SpotMask(:,:,m))
      drawnow
  end
            %}
            
            %mask spot colorvectors with closed holes. Doesnt work for to less data, i.e. too
            %much single pixels
            %{
  SpotMask2=imfill(SpotMask,'holes');
  figure;hi = imagesc(SpotMask2(:,:,1))
  for m=1:255
      set(hi,'CData',SpotMask2(:,:,m))
      drawnow
  end
            %}
            
            %counts how much from each color vector, i.e. generates a valued Mask
            %correction from 255 to 256 and added +1
            %the spotdetectionlink function has also to be changed
            SpotColorImg = zeros(256, 256, 256);
            for n2 = 1:size(S, 1)
                SpotColorImg(S(n2, 1)+1, S(n2, 2)+1, S(n2, 3)+1) = SpotColorImg(S(n2, 1)+1, S(n2, 2)+1, S(n2, 3)+1) + 1;
            end
            %figure;hist(SpotImg(SpotImg ~= 0),1:8)
            
            SpotMaskGauss = imgaussfilt3(double(SpotColorImg), 5);
            
            
            SVAP_index = find(SpotMaskGauss ~= 0);
            %[subscript(:,1),subscript(:,2),subscript(:,3)] = (ind2sub([255,255,255],SVAP_index));
            %SVAP_subscript = uint8(subscript);
            SpotVectorArrayProbability = SpotMaskGauss(SVAP_index);
            
            
            %{
  writerObj1 = VideoWriter('blurred spots','MPEG-4');
  writerObj1.FrameRate = 50;    % to perform realtime movie
  open(writerObj1);
  figure;hi = imagesc(SpotMaskGauss(:,:,1),[0 max(SpotColorImg(:))]);colormap(jet(1000));sng_imfix(2)
  for m = 1:255
      set(hi,'CData',SpotMaskGauss(:,:,m))
      drawnow
      frame = getframe(gcf);
      writeVideo(writerObj1,frame)
  end
  close(writerObj1);
            %}
            
            %sum(SpotMask(:))
            %sum(SpotMask2(:))
            
            %Regionstats = regionprops(SpotMask) %to much single pixels
            %Regionstats = regionprops(boolean(SpotMaskGauss)) %to much single pixels
            
            
            %{
  [X,Y,Z] = meshgrid(1:255,1:255,1:255);
  S2 = double(spotcolorsTT2);
  kk = boundary(S2(:,1),S2(:,2),S2(:,3));
  figure;trisurf(kk,S2(:,1),S2(:,2),S2(:,3),'Facecolor','red','FaceAlpha',0.1)
 
 
  DT = delaunayTriangulation(S2(:,1),S2(:,2),S2(:,3));
  K = convexHull(DT)
  figure;trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),'Facecolor','cyan','FaceAlpha',0.1)
            %}
            
            
            SpotTemplateVar.SpotContrastVector = bestvec;
            SpotTemplateVar.SpotVectorArrayProbability = SpotVectorArrayProbability;
            SpotTemplateVar.SVAP_index = SVAP_index;
            
            if Save_TF
                save([Basepath, '/Template Spot/SpotTemplate.mat'], 'SpotTemplateVar');
                save([TemplatePath3dpf, '/SpotTemplate.mat'], 'SpotTemplateVar');
            end
            
            % % strech image the optain lower variance of spot color
            %{
 
 
            %}
            
        end
        
    end
    methods(Static)
        
        function [cartvec, bestJac, polarvec, index] = FindBestContrastVector(array1, array2, iter, res, image_tf, waitbar_tf)
            %find best contrast vector in a couple of iterations between to 2 classes with 3 features
            %
            
            %
            %Example
            %
            %[cartvec bestJac polarvec index]=FindBestContrastVector(array1,array2,iter,res,imageyn)
            
            
            %{
 
 
 iter = 5;
 res = 5
 imageyn_tf = true;
 waitbar_tf = true
            %}
            
            clear vector Jaccard Jacplot;
            
            
            if waitbar_tf;
                hwb = waitbar(0, 'iteration', 'Name', 'Find Best Contrast Vector');
            end
            
            measurements = res^2;
            
            for k = 1:iter
                
                if waitbar_tf;
                    l = 1;
                    waitbar((((k - 1) * measurements) + l)/(measurements * iter), hwb, ['iteration ', num2str(k), ' of ', num2str(iter)]);
                end
                
                
                %create spherical unit vector in equal spaced positive directions
                if k == 1
                    theta = linspace(0, pi, res);
                    phi = linspace(0, pi, res);
                else
                    rg = pi * 0.5^k; %range
                    theta = linspace(direction(1, mr)-rg, direction(1, mr)+rg, res);
                    phi = linspace(direction(2, mr)-rg, direction(2, mr)+rg, res);
                end
                
                direction = combvec(theta, phi);
                
                vector(1, :) = sin(direction(1, :)) .* cos(direction(2, :));
                vector(2, :) = sin(direction(1, :)) .* sin(direction(2, :));
                vector(3, :) = cos(direction(1, :));
                
                
                for l = 1:measurements
                    [Jaccard(l)] = sng_Rgb2BwContrast(array1, array2, vector(:, l), false);
                    
                    if waitbar_tf && floor(rem((((k - 1) * measurements) + l), (measurements * iter / 100))) == 0
                        
                        waitbar((((k - 1) * measurements) + l)/(measurements * iter));
                    end
                    
                    %(((k-1)*measurements)+l)/(measurements*iter)
                end
                
                [mx, mr] = max(Jaccard);
                
                index(k) = mr
                cartvec(1:3, k) = vector(:, mr)
                bestJac(k) = mx
                polarvec(1:2, k) = direction(1:2, mr)
                
                
                if image_tf
                    Jacplot = reshape(Jaccard, size(theta, 2), size(phi, 2));
                    figure; surf(theta, phi, Jacplot)
                    xlabel('theta'); ylabel('phi')
                    xlim(sng_lim(theta))
                    ylim(sng_lim(phi))
                    colormap colorcube(200)
                end
                
            end
            
            if waitbar_tf;
                close(hwb);
            end
        end
        
        
    end
    
end

