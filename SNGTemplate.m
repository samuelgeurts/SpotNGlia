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
        
        SpotNGliaObject
        
        %Template 
        Template
        Size
        ref_temp
        
        %Midbrain 
        CenterMidBrain
        MidBrainDistanceMap
        %Edge Template
        EdgeTemplate1
        
        %Spot 
        SpotVectorArrayProbability
        SVAP_index
    end
    properties(Transient = true)
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
        BandMidBrain
        distancemapblurm
        skeletonm
        midbrainXCoordList
        midbrainYCoordList
        midbrainInnerContour
        midbrainOuterContour
        
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
    end 
    properties(Transient = true)
        % Brain Images
        
 
    end
    
    methods
        function objt = SNGTemplate
            disp('select template source path');
            objt.sourcePath = uigetdir(pwd, 'select template source path');
            %folder with all 3dpf template information
            objt.TemplatePath3dpf = [objt.sourcePath, '/', 'Template 3 dpf'];
            %folder with all brain rois
            objt.RoiBrainPath = [objt.sourcePath, filesep, 'Roi brain']; 
        end
        function value = get.midbrainXCoordList(objt)
            if isempty(objt.midbrainXCoordList)
                disp('compute Brain parameters');
                objt.BrainTemplate; %compute all Brain paramaters including midbrainXCoordList
            end
            value = objt.midbrainXCoordList;
        end
        function loadTemplate(objt)
            % Mean Fish Registration Template
            %based on older template, can be updated with 50 fish batch but lead
            %probably not to improvement. Source is not sure, TemplateGeneration?
            
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
            save(strcat(objt.sourcePath, filesep, objt.fileName, '.mat'), 'objt');
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
            objt.SpotNGliaObject.LoadParameters('RegObject');
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
        
    end
end

