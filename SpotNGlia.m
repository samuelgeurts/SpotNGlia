classdef SpotNGlia

    properties   
        FishPath = []
        SourcePath = []
        SavePath = []
        
        AnnotatedMidBrainPath = []
        AnnotatedSpotPath = []
        
        SaveName = []
        InfoName = []
               
        fishnumbers = []
        slicenumbers = []
        savedate = []        
        ZFParameters = []        

        CompleteTemplate = []
        
        StackInfo = []
        ImageInfo = []
                        
        InfoResult = []
        
        Annotations = []
        
        %CorrectedFish = []
        %CombinedFish = []
        %AlignedFish = []
    end
      
    properties

        BrainInfo = []
        BrainStats = []

        SpotInfo = []
        SpotStats = []
        %SpotSelection = []      
        
        SpotBrainInfo = []
        SpotBrainStats = []    
    end 
   
    properties
       ImageInfoChecked_TF = false
       Sorting = []
       
    end

    
    %properties (Access = private)
    %    %stackinfo
    %end
    
    methods           
        
        %constructor function
        function obj = SpotNGlia(mode1)
       
            if exist('mode1','var') && mode1 == 0 %training folder
                if strcmp(getenv('OS'),'Windows_NT') && strcmp(getenv('username'),'260018')
                    BasePath = 'C:\Users\260018\Dropbox';        
                elseif strcmp(getenv('OS'),'Windows_NT') && strcmp(getenv('username'),'SNGeu')
                    BasePath = 'C:\Users\SNGeu\Dropbox';
                elseif  strcmp(getenv('USER'),'samuelgeurts')
                    BasePath = '/Users/samuelgeurts/Dropbox';
                elseif  strcmp(getenv('USER'),'Anouk')
                    BasePath = '/Users/Anouk/Dropbox';
                else
                    error('unknown platform and user');
                end
                obj.FishPath = [BasePath,'/','20170327_3dpf'];
                obj.SavePath = [BasePath,'/','SpotNGlia Destination'];                           
                obj.SourcePath = [BasePath,'/','SpotNGlia Source'];                       
                obj.AnnotatedMidBrainPath = [BasePath,'/','SpotNGlia Source','/','Roi brain'];    
                obj.AnnotatedSpotPath = [BasePath,'/','SpotNGlia Source','/','Roi microglia'];    
            else
                disp('select image folder to process')
                obj.FishPath = uigetdir([],'select image folder to process');
                disp('select destination folder')
                obj.SavePath = uigetdir([],'select destination folder');
                disp('select program source folder')
                obj.SourcePath = uigetdir([],'select program source folder');                
            end            
            
            [~,folder] = fileparts(obj.FishPath);
            obj.SaveName = strcat('SNG_',folder);       
            obj.InfoName = strcat('INFO_',folder);        
            
            load([obj.SourcePath,'/','zfinput.mat']);   
            obj.ZFParameters = zfinput;
                        
            %file which contains all batch specific computation info
            matObj = matfile([obj.SavePath,'/',obj.InfoName,'.mat']);
            obj.InfoResult = whos(matObj);
        end   
        
        function obj = SliceCombination(obj,slicenumbers)
            %computes imageinfo and stackinfo
            %       sorting = ['date'/'name']
            %       CheckImageInfo_TF [true/false]
            %Example:   obj = SliceCombination(obj)
            %Example:   obj = SliceCombination(obj,'date',true)
            %Example:   obj = SliceCombination(obj,'name',true)

            if ~isempty(obj.ImageInfo)
                answer = questdlg('Are you sure to overwrite ImageInfo and StackInfo','','Ok','Cancel','Cancel');
                if strcmp(answer,'Ok')
                    obj.ImageInfo = [];
                    obj.StackInfo = [];
                    obj.ImageInfoChecked_TF = false;
                    obj.Sorting = [];
                else
                    return
                end
            end
            
            if isempty(obj.Sorting)                            
                % set imagesorting to date or name if not given               
                obj.Sorting = questdlg('Which sorting algorithm do you want?','Choose sorting method','Date','Name','Date');
                if strcmp(obj.Sorting,'Date')
                    obj.ZFParameters = sng_zfinput(obj.ZFParameters,0,'imageinfo','stackselection','sorting','Date','high'); %date or name  
                elseif strcmp(obj.Sorting,'Name')
                    obj.ZFParameters = sng_zfinput(obj.ZFParameters,0,'imageinfo','stackselection','sorting','Name','high'); %date or name                  
                end
            end
                
            if isempty(obj.ImageInfo)                    
                % compute ImageInfo and StackInfo
                dirinfo =  dir([obj.FishPath,'/*.','tif']);
                imageinfotemp = rmfield(dirinfo,{'isdir','datenum','bytes'}); %removes unimportant fields
                               
                if ~isfield(imageinfotemp,'folder')
                    [imageinfotemp.folder] = deal(obj.FishPath)
                end                

                % make en selection of fishslices if fisnumbers is given as input
                if ~exist('slicenumbers','var')
                    obj.ImageInfo = imageinfotemp;
                    slicenumbers = 1:numel(imageinfotemp);
                elseif (max(slicenumbers) <= numel(imageinfotemp))
                    obj.ImageInfo = imageinfotemp(slicenumbers);    
                else 
                    error('at least one fish does not exist in StackInfo')
                end
                obj.slicenumbers.slicecombination = slicenumbers;
                    
                [obj.ImageInfo] = ImageInfoLink3(obj.ImageInfo,obj.ZFParameters);
                [obj.StackInfo] = StackInfoLink2(obj.ImageInfo);
                
                ImageInfo = obj.ImageInfo; %#ok<NASGU,PROPLC>
                StackInfo = obj.StackInfo;  %#ok<NASGU,PROPLC>
                
                save([obj.SavePath,'/',obj.InfoName,'.mat'],'ImageInfo')
                save([obj.SavePath,'/',obj.InfoName,'.mat'],'StackInfo','-append')
            end
            
            if ~obj.ImageInfoChecked_TF
                
                %openvar('ImageInfo')
                %openvar('StackInfo')                
                
                str = sprintf(['Check slice combination in ImageInfo and StackInfo.',...
                    '\nApply corrections in "imageinfo.CorNextStack".',...
                    '\n   Set value "1" for new fish',...
                    '\n   Set value "2" for removing image from stack.',...
                    '\n   Set value "0" if slice belongs to previous slice']);               
                msgbox(str);
            end
            obj.saveit
            
            
        end   
        
        function obj = PreProcession(obj,fishnumbers)
            
            h = waitbar(0,'Preprocession','Name','SpotNGlia');  
         
            if ~isempty(obj.ImageInfo)
                [obj.StackInfo] = StackInfoLink2(obj.ImageInfo);
            else
                error('First run SliceCombination');
            end
            
            if ~obj.ImageInfoChecked_TF
                answer = questdlg('Do you confirm ImageInfo?','','Yes','No','No');
                 if strcmp(answer,'No')
                    msgbox('First run SliceCombination');
                 elseif strcmp(answer,'Yes')
                    obj.ImageInfoChecked_TF = true;
                 end  
            end
            
            if obj.ImageInfoChecked_TF
                obj.slicenumbers.stackinfo = find([obj.ImageInfo.CorNextStack] ~= 2);
                
                if ~exist('fishnumbers','var')
                    fishnumbers = 1:numel(obj.StackInfo);
                elseif (max(fishnumbers) > numel(obj.StackInfo))
                    error('at least one fish does not exist in StackInfo')
                end
                nfishes = numel(fishnumbers);                    
                obj.fishnumbers.preprocession = fishnumbers;   

                %only preallocate if more than 1 fish has to be processed
                if nfishes > 1
                    PreprocessionInfo(nfishes) = struct('ColorWarp',[],...
                        'SliceWarp',[],...
                        'NNC_Img_before',[],...
                        'NNC_Img_after',[],...
                        'NNC_BP_before',[],...
                        'NNC_BP_after',[]);
                end
                
                %folder to save CorrectedFish
                TempFolderName = ([obj.SavePath,'/','CorrectedFish']);
                if ~exist(TempFolderName,'dir')
                    mkdir(TempFolderName)
                end
                
                for k1 = 1:nfishes                
                    
                    waitbar(k1/nfishes,h,'Preprocession')                
                    
                    %select slices
                    fn = fishnumbers(k1);
                                        
                    ImageSlice = cell(1,obj.StackInfo(fn).stacksize);%preallocate for every new slice
                    for k2 = 1:obj.StackInfo(fn).stacksize
                        ImageSlice{k2} = imread([obj.FishPath,'/',obj.StackInfo(fn).imagenames{k2}]);
                    end
                    
%                    if savefig1_TF
%                        sng_SaveCell2TiffStack(ImageSlice,[SavePath,'/stack/',stackinfo(k1).stackname,'.tif'])
%                    end

                    %[CorrectedSlice,ppoutput(k1)] = PreprocessionLink(ImageSlice,obj.ZFParameters);
                    if k1 == 1
                        [CorrectedFishTemp,PreprocessionInfo] = PreprocessionLink(ImageSlice,obj.ZFParameters);               
                    else
                        [CorrectedFishTemp,PreprocessionInfo(k1)] = PreprocessionLink(ImageSlice,obj.ZFParameters);
                    end

                    sng_SaveCell2TiffStack(CorrectedFishTemp,[TempFolderName,'/',obj.StackInfo(k1).stackname,'.tif'])
                                       
                    %obj.CorrectedFish(k1).image = CorrectedFishTemp; <- %save to object takes to much space
        
%TODO for imaging   [CorrectedSlice,ppoutput,ImageSliceCor,FiltIm] = PreprocessionLink(ImageSlice,zfinput)
%                    if savefig2_TF
%                        sng_SaveCell2TiffStack(CorrectedSlice,[obj.SavePath,'/',StackInfo(k1).stackname,'.tif'])
%                    end
                end               
            end       
            obj.saveit
            delete(h)           
            save([obj.SavePath,'/',obj.InfoName,'.mat'],'PreprocessionInfo','-append')  
        end   
        
        function obj = ExtendedDeptOfField(obj,fishnumbers)
            
            h = waitbar(0,'Extended Dept of Field','Name','SpotNGlia');  
            
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif max(fishnumbers) > numel(obj.StackInfo)
                error('at least one fish does not exist in PreProcessionInfo')                   
            end
            nfishes = numel(fishnumbers);             

            obj.fishnumbers.extendeddeptoffield = fishnumbers; 
            %obj.fishnumbers.extendeddeptoffield = obj.fishnumbers.preprocession(fishnumbers);
            %fix later fishnumber
             
            %{
             %Stack has to be added but than input parameters are needed,
             so the function PreProcessionLink has to be separated in that
             case. For now we choose to store every image the object. If it
             leads to memory issue on other obtion has to be made. For
             example save the images in a folder outsite the object.
             
             for k1 = 1:nfishes                    
                 fn = fishnumbers(k1);
                 ImageSlice = cell(1,obj.StackInfo(fn).stacksize);%preallocate for every new slice
                 for k2 = 1:numel(ImageSlice)
                 ImageSlice = imread([obj.FishPath,'/',obj.StackInfo(fn).imagenames{k2}]);
                 [ImageSliceCor{k2},~] = sng_RGB_IATwarp2(ImageSlice(:,:,1:3),obj.PreprocessionInfo(k1).ColorWarp{1});
                 %[cellImg1{k2}, ~] = iat_inverse_warping(ImageSliceCor{k2}, obj.PreprocessionInfo(k1).SliceWarp{k2}, par.transform, 1:N, 1:M);
                 end
             end
             %}                
            if nfishes > 1
                ExtendedDeptOfFieldInfo(nfishes) = struct('IndexMatrix',[],...
                    'variance_sq',[]);
            end
            
            %folder to save CombinedFish
            TempFolderName = ([obj.SavePath,'/','CombinedFish']);
            if ~exist(TempFolderName,'dir')
                mkdir(TempFolderName)
            end
                
            for k1 = 1:nfishes
                waitbar(k1/nfishes,h)                
                fn = fishnumbers(k1);
                
                CorrectedFish = sng_openimstack2([obj.SavePath,'/','CorrectedFish','/',obj.StackInfo(fn).stackname,'.tif']);             
                [CombinedFishTemp,ExtendedDeptOfFieldInfo(k1)] = ExtendedDeptofFieldLink2(CorrectedFish,obj.ZFParameters);                
                %obj.CombinedFish(k1).image = CombinedFishTemp; because it takes to much space               
                imwrite(uint8(CombinedFishTemp),[TempFolderName,'/',obj.StackInfo(k1).stackname,'.tif'],...
                    'WriteMode','overwrite','Compression','none');
            end
            
            obj.saveit
            save([obj.SavePath,'/',obj.InfoName,'.mat'],'ExtendedDeptOfFieldInfo','-append')  
            delete(h)
        end   
        
        function obj = Registration(obj,fishnumbers)
            
            h = waitbar(0,'Registration','Name','SpotNGlia');  
            
            if isempty(obj.CompleteTemplate)
                obj.CompleteTemplate = LoadTemplateLink3([obj.SourcePath,'/','Template 3 dpf']);                         
            end
                        
            if ~exist('fishnumbers','var')
                fishnumbers = numel(obj.StackInfo);
            elseif max(fishnumbers) > numel(obj.StackInfo)
                error('at least one fish does not exist in ExtendedDeptOfFieldInfo')                   
            end
            nfishes = numel(fishnumbers);             
            obj.fishnumbers.registration = fishnumbers;
            
            %folder to save AlignedFish
            TempFolderName = ([obj.SavePath,'/','AlignedFish']);
            if ~exist(TempFolderName,'dir')
                mkdir(TempFolderName)
            end
                                        
            RegistrationInfo = cell(nfishes,1); 
            for k1 = 1:nfishes
                waitbar(k1/nfishes,h)                                
                fn = fishnumbers(k1);                
                CombinedFish = imread([obj.SavePath,'/','CombinedFish','/',obj.StackInfo(fn).stackname,'.tif']);          
                [AlignedFishTemp,RegistrationInfo{k1,1}] = AllignmentLink5(CombinedFish,obj.CompleteTemplate,obj.ZFParameters);
                %obj.AlignedFish(k1).image = AlignedFishTemp; %because it
                %takes to much space                                
                imwrite(uint8(AlignedFishTemp),[TempFolderName,'/',obj.StackInfo(k1).stackname,'.tif'],...
                    'WriteMode','overwrite','Compression','none');
            end

            obj.saveit                        
            save([obj.SavePath,'/',obj.InfoName,'.mat'],'RegistrationInfo','-append')  
            delete(h)            
        end   
        
        function obj = BrainSegmentation(obj,fishnumbers)
            h = waitbar(0,'Brain Segmentation','Name','SpotNGlia');  
            
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif max(fishnumbers) > numel(obj.StackInfo)
                error('at least one fish does not exist in RegistrationInfo')                   
            end
            
            nfishes = numel(fishnumbers);             
            obj.fishnumbers.brainsegmentation = fishnumbers; 

            if nfishes > 1
            	BrainSegmentationInfo(nfishes) = struct('EdgeFilterWidth',[],...
                    'ShortestPath',[],...
                    'ShortestPathValue',[],...
                    'BrainEdge',[]);
            end            
            
            for k1 = 1:nfishes                    
                fn = fishnumbers(k1);
                waitbar(k1/nfishes,h)
                AlignedFish = imread([obj.SavePath,'/','AlignedFish','/',obj.StackInfo(fn).stackname,'.tif']);                         
                [~,BrainSegmentationInfo(k1)] = MidBrainDetectionLink3(AlignedFish,obj.CompleteTemplate,obj.ZFParameters);      
                %obj.BrainFish(k1).image = BrainFishTemp;
            end
            obj.saveit            
            save([obj.SavePath,'/',obj.InfoName,'.mat'],'BrainSegmentationInfo','-append')  

            delete(h)            
        end            

        function obj = SpotDetection(obj,fishnumbers)
            
            load([obj.SavePath,'/',obj.InfoName,'.mat'],'BrainSegmentationInfo')
        
            h = waitbar(0,'SpotDetection','Name','SpotNGlia');  
            
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif max(fishnumbers) > numel(obj.StackInfo)
                error('at least one fish does not exist in BrainInfo')                   
            end
            
            nfishes = numel(fishnumbers);             
            obj.fishnumbers.spotdetection = fishnumbers; 

            %if nfishes > 1
            %	spoutput(nfishes) = struct('EdgeFilterWidth',[],...
            %        'ShortestPath',[],...
            %        'ShortestPathValue',[],...
            %        'BrainEdge',[]);
            %end
            SpotsDetected = cell(nfishes,1);
            SpotParameters = cell(nfishes,1);
           
            
            if nfishes > 1
            	SpotDetectionInfo(nfishes) = struct('Multiproduct',[],...
                    'MultiproductThreshold',[]);
            end
            
            for k1 = 1:nfishes                    
                fn = fishnumbers(k1);
                waitbar(k1/nfishes,h)
                AlignedFish = imread([obj.SavePath,'/','AlignedFish','/',obj.StackInfo(fn).stackname,'.tif']);                         
                %[~,SpotDetectionInfo{k1,1}] = SpotDetectionLink2(AlignedFish,obj.CompleteTemplate,BrainSegmentationInfo(k1).BrainEdge,obj.ZFParameters);                   
                [SpotsDetected{k1},SpotParameters{k1},SpotDetectionInfo(k1)] = SpotDetectionSNG(AlignedFish,obj.CompleteTemplate,BrainSegmentationInfo(k1).BrainEdge,obj.ZFParameters);
            end
            obj.saveit
            save([obj.SavePath,'/',obj.InfoName,'.mat'],'SpotDetectionInfo','-append')  
            save([obj.SavePath,'/',obj.InfoName,'.mat'],'SpotsDetected','-append')  
            save([obj.SavePath,'/',obj.InfoName,'.mat'],'SpotParameters','-append')  

            delete(h)            
        end            
            
        function obj = CompleteProgram(obj,fishnumbers)
            
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif (max(fishnumbers) > numel(obj.StackInfo))
                error('at least one fish does not exist in StackInfo')
            end        

            obj = obj.PreProcession(fishnumbers);
            obj = obj.ExtendedDeptOfField(fishnumbers);
            obj = obj.Registration(fishnumbers);
            obj = obj.BrainSegmentation(fishnumbers);
            obj = obj.SpotDetection(fishnumbers);                        
        end
                    
        function obj = LoadAnnotations(obj)

            if ~isempty('obj.AnnotatedMidBrainPath') || ~isempty('obj.AnnotatedSpotPath')
                nfishes = numel(obj.StackInfo);            
                tform_1234 = cell(nfishes,1);
                %annotated brain and spots
                load([obj.SavePath,'/',obj.InfoName,'.mat'],'RegistrationInfo')
                for k1 = 1:nfishes
                    tform_1234{k1} = RegistrationInfo{k1}(strcmp({RegistrationInfo{k1}.name},'tform_complete')).value; %#ok<USENS>
                end
            end
            
            %annotated midbrain roi            
            if isempty('obj.AnnotatedMidBrainPath')
                disp('Select Annotated MidBrain Path')
                obj.AnnotatedMidBrainPath = uigetdir([],'Select Annotated MidBrain Path');
            end
            
            if ~isempty('obj.AnnotatedMidBrainPath')         
                ambr = cell(nfishes,1);               
                for k1 = 1:nfishes
                    RoiBrain = ReadImageJROI([obj.AnnotatedMidBrainPath,'/',obj.StackInfo(k1).stackname,'.zip']);

                    amb = sng_roicell2poly(RoiBrain,1);   
                    [ambr{k1}(:,1),ambr{k1}(:,2)] = transformPointsForward(tform_1234{k1},amb(:,1),amb(:,2));
                    ambr{k1} = double(ambr{k1});           
                end                                    
                [obj.Annotations(1:nfishes).MidBrain] = ambr{:};    
            end
            
            %annotated spots            
            if isempty('obj.AnnotatedSpotPath')
                disp('Select Annotated Spot Path')
                obj.AnnotatedSpotPath = uigetdir([],'Select Annotated Spot Path');
            end
            
            if ~isempty('obj.AnnotatedSpotPath')
                     SpotAnn = cell(nfishes,1);
                for k1 = 1:nfishes
                    RoiMicroglia = ReadImageJROI([obj.AnnotatedSpotPath,'/',obj.StackInfo(k1).stackname,'.roi']);

                    SpotAnn{k1} = zeros(size(RoiMicroglia.mfCoordinates,1),2);
                    [SpotAnn{k1}(:,1),SpotAnn{k1}(:,2)] = transformPointsForward(tform_1234{k1},...
                        RoiMicroglia.mfCoordinates(:,1),RoiMicroglia.mfCoordinates(:,2));                
                end
                [obj.Annotations(1:nfishes).Spots] = SpotAnn{:};
            end
                       
            %load([obj.SavePath,'/',obj.InfoName,'.mat'],'BrainSegmentationInfo')
            %load([obj.SavePath,'/',obj.InfoName,'.mat'],'SpotDetectionInfo')
            %
            %computed brain and spots
            %for k1 = 1:numel(obj.StackInfo)
            %    cmbr{k1} = BrainSegmentationInfo(k1).BrainEdge;
            %    obj.SpotParameters{k1} = SpotDetectionInfo{k1}(strcmp({SpotDetectionInfo{k1}.name},'SpotParameters')).value;
            %end
        end
 
        function obj = BrainVal(obj)
            
            load([obj.SavePath,'/',obj.InfoName,'.mat'],'BrainSegmentationInfo')
                     
            nfishes = numel(obj.StackInfo);
            
            Jaccard = zeros(nfishes,1);
            Dice = zeros(nfishes,1);
            
            cmbr = {BrainSegmentationInfo.BrainEdge};
            ambr = {obj.Annotations.MidBrain};      
                       
            for k1 = 1:nfishes

                mx = round(max([cmbr{k1};fliplr(ambr{k1})]));
                
                a = poly2mask(cmbr{k1}(:,2),cmbr{k1}(:,1),mx(1),mx(2));
                b = poly2mask(ambr{k1}(:,1),ambr{k1}(:,2),mx(1),mx(2));
                
                intersection = (a & b);
                union = (a | b);

                Jaccard(k1) = sum(intersection(:))/sum(union(:));
                %Overlap(k1) = sum(intersection(:)) / min(sum(a(:)),sum(b(:)));
                Dice(k1) = (2*sum(intersection(:)))/(sum(union(:))+sum(intersection(:)));
                
            end
            
            %fill the BrainInfo Structure (per fish)
            temp = num2cell(Jaccard);[obj.BrainInfo(1:nfishes).Jaccard] = temp{:};
            temp = num2cell(Dice);[obj.BrainInfo(1:nfishes).Dice] = temp{:};
            
            %fill BrainStats Structure, statistical values
            obj.BrainStats.FiveNumberSummaryJaccard = sng_FiveNumberSum(Jaccard);
            obj.BrainStats.MeanJaccard = mean(Jaccard);
            
            obj.BrainStats.FiveNumberSummaryDice = sng_FiveNumberSum(Dice);
            obj.BrainStats.MeanDice = mean(Dice);       
            
        end     
                
        function obj = SpotVal(obj)
            
            load([obj.SavePath,'/',obj.InfoName,'.mat'],'SpotParameters')
                        
            nfishes = numel(obj.StackInfo);
            
            LinkDistance = cell(nfishes,1);
            CorrectSpots = cell(nfishes,1);
            nCorrect = zeros(nfishes,1);
            FalsePosSpots = cell(nfishes,1);
            nFalsePos = zeros(nfishes,1);
            FalseNegSpots = cell(nfishes,1);
            nFalseNeg = zeros(nfishes,1);
            Precision = zeros(nfishes,1);
            Recall = zeros(nfishes,1);
            F1score = zeros(nfishes,1);
            AbsDifference = zeros(nfishes,1);
            RelDifference = zeros(nfishes,1);
            
            SpotCom = cell(nfishes,1);
            ambsS = cell(nfishes,1);
                        
            %value added to prevent for substraction by zero        
            a = 0.001; 
            
            

            for k1 = 1:nfishes

                Spotpar = SpotParameters{k1};
                ambr = obj.Annotations(k1).MidBrain;
                ambs = obj.Annotations(k1).Spots;  

                
                if ~isempty(Spotpar)
                    
                    %select spot insite annotated brainregion
                    [spotcentroids] = reshape([Spotpar.Centroid],2,numel(Spotpar))';                               
                    [ins,~] = inpolygon(spotcentroids(:,1),spotcentroids(:,2),ambr(:,1),ambr(:,2));

                     SpotparS = Spotpar(...
                        ins' &...
                        [Spotpar.MinProbability] &... %selection of spotinfo based on conditions
                        [Spotpar.LargerThan] &...
                        [Spotpar.SmallerThan]);    

                    SpotCom{k1} = reshape([SpotparS.Centroid],2,numel(SpotparS))';
                else
                    SpotCom{k1} = [];
                end
                
                %because the annotated midbrain does not corresponds
                %exactly with the annotated spots, de annotated spots are
                %also filtered to compute the performance of the
                %spotdetection only
                [ins2,~] = inpolygon(ambs(:,1),ambs(:,2),ambr(:,1),ambr(:,2));
                ambsx = ambs(:,1);ambsx = ambsx(ins2);
                ambsy = ambs(:,2);ambsy = ambsy(ins2);                                                
                ambsS{k1} = [ambsx,ambsy];                

                [Correct,FalsePos,FalseNeg,link] = sng_CoordinateMatching...
                    (SpotCom{k1},ambsS{k1},10);                
                
                LinkDistance{k1} = link;

                if exist('Correct','var')
                    CorrectSpots{k1} = Correct;
                    nCorrect(k1) = size(Correct,1);
                else
                    CorrectSpots{k1} = [];
                    nCorrect(k1) = 0;
                end
                if exist('FalsePos','var')
                    FalsePosSpots{k1} = FalsePos;
                    nFalsePos(k1) = size(FalsePos,1);
                else
                    FalsePosSpots{k1} = [];
                    nFalsePos(k1) = 0;
                end
                if exist('FalseNeg','var')
                    FalseNegSpots{k1} = FalseNeg;
                    nFalseNeg(k1) = size(FalseNeg,1);
                else
                    FalseNegSpots{k1} = [];
                    nFalseNeg(k1) = 0;
                end

                Precision(k1) = nCorrect(k1)/(nCorrect(k1) + nFalsePos(k1)+ a);
                Recall(k1) = nCorrect(k1)/(nCorrect(k1) + nFalseNeg(k1)+ a);
                F1score(k1) = 2 * (Precision(k1)*Recall(k1))/(Precision(k1) + Recall(k1) + a);
                AbsDifference(k1) = nFalsePos(k1) - nFalseNeg(k1);
                RelDifference(k1) = (nFalsePos(k1) - nFalseNeg(k1)) / nCorrect(k1);
                

                obj.SpotStats.FiveNumberSummaryPrecision = sng_FiveNumberSum(Precision);
                obj.SpotStats.MeanPrecision = mean(Precision);

                obj.SpotStats.FiveNumberSummaryRecall = sng_FiveNumberSum(Recall);
                obj.SpotStats.MeanRecall = mean(Recall);                

                obj.SpotStats.FiveNumberSummaryF1score = sng_FiveNumberSum(F1score);
                obj.SpotStats.MeanF1score = mean(F1score); 
                
                obj.SpotStats.FiveNumberSummaryAbsDifference  = sng_FiveNumberSum(AbsDifference);
                obj.SpotStats.MeanAbsDifference = mean(AbsDifference); 

                obj.SpotStats.FiveNumberSummaryRelDifference  = sng_FiveNumberSum(RelDifference);
                obj.SpotStats.MeanRelDifference = mean(RelDifference); 


            end

            %store variables in obj.SpotInfo
            [obj.SpotInfo(1:nfishes).LinkDistance] = LinkDistance{:};
            [obj.SpotInfo(1:nfishes).CorrectSpots] = CorrectSpots{:};
            temp = num2cell(nCorrect);[obj.SpotInfo(1:nfishes).nCorrect] = temp{:};
            [obj.SpotInfo(1:nfishes).FalsePosSpots] = FalsePosSpots{:};
            temp = num2cell(nFalsePos); [obj.SpotInfo(1:nfishes).nFalsePos] = temp{:};
            [obj.SpotInfo(1:nfishes).FalseNegSpots] = FalseNegSpots{:};
            temp = num2cell(nFalseNeg); [obj.SpotInfo(1:nfishes).nFalseNeg] = temp{:};
            temp = num2cell(Precision); [obj.SpotInfo(1:nfishes).Precision] = temp{:};
            temp = num2cell(Recall); [obj.SpotInfo(1:nfishes).Recall] = temp{:};
            temp = num2cell(F1score); [obj.SpotInfo(1:nfishes).F1score] = temp{:};            
            temp = num2cell(AbsDifference); [obj.SpotInfo(1:nfishes).AbsDifference] = temp{:};
            temp = num2cell(RelDifference); [obj.SpotInfo(1:nfishes).RelDifference] = temp{:};
            %[obj.SpotSelection(1:nfishes).AnnotatedSpots] = ambsS{:};                
            %[obj.SpotSelection(1:nfishes).ComputedSpots] = SpotCom{:};
                       
        end
        
        function obj = SpotBrainVal(obj)
            
            %load([obj.SavePath,'/',obj.InfoName,'.mat'],'SpotParameters')
            load([obj.SavePath,'/',obj.InfoName,'.mat'],'SpotsDetected')
                                   
            
            nfishes = numel(obj.StackInfo);
        
            
            LinkDistance = cell(nfishes,1);
            CorrectSpots = cell(nfishes,1);
            nCorrect = zeros(nfishes,1);              
            FalsePosSpots = cell(nfishes,1);
            nFalsePos = zeros(nfishes,1);              
            FalseNegSpots = cell(nfishes,1);
            nFalseNeg = zeros(nfishes,1);
            Precision = zeros(nfishes,1);
            Recall = zeros(nfishes,1);
            F1score = zeros(nfishes,1);
            AbsDifference = zeros(nfishes,1);
            RelDifference = zeros(nfishes,1);            

            %SpotCom = cell(nfishes,1);

            %value added to prevent for substraction by zero        
            a = 0.001; 

            for k1 = 1:nfishes

                Spotpar = SpotsDetected{k1};
                %ambr = obj.BrainInfo(k1).AnnotatedMidBrainRegistrated;              

                SpotCom = reshape([Spotpar.Centroid],2,numel(Spotpar))';
                ambs = obj.Annotations(k1).Spots;  


                [Correct,FalsePos,FalseNeg,link] = sng_CoordinateMatching...
                    (SpotCom,ambs,10);

                LinkDistance{k1} = link;

                if exist('Correct','var')
                    CorrectSpots{k1} = Correct;
                    nCorrect(k1) = size(Correct,1);
                else
                    CorrectSpots{k1} = [];
                    nCorrect(k1) = 0;
                end
                if exist('FalsePos','var')
                    FalsePosSpots{k1} = FalsePos;
                    nFalsePos(k1) = size(FalsePos,1);
                else
                    FalsePosSpots{k1} = [];
                    nFalsePos(k1) = 0;
                end
                if exist('FalseNeg','var')
                    FalseNegSpots{k1} = FalseNeg;
                    nFalseNeg(k1) = size(FalseNeg,1);
                else
                    FalseNegSpots{k1} = [];
                    nFalseNeg(k1) = 0;
                end

                Precision(k1) = nCorrect(k1)/(nCorrect(k1) + nFalsePos(k1)+ a);
                Recall(k1) = nCorrect(k1)/(nCorrect(k1) + nFalseNeg(k1)+ a);
                F1score(k1) = 2 * (Precision(k1)*Recall(k1))/(Precision(k1) + Recall(k1) + a);
                AbsDifference(k1) = nFalsePos(k1) - nFalseNeg(k1);
                RelDifference(k1) = (nFalsePos(k1) - nFalseNeg(k1)) / nCorrect(k1);

                obj.SpotBrainStats.FiveNumberSummaryPrecision = sng_FiveNumberSum(Precision);
                obj.SpotBrainStats.MeanPrecision = mean(Precision);

                obj.SpotBrainStats.FiveNumberSummaryRecall = sng_FiveNumberSum(Recall);
                obj.SpotBrainStats.MeanRecall = mean(Recall);                

                obj.SpotBrainStats.FiveNumberSummaryF1score = sng_FiveNumberSum(F1score);
                obj.SpotBrainStats.MeanF1score = mean(F1score);
                
                obj.SpotBrainStats.FiveNumberSummaryAbsDifference  = sng_FiveNumberSum(AbsDifference);
                obj.SpotBrainStats.MeanAbsDifference = mean(AbsDifference); 

                obj.SpotBrainStats.FiveNumberSummaryRelDifference  = sng_FiveNumberSum(RelDifference);
                obj.SpotBrainStats.MeanRelDifference = mean(RelDifference); 

            end

            %store variables in obj.SpotInfo
            [obj.SpotBrainInfo(1:nfishes).LinkDistance] = LinkDistance{:};
            [obj.SpotBrainInfo(1:nfishes).CorrectSpots] = CorrectSpots{:};
            temp = num2cell(nCorrect);[obj.SpotBrainInfo(1:nfishes).nCorrect] = temp{:};
            [obj.SpotBrainInfo(1:nfishes).FalsePosSpots] = FalsePosSpots{:};
            temp = num2cell(nFalsePos); [obj.SpotBrainInfo(1:nfishes).nFalsePos] = temp{:};
            [obj.SpotBrainInfo(1:nfishes).FalseNegSpots] = FalseNegSpots{:};
            temp = num2cell(nFalseNeg); [obj.SpotBrainInfo(1:nfishes).nFalseNeg] = temp{:};
            temp = num2cell(Precision); [obj.SpotBrainInfo(1:nfishes).Precision] = temp{:};
            temp = num2cell(Recall); [obj.SpotBrainInfo(1:nfishes).Recall] = temp{:};
            temp = num2cell(F1score); [obj.SpotBrainInfo(1:nfishes).F1score] = temp{:};
            temp = num2cell(AbsDifference); [obj.SpotBrainInfo(1:nfishes).AbsDifference] = temp{:};
            temp = num2cell(RelDifference); [obj.SpotBrainInfo(1:nfishes).RelDifference] = temp{:};                       
        end

        function ShowBoxPlot(obj,exportit)
            %Example show and save
            %   show(obj,1)
            %Example only show
            %   obj.show
            
            %%% brainval spotval
            fsx = 6;fsy = 10;
                      
            if isfield(obj.BrainInfo,'Jaccard')
                [h1,g1] = setfigax1;                    
                boxplot(g1,[...
                    obj.BrainInfo.Jaccard;...
                    obj.BrainInfo.Dice]',{'Jaccard','Dice'});
                title('Brain Validation')
                setfigax2(h1,g1)            
            end
            if isfield(obj.SpotInfo,'Precision')
                [h2,g2] = setfigax1;
                boxplot(g2,[...
                    obj.SpotInfo.Precision;...
                    obj.SpotInfo.Recall;...
                    obj.SpotInfo.F1score]',{'Precision','Recall','F1score'});
                title('Spot Validation')            
                setfigax2(h2,g2)
            end
            
            if isfield(obj.SpotBrainInfo,'Precision')
                [h3,g3] = setfigax1;      
                boxplot(g3,[...
                    obj.SpotBrainInfo.Precision;...
                    obj.SpotBrainInfo.Recall;...
                    obj.SpotBrainInfo.F1score]',{'Precision','Recall','F1score'});
                title('Spot Validation on Computed Brain')          
                setfigax2(h3,g3)
            end

            if isfield(obj.SpotInfo,'AbsDifference')
                [h4,g4] = setfigax1;      
                boxplot(g4,[...
                    obj.SpotInfo.AbsDifference]',{'AbsDifference'});
                title('Spot Validation')        
                
                set(g4,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
                set(g4,'Units','centimeters','Position',[1.2 1.2 fsx-1.7 fsy-1.7])
                
                if exist('exportit','var') && exportit
                    export_fig(h4 ,['/Users/samuelgeurts/Desktop/','boxplot',num2str(h4.Number)], '-png', '-r600', '-nocrop');
                end
                
                ScaledFigure.calibrateDisplay(113.6) %113.6 for this screen
                ScaledFigure(h4,'reuse')
                set(h4,'Units','Centimeters')
                set(h4,'Position',(get(h4,'Position') + [fsx 0 0 0]))  
            end

            if isfield(obj.SpotBrainInfo,'AbsDifference')
                [h5,g5] = setfigax1;      
                boxplot(g5,[...
                    obj.SpotBrainInfo.AbsDifference]',{'AbsDifference'});
                title('Spot Validation on Computed Brain')        
                
                set(g5,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
                set(g5,'Units','centimeters','Position',[1.2 1.2 fsx-1.7 fsy-1.7])
                
                if exist('exportit','var') && exportit
                    export_fig(h5 ,['/Users/samuelgeurts/Desktop/','boxplot',num2str(h5.Number)], '-png', '-r600', '-nocrop');
                end
                
                ScaledFigure.calibrateDisplay(113.6) %113.6 for this screen
                ScaledFigure(h5,'reuse')
                set(h5,'Units','Centimeters')
                set(h5,'Position',(get(h5,'Position') + [fsx 0 0 0]))  
            end
            
            if isfield(obj.SpotInfo,'RelDifference')
                [h6,g6] = setfigax1;      
                boxplot(g6,[...
                    obj.SpotInfo.RelDifference]',{'RelDifference'});
                title('Spot Validation')        
                
                set(g6,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
                set(g6,'Units','centimeters','Position',[1.2 1.2 fsx-1.7 fsy-1.7])
                
                if exist('exportit','var') && exportit
                    export_fig(h6 ,['/Users/samuelgeurts/Desktop/','boxplot',num2str(h6.Number)], '-png', '-r600', '-nocrop');
                end
                
                ScaledFigure.calibrateDisplay(113.6) %113.6 for this screen
                ScaledFigure(h6,'reuse')
                set(h6,'Units','Centimeters')
                set(h6,'Position',(get(h6,'Position') + [fsx 0 0 0]))  
            end
            
            if isfield(obj.SpotBrainInfo,'RelDifference')
                [h7,g7] = setfigax1;      
                boxplot(g7,[...
                    obj.SpotBrainInfo.RelDifference]',{'RelDifference'});
                title('Spot Validation on Computed Brain')        
                
                set(g7,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
                set(g7,'Units','centimeters','Position',[1.2 1.2 fsx-1.7 fsy-1.7])
                
                if exist('exportit','var') && exportit
                    export_fig(h7 ,['/Users/samuelgeurts/Desktop/','boxplot',num2str(h7.Number)], '-png', '-r600', '-nocrop');
                end
                
                ScaledFigure.calibrateDisplay(113.6) %113.6 for this screen
                ScaledFigure(h7,'reuse')
                set(h7,'Units','Centimeters')
                set(h7,'Position',(get(h7,'Position') + [fsx 0 0 0]))  
            end            
            
            %%
            function [figurehandle,axishandle] = setfigax1
                %fsx = 5;fsy = 10;
                figurehandle = figure('PaperUnits','centimeters','Color',[1 1 1]);
                axishandle = gca;
                sng_figcm(fsx,fsy);
            end
            function setfigax2(figurehandle,axishandle)
                set(axishandle,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
                set(axishandle,'Units','centimeters','Position',[1.2 1.2 fsx-1.7 fsy-1.7]);
                set(axishandle,'YLim',[-0.02,1.02]);
                
                if exist('exportit','var') && exportit
                    export_fig(figurehandle ,['/Users/samuelgeurts/Desktop/','boxplot',num2str(figurehandle.Number)], '-png', '-r600', '-nocrop');
                end
                
                ScaledFigure.calibrateDisplay(113.6); %113.6 for this screen
                ScaledFigure(figurehandle,'reuse');
                set(figurehandle,'Units','Centimeters');
                set(figurehandle,'Position',(get(figurehandle,'Position') + [fsx 0 0 0]));             
            end

                             
        end        
        
        function ShowFishVal(obj,fishnumbers,ComOrAnn)
            %if ComOrAnn is not given, the computed brain is used
            
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif (max(fishnumbers) > numel(obj.StackInfo))
                error('at least one fish does not exist in StackInfo')
            end        
            
            %CompOrAnn determines if the computed or the annotated brain is
            %displayed inlcudin ghte Correct/FalsePos/FalseNeg
                
            load([obj.SavePath,'/',obj.InfoName,'.mat'],'BrainSegmentationInfo')

            if ~isempty({obj.Annotations.MidBrain})
                ambr = {obj.Annotations.MidBrain};
            end
            
            
            for k1 = 1:numel(fishnumbers)               
                fn = fishnumbers(k1);
                
                cmbr = BrainSegmentationInfo(fn).BrainEdge;
                AlignedFish = imread([obj.SavePath,'/','AlignedFish','/',obj.StackInfo(fn).stackname,'.tif']);                         

%%                    
                figure;imshow(uint8(AlignedFish))                    
                hold on;
                plot(cmbr(:,2),cmbr(:,1),'Color',[255 75 75]/255,'LineWidth',2);
                
                if exist('ambr','var')
                    plot(ambr{fn}(:,1),ambr{fn}(:,2),'Color',[75 75 255]/255,'LineWidth',2);
                end
                    
                if ~exist('ComOrAnn','var') || strcmp(ComOrAnn,'Com')
                    Correct = obj.SpotBrainInfo(fn).CorrectSpots;
                    FalsePos = obj.SpotBrainInfo(fn).FalsePosSpots;
                    FalseNeg = obj.SpotBrainInfo(fn).FalseNegSpots;
                    
                    ps = obj.SpotBrainInfo(fn).Precision; %used for textbox
                    rs = obj.SpotBrainInfo(fn).Recall;
                    fs = obj.SpotBrainInfo(fn).F1score;
                                        
                elseif strcmp(ComOrAnn,'Ann')
                    Correct = obj.SpotInfo(fn).CorrectSpots;
                    FalsePos = obj.SpotInfo(fn).FalsePosSpots;
                    FalseNeg = obj.SpotInfo(fn).FalseNegSpots;
                                      
                    ps = obj.SpotInfo(fn).Precision; %used for textbox
                    rs = obj.SpotInfo(fn).Recall;
                    fs = obj.SpotInfo(fn).F1score;
                    
                end
                                               
                %scatter(Correct(:,1), Correct(:,2),400,'LineWidth',2,'MarkerEdgeColor',1/255*[75 255 75]);
                %scatter(FalsePos(:,1), FalsePos(:,2),400,'LineWidth',2,'MarkerEdgeColor',1/255*[255 75 75]);
                %scatter(FalseNeg(:,1), FalseNeg(:,2),400,'LineWidth',2,'MarkerEdgeColor',1/255*[75 75 255]);
                %legend('Computed brain','Annotated brain','Correct','False Positive','False Negative')

                scatter([Correct(:,1);FalsePos(:,1)], [Correct(:,2);FalsePos(:,2)],400,'LineWidth',2,'MarkerEdgeColor',1/255*[255 75 75]);             
                scatter([Correct(:,1);FalseNeg(:,1)], [Correct(:,2);FalseNeg(:,2)],270,'LineWidth',2,'MarkerEdgeColor',1/255*[75 75 255]);             
                legend('Computed brain','Annotated brain','Computed spots','Annotated spots')
                
                
                nc = size(Correct,1) + size(FalsePos,1); %number of computed spots
                na = size(Correct,1) + size(FalseNeg,1); %number of annotated spots
                jb = obj.BrainInfo(val(k)).Jaccard; %Jaccard number
                db = obj.BrainInfo(val(k)).Dice; % Dice number
                
                                
                str = sprintf(['BRAIN PARAMETERS',...
                    '\n','Jaccard: %.2f',...
                    '\n','dice: %.2f',...
                    '\n\n','SPOT PARAMETERS',...
                    '\n','Precision: %.2f',...
                    '\n','Recall: %.2f',...
                    '\n','F1score: %.2f',...
                    '\n\n','Computed spots: %.0f',...
                    '\n','Annotated spots: %.0f',...                    
                    '\n\n','Correct spots: %.0f',...
                    '\n','False Positives: %.0f',...
                    '\n','False Negatives: %.0f'],...
                    jb,db,ps,rs,fs,nc,na,...
                    size(Correct,1),...
                    size(FalsePos,1),...
                    size(FalseNeg,1));             
                
                annotation('textbox',[.1 .63 .7 .3],'String',str,'FitBoxToText','on');
           
            end 
        end
        
        function ShowFish(obj,fishnumbers)
            %if ComOrAnn is not given, the computed brain is used
            
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif (max(fishnumbers) > numel(obj.StackInfo))
                error('at least one fish does not exist in StackInfo')
            end        
            
            %CompOrAnn determines if the computed or the annotated brain is
            %displayed inlcudin ghte Correct/FalsePos/FalseNeg
                
            load([obj.SavePath,'/',obj.InfoName,'.mat'],'BrainSegmentationInfo')
            load([obj.SavePath,'/',obj.InfoName,'.mat'],'SpotsDetected')

            
            for k1 = 1:numel(fishnumbers)               
                fn = fishnumbers(k1);
                
                cmbr = BrainSegmentationInfo(fn).BrainEdge;
                AlignedFish = imread([obj.SavePath,'/','AlignedFish','/',obj.StackInfo(fn).stackname,'.tif']);                
                cmbs = reshape([SpotsDetected{fn}.Centroid],2,numel(SpotsDetected{fn}))';
%%                    
                figure;imshow(uint8(AlignedFish))                    
                hold on;
                plot(cmbr(:,2),cmbr(:,1),'Color',[255 75 75]/255,'LineWidth',2);
                                               
                scatter(cmbs(:,1), cmbs(:,2),400,'LineWidth',2,'MarkerEdgeColor',1/255*[255 75 75]);             
                               
                nc = size(cmbs,1) %number of computed spots                             
                str = sprintf('Computed Spots: %.0f',nc)         
                a = annotation('textbox',[.1 .63 .7 .3],'String',str,'FitBoxToText','on','FontSize',15);
           
            end 
        end
        
        function saveit(obj)
            %obj.savedate = datetime;
            %dt = strcat(string(year(date)),string(month(date)),string(day(date)));                     
            %firstim = obj.ImageInfo(1).name;
            %firstn = char(regexp(firstim,'\d+.tif','match'));  
            %firstn = strrep(firstn,'.tif','');
            %name = strrep(firstim, [firstn,'.tif'],'');

            save(strcat(obj.SavePath,'/',obj.SaveName,'.mat'),'obj');    
        end
        
       
    end
end

