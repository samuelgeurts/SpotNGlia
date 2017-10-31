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
%{      
    properties
        Numb = []        
        date = []       
        
        SpotParameters = []

      

        BrainInfo = []
        BrainStats = []

        SpotStats = []
        SpotSelection = []      
        
        SpotBrainInfo = []
        SpotBrainStats = []    
    end 
%}
    
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
                
                ImageInfo = obj.ImageInfo;
                StackInfo = obj.StackInfo;
                
                save([obj.SavePath,'/',obj.InfoName,'.mat'],'ImageInfo')
                save([obj.SavePath,'/',obj.InfoName,'.mat'],'StackInfo','-append')
            end
            
            if ~obj.ImageInfoChecked_TF
                
                openvar('obj.ImageInfo')
                openvar('obj.StackInfo')                
                
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
            
            h = waitbar(0,'Combine fish slices','Name','SpotNGlia');  
            
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.PreprocessionInfo);
            elseif max(fishnumbers) > numel(obj.PreprocessionInfo)
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
                waitbar(k1/nfishes,h,'Extended Dept of Field')                
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
            
            h = waitbar(0,'Combine fish slices','Name','SpotNGlia');  
            
            if isempty(obj.CompleteTemplate)
                obj.CompleteTemplate = LoadTemplateLink3([obj.SourcePath,'/','Template 3 dpf']);                         
            end
                        
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.ExtendedDeptOfFieldInfo);
            elseif max(fishnumbers) > numel(obj.ExtendedDeptOfFieldInfo)
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
                waitbar(k1/nfishes,h,'Align to template')                                
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
            h = waitbar(0,'Combine fish slices','Name','SpotNGlia');  
            
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.RegistrationInfo);
            elseif max(fishnumbers) > numel(obj.RegistrationInfo)
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
                waitbar(k1/nfishes,h,'Brain Segmentation')
                AlignedFish = imread([obj.SavePath,'/','AlignedFish','/',obj.StackInfo(fn).stackname,'.tif']);                         
                [~,BrainSegmentationInfo(k1)] = MidBrainDetectionLink3(AlignedFish,obj.CompleteTemplate,obj.ZFParameters);      
                %obj.BrainFish(k1).image = BrainFishTemp;
            end
            obj.saveit            
            save([obj.SavePath,'/',obj.InfoName,'.mat'],'BrainSegmentationInfo','-append')  

            delete(h)            
        end            

        function obj = SpotDetection(obj,fishnumbers)
            h = waitbar(0,'Combine fish slices','Name','SpotNGlia');  
            
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.BrainSegmentationInfo);
            elseif max(fishnumbers) > numel(obj.BrainSegmentationInfo)
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
            
            for k1 = 1:nfishes                    
                fn = fishnumbers(k1);
                waitbar(k1/nfishes,h,'Spot Detection')
                AlignedFish = imread([obj.SavePath,'/','AlignedFish','/',obj.StackInfo(fn).stackname,'.tif']);                         
                [~,SpotDetectionInfo{k1,1}] = SpotDetectionLink2(AlignedFish,obj.CompleteTemplate,obj.BrainSegmentationInfo(k1).BrainEdge,obj.ZFParameters);   
                %obj.BrainFish(k1).image = BrainFishTemp;
            end
            obj.saveit
            save([obj.SavePath,'/',obj.InfoName,'.mat'],'SpotDetectionInfo','-append')  

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
            
        function saveit(obj)
            %obj.savedate = datetime;
            %dt = strcat(string(year(date)),string(month(date)),string(day(date)));                     
            %firstim = obj.ImageInfo(1).name;
            %firstn = char(regexp(firstim,'\d+.tif','match'));  
            %firstn = strrep(firstn,'.tif','');
            %name = strrep(firstim, [firstn,'.tif'],'');

            save(strcat(obj.SavePath,'/',obj.SaveName,'.mat'),'obj');    
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
            if ~isempty('obj.AnnotatedMidBrainPath')
                disp('Select Annotated MidBrain Path')
                obj.AnnotatedMidBrainPath = uigetdir([],'Select Annotated MidBrain Path');
                ambr = cell(nfishes,1);               
                for k1 = 1:nfishes
                    RoiBrain = ReadImageJROI([obj.SourcePath,'/Roi brain/',obj.StackInfo(k1).stackname,'.zip']);

                    amb = sng_roicell2poly(RoiBrain,1);   
                    [ambr{k1}(:,1),ambr{k1}(:,2)] = transformPointsForward(tform_1234{k1},amb(:,1),amb(:,2));
                    ambr{k1} = double(ambr{k1});           
                end                                    
                [obj.Annotations(1:nfishes).MidBrain] = ambr{:};    
            end
            
            %annotated spots            
            if ~isempty('obj.AnnotatedSpotPath')
                disp('Select Annotated Spot Path')
                obj.AnnotatedSpotPath = uigetdir([],'Select Annotated Spot Path');
                     SpotAnn = cell(nfishes,1);
                for k1 = 1:nfishes
                    RoiMicroglia = ReadImageJROI([obj.SourcePath,'/Roi microglia/',obj.StackInfo(k1).stackname,'.roi']);

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
        
 %{
       function ShowBrain(obj,fishnumbers,ComOrAnn)
            
            %CompOrAnn determines if the computed or the annotated brain is
            %displayed inlcudin ghte Correct/FalsePos/FalseNeg
                
            load([obj.SavePath,'/',obj.InfoName,'.mat'],'BrainSegmentationInfo')
            
            temp = load([obj.SpotPath,'/stackinfo.mat']);
            stackinfo = temp.stackinfo;    
            
            for k = 1:numel(fishnumbers)               
                fn = fishnumers(k);
                
                %tform_1234 = stackinfo(k1).Registration(strcmp({stackinfo(k1).Registration.name},'tform_complete')).value; 

% !!!!!       %%%%Ga verder met checken of Ialigned bestaat en anders berekenen.
%                if ~exist([obj.SavePath,'/AlignedFish'],'dir')
%                    if exist([obj.SavePath,'/Fish'],'dir')
%                
%                CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',[obj.ImageInfo(k1).Name],'.tif']);           
%                Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k1).ExtendedDeptOfField.IndexMatrix);
%                Ialigned = imwarp(Icombined,tform_1234,'FillValues',255,'OutputView',CompleteTemplate.ref_temp);
%
% misschien nog eens toepassen maar dan moet CombinedInfo kleiner en vanaf
% original fishes berekenen, wel mooi idee ->github

                cmbr = obj.BrainSegmentationInfo(k1).ComputedMidBrain;
                ambr = obj.BrainInfo(k1).AnnotatedMidBrainRegistrated;                    
                
                    %figure;imshow(uint8(Icombined))
                    %hold on;plot(mbr(:,1),mbr(:,2))
%%                    
                figure;imshow(uint8(Ialigned))                    
                hold on;
                plot(cmbr(:,2),cmbr(:,1),'Color',[255 75 75]/255,'LineWidth',2);                    
                plot(ambr(:,1),ambr(:,2),'Color',[75 75 255]/255,'LineWidth',2);   
                if strcmp(ComOrAnn,'Com')
                    Correct = obj.SpotBrainInfo(k1).CorrectSpots;
                    FalsePos = obj.SpotBrainInfo(k1).FalsePosSpots;
                    FalseNeg = obj.SpotBrainInfo(k1).FalseNegSpots;
                    
                    ps = obj.SpotBrainInfo(val(k)).Precision; %used for textbox
                    rs = obj.SpotBrainInfo(val(k)).Recall;
                    fs = obj.SpotBrainInfo(val(k)).F1score;
                                        
                elseif strcmp(ComOrAnn,'Ann')
                    Correct = obj.SpotInfo(k1).CorrectSpots;
                    FalsePos = obj.SpotInfo(k1).FalsePosSpots;
                    FalseNeg = obj.SpotInfo(k1).FalseNegSpots;
                                      
                    ps = obj.SpotInfo(val(k)).Precision; %used for textbox
                    rs = obj.SpotInfo(val(k)).Recall;
                    fs = obj.SpotInfo(val(k)).F1score;
                    
                else
                    Error('choose between choose between Ann and Com, for displaying correct/falsepos/falseneg')
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
        %}
        
    end
end

