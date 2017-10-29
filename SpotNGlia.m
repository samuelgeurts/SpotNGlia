classdef SpotNGlia

    properties   
        FishPath = []
        SourcePath = []
        SavePath = []
        
        fishnumbers = []
        slicenumbers = []
        savedate = []
        
        ZFParameters = []        

        CompleteTemplate = []
        RoiMicrogliaPath = []    
        RoiBrainPath = []  
        
        StackInfo = []
        ImageInfo = []                        
        PreprocessionInfo = []        
        ExtendedDeptOfFieldInfo = []        
        RegistrationInfo = []
        BrainSegmentationInfo = []
        SpotInfo = []
                        
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
            else
                disp('select image folder to process')
                obj.FishPath = uigetdir([],'select image folder to process');
                disp('select destination folder')
                obj.SavePath = uigetdir([],'select destination folder');
                disp('select program source folder')
                obj.SourcePath = uigetdir([],'select program source folder');                
            end            
            
            load([obj.SourcePath,'/','zfinput.mat']);   
            obj.ZFParameters = zfinput;
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
                obj.slicenumbers.stackinfo = find([obj.ImageInfo.CorNextStack] ~= 2)
                
                if ~exist('fishnumbers','var')
                    fishnumbers = 1:numel(obj.StackInfo);
                elseif (max(fishnumbers) > numel(obj.StackInfo))
                    error('at least one fish does not exist in StackInfo')
                end
                nfishes = numel(fishnumbers);                    
                obj.fishnumbers.preprocession = fishnumbers;   

                %only preallocate if more than 1 fish has to be processed
                if nfishes > 1
                    ppoutput(nfishes) = struct('ColorWarp',[],...
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
                        [CorrectedFishTemp,ppoutput] = PreprocessionLink(ImageSlice,obj.ZFParameters);               
                    else
                        [CorrectedFishTemp,ppoutput(k1)] = PreprocessionLink(ImageSlice,obj.ZFParameters);
                    end

                    sng_SaveCell2TiffStack(CorrectedFishTemp,[TempFolderName,'/',obj.StackInfo(k1).stackname,'.tif'])
                                       
                    %obj.CorrectedFish(k1).image = CorrectedFishTemp; <- %save to object takes to much space
        
%TODO for imaging   [CorrectedSlice,ppoutput,ImageSliceCor,FiltIm] = PreprocessionLink(ImageSlice,zfinput)
%                    if savefig2_TF
%                        sng_SaveCell2TiffStack(CorrectedSlice,[obj.SavePath,'/',StackInfo(k1).stackname,'.tif'])
%                    end
                end
                
                obj.PreprocessionInfo = ppoutput;
                
            end       
            obj.saveit
            delete(h)
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
                edoutput(nfishes) = struct('IndexMatrix',[],...
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
                [CombinedFishTemp,edoutput(k1)] = ExtendedDeptofFieldLink2(CorrectedFish,obj.ZFParameters);                
                %obj.CombinedFish(k1).image = CombinedFishTemp; because it takes to much space               
                imwrite(uint8(CombinedFishTemp),[TempFolderName,'/',obj.StackInfo(k1).stackname,'.tif'],...
                    'WriteMode','overwrite','Compression','none');
            end
            obj.ExtendedDeptOfFieldInfo = edoutput;
            obj.saveit
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
                                        
            rgoutput = cell(nfishes,1); 
            for k1 = 1:nfishes
                waitbar(k1/nfishes,h,'Align to template')                                
                fn = fishnumbers(k1);                
                CombinedFish = imread([obj.SavePath,'/','CombinedFish','/',obj.StackInfo(fn).stackname,'.tif']);          
                [AlignedFishTemp,rgoutput{k1,1}] = AllignmentLink5(CombinedFish,obj.CompleteTemplate,obj.ZFParameters);
                %obj.AlignedFish(k1).image = AlignedFishTemp; %because it
                %takes to much space                                
                imwrite(uint8(AlignedFishTemp),[TempFolderName,'/',obj.StackInfo(k1).stackname,'.tif'],...
                    'WriteMode','overwrite','Compression','none');
            end
            obj.RegistrationInfo = rgoutput;
            obj.saveit
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
            	broutput(nfishes) = struct('EdgeFilterWidth',[],...
                    'ShortestPath',[],...
                    'ShortestPathValue',[],...
                    'BrainEdge',[]);
            end            
            
            for k1 = 1:nfishes                    
                fn = fishnumbers(k1);
                waitbar(k1/nfishes,h,'Brain Segmentation')
                AlignedFish = imread([obj.SavePath,'/','AlignedFish','/',obj.StackInfo(fn).stackname,'.tif']);                         
                [~,broutput(k1)] = MidBrainDetectionLink3(AlignedFish,obj.CompleteTemplate,obj.ZFParameters);      
                %obj.BrainFish(k1).image = BrainFishTemp;
            end
            obj.BrainInfo = broutput;
            obj.saveit
            delete(h)            
        end            

        function obj = SpotDetection(obj,fishnumbers)
            h = waitbar(0,'Combine fish slices','Name','SpotNGlia');  
            
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.BrainInfo);
            elseif max(fishnumbers) > numel(obj.BrainInfo)
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
                [~,spoutput{k1,1}] = SpotDetectionLink2(AlignedFish,obj.CompleteTemplate,obj.BrainInfo(k1).BrainEdge,obj.ZFParameters);   
                %obj.BrainFish(k1).image = BrainFishTemp;
            end
            obj.SpotInfo = spoutput;
            obj.saveit            
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
            obj = obj.SpotInfo(fishnumbers);            
            obj.saveit
            
        end
            
        function saveit(obj)
            obj.savedate = datetime;
            
            firstim = obj.ImageInfo(1).name;
            firstn = char(regexp(firstim,'\d+.tif','match'));  
            firstn = strrep(firstn,'.tif','');
            name = strrep(firstim, [firstn,'.tif'],'');
            dt = strcat(string(year(date)),string(month(date)),string(day(date)));            
            save(strcat(obj.SavePath,'/SNG_',dt,'_',name,'.mat'),'obj');      
        end        
    end
end

