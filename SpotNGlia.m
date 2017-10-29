classdef SpotNGlia

    properties   
        OriginalPath = []
        TemplatePath3dpf = []
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
                        
        CorrectedFish = []
        CombinedFish = []
        AlignedFish = []
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
                [obj.OriginalPath,BasePath] = PathsComplete('or','bp');            
                 obj.SavePath = [BasePath,'/','SpotNGlia'];                           
                [obj.RoiMicrogliaPath,obj.RoiBrainPath] = PathsComplete('rois','roib');
            elseif exist('mode1','var') && mode1 == 9 %temporary mode with temp folder
                [obj.OriginalPath,BasePath] = PathsComplete('or','bp');            
                 obj.SavePath = [BasePath,'/','SpotNGlia'];
                 obj.OriginalPath = [BasePath,'/','originalTEMP'];               
                [obj.RoiMicrogliaPath,obj.RoiBrainPath] = PathsComplete('rois','roib');
            else %open new folder
                disp('select source folder')
                obj.OriginalPath = uigetdir;
                disp('select save folder')
                obj.SavePath = uigetdir;
            end
            
            obj.TemplatePath3dpf =  PathsComplete('tp3');            
                      
            zfinput = struct('stage',[],'substage',[],'name',[]','value',[],'sensitivity',[]);            
            zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','levels',1,'low'); %number of scaled levels correlation is performed
            zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','iterations',10,'low'); %number of iterations iat is performed
            zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','scale',1/16,'low'); %lowered initial scale to increase computation speed for correlation
            zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','threshold',0.96,'severe'); %threshold for correlation selection
            zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','threshold2',2,'severe'); %threshold for warp modulus
            zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','sorting','','high'); %date or name  
            zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','onoff',false,''); %sigma for bandpassfilter (lowpass)    
            zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','sigmalp',1,'moderate'); %sigma for bandpassfilter (lowpass)
            zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','sigmahp',4,'moderate'); %for bandpassfilter (highpass)
            zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','scaleC',1/4,'low'); %lowered initial scale to increase computation speed for correlation
            zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','levelsC',2,'low'); %number of scaled levels correlation is performed
            zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','iterationsC',10,'low'); %number of iterations iat is performed
            zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','scaleS',1/4,'low'); %lowered initial scale to increase computation speed for correlation
            zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','levelsS',3,'low'); %number of scaled levels correlation is performed
            zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','iterationsS',20,'low'); %iterations iat is performed
            zfinput = sng_zfinput(zfinput,0,'ExtendedDeptOfField','edof','variancedisksize',7,'moderate'); %sigma for bandpassfilter (lowpass)    
            zfinput = sng_zfinput(zfinput,0,'Registration','BackgroundRemoval','Method','TriangleSmooth','?');
            zfinput = sng_zfinput(zfinput,0,'Registration','BackgroundRemoval','Smooth',1,'?');  
            zfinput = sng_zfinput(zfinput,0,'Registration','BackgroundRemoval','ChannelMethod','cuboid','?');  
            zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','AngleSteps1',100,'?');  
            zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','Scale1',1/16,'?');  
            zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','AngleRange2',0.01,'?');  
            zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','AngleSteps2',50,'?');
            zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','FinalCenteredYCrop',1000,'?');
            zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','ScaleRange',[0.6,1.3],'?');      
            zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','ScaleSteps',200,'?');      
            zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','Scale2',1/16,'?');      
            zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','RemoveTail',600,'');   
            zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','Scale3',1/4,''); %lowered initial scale to increase computation speed for correlation
            zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','AffineMethod','translation',''); %affine or translation
            zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','Levels',3,''); %number of scaled levels correlation is performed
            zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','Iterations',40,''); %number of iterations iat is performed
            zfinput = sng_zfinput(zfinput,0,'BrainSegmentation','','Method','TriangleSmooth','?');
            zfinput = sng_zfinput(zfinput,0,'BrainSegmentation','','Method2',4,'?');            
            zfinput = sng_zfinput(zfinput,0,'SpotDetection','RgbToGray','ColorToGrayVector',[0;1;0],'');  %select color channel [0 1 0], 
            zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleBase',0.5,'');  
            zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleLevels',7,'');  
            zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','Kthreshold',0,'');  
            zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPlevels',5:7,'');  
            zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPthreshold',200,'');  
            zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MinSpotSize',8.4,''); % size selection 
            zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MaxSpotSize',458,''); % size selection 
            zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MinProbability',0.066,''); %color selection  
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
                dirinfo =  dir([obj.OriginalPath,'/*.','tif']);
                imageinfotemp = rmfield(dirinfo,{'isdir','datenum','bytes'}); %removes unimportant fields

                % make en selection of fishslices if fisnumbers is given as input
                if ~exist('fishnumbers','var')
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
        end
        
        function obj = PreProcession(obj,fishnumbers)
            
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
           
                for k1 = 1:nfishes 
                    %select slices
                    fn = fishnumbers(k1);
                                        
                    ImageSlice = cell(1,obj.StackInfo(fn).stacksize);%preallocate for every new slice
                    for k2 = 1:obj.StackInfo(fn).stacksize
                        ImageSlice{k2} = imread([obj.OriginalPath,'/',obj.StackInfo(fn).imagenames{k2}]);
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
                    
                    obj.CorrectedFish(k1).image = CorrectedFishTemp;
        
%TODO for imaging   [CorrectedSlice,ppoutput,ImageSliceCor,FiltIm] = PreprocessionLink(ImageSlice,zfinput)

%                    if savefig2_TF
%                        sng_SaveCell2TiffStack(CorrectedSlice,[obj.SavePath,'/',StackInfo(k1).stackname,'.tif'])
%                    end
                end
                
                obj.PreprocessionInfo = ppoutput;
                
            end       
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
                 ImageSlice = imread([obj.OriginalPath,'/',obj.StackInfo(fn).imagenames{k2}]);
                 [ImageSliceCor{k2},~] = sng_RGB_IATwarp2(ImageSlice(:,:,1:3),obj.PreprocessionInfo(k1).ColorWarp{1});
                 %[cellImg1{k2}, ~] = iat_inverse_warping(ImageSliceCor{k2}, obj.PreprocessionInfo(k1).SliceWarp{k2}, par.transform, 1:N, 1:M);
                 end
             end
             %}                
            if nfishes > 1
                    edoutput(nfishes) = struct('IndexMatrix',[],...
                        'variance_sq',[]);
            end
             
            for k1 = 1:nfishes
                waitbar(nfishes/k1,h,'Extended Dept of Field')                
                fn = fishnumbers(k1);
                [CombinedFishTemp,edoutput(k1)] = ExtendedDeptofFieldLink2(obj.CorrectedFish(fn).image,obj.ZFParameters);                
                obj.CombinedFish(k1).image = CombinedFishTemp;
            end
            obj.ExtendedDeptOfFieldInfo = edoutput;

        end
        
        function obj = Registration(obj,fishnumbers)
            
            if isempty(obj.CompleteTemplate)
                obj.CompleteTemplate = LoadTemplateLink3(obj.TemplatePath3dpf);                         
            end
            
            
            if ~exist('fishnumbers','var')
                fishnumbers = 1:numel(obj.ExtendedDeptOfFieldInfo);
            elseif max(fishnumbers) > numel(obj.ExtendedDeptOfFieldInfo)
                error('at least one fish does not exist in ExtendedDeptOfFieldInfo')                   
            end
            
            nfishes = numel(fishnumbers);             
            obj.fishnumbers.registration = fishnumbers;
            
            rgoutput = cell(nfishes,1); 
            for k1 = 1:nfishes                    
                fn = fishnumbers(k1);
                [AlignedFishTemp,rgoutput{k1,1}] = AllignmentLink5(obj.CombinedFish(fn).image,obj.CompleteTemplate,obj.ZFParameters);
                obj.AlignedFish(k1).image = AlignedFishTemp;
            end
            obj.RegistrationInfo = rgoutput;
        end
        
        function obj = BrainSegmentation(obj,fishnumbers)
            
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
                [~,broutput(k1)] = MidBrainDetectionLink3(obj.AlignedFish(fn).image,obj.CompleteTemplate,obj.ZFParameters);      
                %obj.BrainFish(k1).image = BrainFishTemp;
            end
            obj.BrainInfo = broutput;
        end            

        function obj = SpotDetection(obj,fishnumbers)
            
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
                [~,spoutput{k1,1}] = SpotDetectionLink2(obj.AlignedFish(fn).image,obj.CompleteTemplate,obj.BrainInfo(k1).BrainEdge,obj.ZFParameters);   
                %obj.BrainFish(k1).image = BrainFishTemp;
            end
            obj.SpotInfo = spoutput;
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
            name = strrep(firstim, [firstn,'.tif'],'')
            dt = strcat(string(year(date)),string(month(date)),string(day(date)));            
            save(strcat(obj.SavePath,'/SNG_',dt,'_',name,'.mat'),'obj');      
        end
                
        
    end
end

