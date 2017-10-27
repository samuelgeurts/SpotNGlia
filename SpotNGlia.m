classdef SpotNGlia

    properties   
        OriginalPath = []
        SavePath = []
        
        fishnumbers = []
        
        ZFParameters = []        
        
        StackInfo = []
        ImageInfo = []                        
        PreprocessionInfo = []
 
    end
    
    properties
        Numb = []        
        date = []       
        
        SpotParameters = []

        RoiMicrogliaPath = []    
        RoiBrainPath = []        

        BrainInfo = []
        BrainStats = []

        SpotInfo = []
        SpotStats = []
        SpotSelection = []      
        
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
                [obj.OriginalPath,BasePath] = PathsComplete('or','bp');            
                 obj.SavePath = [BasePath,'StackObjects'];                           
                [obj.RoiMicrogliaPath,obj.RoiBrainPath] = PathsComplete('rois','roib');
            elseif exist('mode1','var') && mode1 == 9 %temporary mode with temp folder
                [obj.OriginalPath,BasePath] = PathsComplete('or','bp');            
                 obj.SavePath = [BasePath,'/','StackObjects'];
                 obj.OriginalPath = [BasePath,'/','originalTEMP'];               
                [obj.RoiMicrogliaPath,obj.RoiBrainPath] = PathsComplete('rois','roib');
            else %open new folder
                disp('select source folder')
                obj.OriginalPath = uigetdir;
                disp('select save folder')
                obj.SavePath = uigetdir;
            end
                        
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
            
            %{
                             
            if exist('SpotPath','var')
                obj.SpotPath = SpotPath;
            end
                        
            temp = load([obj.SpotPath,'/stackinfo.mat']);
            stackinfo = temp.stackinfo;
                      
            temp = load([obj.SpotPath,'/zfinput.mat']);
            obj.zfinput = temp.zfinput;
            
            obj.Numb = numel(stackinfo);

            %open and store roi annotated and computed info
            cmbr = cell(obj.Numb,1);
            ambr = cell(obj.Numb,1);
            SpotAnn = cell(obj.Numb,1);
            
            for k1 = 1:obj.Numb
                RoiBrain = ReadImageJROI([obj.RoiBrainPath,'/',stackinfo(k1).stackname,'.zip']);
                RoiMicroglia = ReadImageJROI([obj.RoiMicrogliaPath,'/',stackinfo(k1).stackname,'.roi']);

                tform_1234 = stackinfo(k1).Registration(strcmp({stackinfo(k1).Registration.name},'tform_complete')).value; 

                %computed midbrain roi
                cmbr{k1} = stackinfo(k1).BrainSegmentation.BrainEdge;
                
                %annotated midbrain roi
                amb = sng_roicell2poly(RoiBrain,1);   
                [ambr{k1}(:,1),ambr{k1}(:,2)] = transformPointsForward(tform_1234,amb(:,1),amb(:,2));
                ambr{k1} = double(ambr{k1});          
            
                %computed spots
                obj.SpotParameters{k1} =  stackinfo(k1).SpotSelection(strcmp({stackinfo(k1).SpotSelection.name},'SpotParameters')).value;

                %annotated spots
                SpotAnn{k1} = zeros(size(RoiMicroglia.mfCoordinates,1),2);
                [SpotAnn{k1}(:,1),SpotAnn{k1}(:,2)] = transformPointsForward(tform_1234,...
                    RoiMicroglia.mfCoordinates(:,1),RoiMicroglia.mfCoordinates(:,2));                
            end

            temp = {stackinfo.stackname};[obj.ImageInfo(1:obj.Numb).Name] = temp{:};
            [obj.BrainInfo(1:obj.Numb).ComputedMidBrain] = cmbr{:};
            [obj.BrainInfo(1:obj.Numb).AnnotatedMidBrainRegistrated] = ambr{:};           
            [obj.SpotInfo(1:obj.Numb).AnnotatedSpots] = SpotAnn{:};
            
            %}
                  
        end
        
        function obj = SliceCombination(obj,fishnumbers)
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
                elseif (max(fishnumbers) <= numel(imageinfotemp))
                    obj.ImageInfo = imageinfotemp(fishnumbers);    
                else 
                    error('at least one fish does not exist in StackInfo')
                end
                obj.fishnumbers.slice2stack = fishnumbers;
                    
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
                if ~exist('fishnumbers','var')
                    fishnumbers = 1:numel(obj.StackInfo);
                elseif (max(fishnumbers) <= numel(obj.StackInfo))
                    nfishes = numel(fishnumbers);                    
                else
                    error('at least one fish does not exist in StackInfo')
                end
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
                        [~,ppoutput] = PreprocessionLink(ImageSlice,obj.ZFParameters);
                    else
                        [~,ppoutput(k1)] = PreprocessionLink(ImageSlice,obj.ZFParameters);
                    end
                    
%TODO for imaging   [CorrectedSlice,ppoutput,ImageSliceCor,FiltIm] = PreprocessionLink(ImageSlice,zfinput)

%                    if savefig2_TF
%                        sng_SaveCell2TiffStack(CorrectedSlice,[obj.SavePath,'/',StackInfo(k1).stackname,'.tif'])
%                    end
                end
                                
                obj.PreprocessionInfo = ppoutput;
               
            end       
        end
        
        function obj = ExtendedDeptOfField(obj,fishnumbers)
            
             if ~exist('fishnumbers','var')
                    fishnumbers = 1:numel(obj.PreprocessionInfo);
             elseif (max(fishnumbers) <= numel(obj.PreprocessionInfo))
                     nfishes = numel(fishnumbers);
             else
                    error('at least one fish does not exist in PreProcession')                   
             end
             obj.fishnumbers.preprocession = fishnumbers;   

             %{
             for j = 2:numel(cellImg1)
                [cellImg1{j}, supportECC] = iat_inverse_warping(cellImg1{j}, ECCWarp{j}, par.transform, 1:N, 1:M);
             support = support & logical(supportECC);
             end
             %}
             
             
             
             
             
             
             
             %%%GA HIER VERDER
             
             [Icombined,edoutput{k2}] = ExtendedDeptofFieldLink2(CorrectedSlice,obj.ZFParameters);

             
             
        end
    end
end


%{        
        
        function obj = BrainVal(obj)
            
            Jaccard = zeros(obj.Numb,1);
            Dice = zeros(obj.Numb,1);
            
            for k1 = 1:obj.Numb

                mx = round(max([obj.BrainInfo(k1).ComputedMidBrain;...
                fliplr(obj.BrainInfo(k1).AnnotatedMidBrainRegistrated)]));
                
                a = poly2mask(obj.BrainInfo(k1).ComputedMidBrain(:,2),...
                    obj.BrainInfo(k1).ComputedMidBrain(:,1),mx(1),mx(2));
                b = poly2mask(obj.BrainInfo(k1).AnnotatedMidBrainRegistrated(:,1),...
                    obj.BrainInfo(k1).AnnotatedMidBrainRegistrated(:,2),mx(1),mx(2));
                
                intersection = (a & b);
                union = (a | b);

                Jaccard(k1) = sum(intersection(:))/sum(union(:));
                %Overlap(k1) = sum(intersection(:)) / min(sum(a(:)),sum(b(:)));
                Dice(k1) = (2*sum(intersection(:)))/(sum(union(:))+sum(intersection(:)));
                
            end
            
            %fill the BrainInfo Structure (per fish)
            temp = num2cell(Jaccard);[obj.BrainInfo(1:obj.Numb).Jaccard] = temp{:};
            temp = num2cell(Dice);[obj.BrainInfo(1:obj.Numb).Dice] = temp{:};
            
            %fill BrainStats Structure, statistical values
            obj.BrainStats.FiveNumberSummaryJaccard = sng_FiveNumberSum(Jaccard);
            obj.BrainStats.MeanJaccard = mean(Jaccard);
            
            obj.BrainStats.FiveNumberSummaryDice = sng_FiveNumberSum(Dice);
            obj.BrainStats.MeanDice = mean(Dice);       
            
        end
        
        function obj = SpotVal(obj)
            
            LinkDistance = cell(obj.Numb,1);
            CorrectSpots = cell(obj.Numb,1);
            nCorrect = zeros(obj.Numb,1);              
            FalsePosSpots = cell(obj.Numb,1);
            nFalsePos = zeros(obj.Numb,1);              
            FalseNegSpots = cell(obj.Numb,1);
            nFalseNeg = zeros(obj.Numb,1);
            Precision = zeros(obj.Numb,1);
            Recall = zeros(obj.Numb,1);
            F1score = zeros(obj.Numb,1);
            AbsDifference = zeros(obj.Numb,1);
            RelDifference = zeros(obj.Numb,1);
            SpotCom = cell(obj.Numb,1);
            ambsS = cell(obj.Numb,1);     
            
            
            %value added to prevent for substraction by zero        
            a = 0.001; 

            for k1 = 1:obj.Numb

                Spotpar = obj.SpotParameters{k1};
                ambr = obj.BrainInfo(k1).AnnotatedMidBrainRegistrated;              
                
                if ~isempty(Spotpar)
    %{
                    SpotparS = Spotpar(...
                        [Spotpar.MinProbability] &... %selection of spotinfo based on conditions
                        [Spotpar.Insite] &...
                        [Spotpar.LargerThan] &...
                        [Spotpar.SmallerThan]);
    %}
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
                ambs = obj.SpotInfo(k1).AnnotatedSpots;  
                [ins2,~] = inpolygon(ambs(:,1),ambs(:,2),ambr(:,1),ambr(:,2));
                ambsx = ambs(:,1);ambsx = ambsx(ins2);
                ambsy = ambs(:,2);ambsy = ambsy(ins2);                                                
                ambsS{k1} = [ambsx,ambsy];                
%{
                [Correct,FalsePos,FalseNeg,link] = sng_CoordinateMatching...
                    (SpotCom,obj.SpotInfo(k1).AnnotatedSpots,10);
%}
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
            [obj.SpotInfo(1:obj.Numb).LinkDistance] = LinkDistance{:};
            [obj.SpotInfo(1:obj.Numb).CorrectSpots] = CorrectSpots{:};
            temp = num2cell(nCorrect);[obj.SpotInfo(1:obj.Numb).nCorrect] = temp{:};
            [obj.SpotInfo(1:obj.Numb).FalsePosSpots] = FalsePosSpots{:};
            temp = num2cell(nFalsePos); [obj.SpotInfo(1:obj.Numb).nFalsePos] = temp{:};
            [obj.SpotInfo(1:obj.Numb).FalseNegSpots] = FalseNegSpots{:};
            temp = num2cell(nFalseNeg); [obj.SpotInfo(1:obj.Numb).nFalseNeg] = temp{:};
            temp = num2cell(Precision); [obj.SpotInfo(1:obj.Numb).Precision] = temp{:};
            temp = num2cell(Recall); [obj.SpotInfo(1:obj.Numb).Recall] = temp{:};
            temp = num2cell(F1score); [obj.SpotInfo(1:obj.Numb).F1score] = temp{:};            
            temp = num2cell(AbsDifference); [obj.SpotInfo(1:obj.Numb).AbsDifference] = temp{:};
            temp = num2cell(RelDifference); [obj.SpotInfo(1:obj.Numb).RelDifference] = temp{:};
            [obj.SpotSelection(1:obj.Numb).AnnotatedSpots] = ambsS{:};                
            [obj.SpotSelection(1:obj.Numb).ComputedSpots] = SpotCom{:};
                       
        end
       
        function obj = SpotBrainVal(obj)
            
            LinkDistance = cell(obj.Numb,1);
            CorrectSpots = cell(obj.Numb,1);
            nCorrect = zeros(obj.Numb,1);              
            FalsePosSpots = cell(obj.Numb,1);
            nFalsePos = zeros(obj.Numb,1);              
            FalseNegSpots = cell(obj.Numb,1);
            nFalseNeg = zeros(obj.Numb,1);
            Precision = zeros(obj.Numb,1);
            Recall = zeros(obj.Numb,1);
            F1score = zeros(obj.Numb,1);
            AbsDifference = zeros(obj.Numb,1);
            RelDifference = zeros(obj.Numb,1);            

            %value added to prevent for substraction by zero        
            a = 0.001; 

            for k1 = 1:obj.Numb

                Spotpar = obj.SpotParameters{k1};
                ambr = obj.BrainInfo(k1).AnnotatedMidBrainRegistrated;              

                SpotparS = Spotpar(...
                    [Spotpar.MinProbability] &... %selection of spotinfo based on conditions
                    [Spotpar.Insite] &...
                    [Spotpar.LargerThan] &...
                    [Spotpar.SmallerThan]);

                SpotCom = reshape([SpotparS.Centroid],2,numel(SpotparS))';
                
                %because the annotated midbrain does not corresponds
                %exactly with the annotated spots, de annotated spots are
                %also filtered to compute the performance of the
                %spotdetection only
                ambs = obj.SpotInfo(k1).AnnotatedSpots;  
                [ins2,~] = inpolygon(ambs(:,1),ambs(:,2),ambr(:,1),ambr(:,2));
                ambsx = ambs(:,1);ambsx = ambsx(ins2);
                ambsy = ambs(:,2);ambsy = ambsy(ins2);                                                
                ambsS{k1} = [ambsx,ambsy];                
%{
                [Correct,FalsePos,FalseNeg,link] = sng_CoordinateMatching...
                    (SpotCom,obj.SpotInfo(k1).AnnotatedSpots,10);
%}
                [Correct,FalsePos,FalseNeg,link] = sng_CoordinateMatching...
                    (SpotCom,ambsS{k1},10);                   
                
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
            [obj.SpotBrainInfo(1:obj.Numb).LinkDistance] = LinkDistance{:};
            [obj.SpotBrainInfo(1:obj.Numb).CorrectSpots] = CorrectSpots{:};
            temp = num2cell(nCorrect);[obj.SpotBrainInfo(1:obj.Numb).nCorrect] = temp{:};
            [obj.SpotBrainInfo(1:obj.Numb).FalsePosSpots] = FalsePosSpots{:};
            temp = num2cell(nFalsePos); [obj.SpotBrainInfo(1:obj.Numb).nFalsePos] = temp{:};
            [obj.SpotBrainInfo(1:obj.Numb).FalseNegSpots] = FalseNegSpots{:};
            temp = num2cell(nFalseNeg); [obj.SpotBrainInfo(1:obj.Numb).nFalseNeg] = temp{:};
            temp = num2cell(Precision); [obj.SpotBrainInfo(1:obj.Numb).Precision] = temp{:};
            temp = num2cell(Recall); [obj.SpotBrainInfo(1:obj.Numb).Recall] = temp{:};
            temp = num2cell(F1score); [obj.SpotBrainInfo(1:obj.Numb).F1score] = temp{:};
            temp = num2cell(AbsDifference); [obj.SpotBrainInfo(1:obj.Numb).AbsDifference] = temp{:};
            temp = num2cell(RelDifference); [obj.SpotBrainInfo(1:obj.Numb).RelDifference] = temp{:};                       
        end
        
        function show(obj,exportit)
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
        
        function showbrain(obj,val,ComOrAnn)
            
            %CompOrAnn determines if the computed or the annotated brain is
            %displayed inlcudin ghte Correct/FalsePos/FalseNeg
                
            
            
            TemplatePath3dpf = '/Users/samuelgeurts/Dropbox/Template 3 dpf';
            PreprocessionPath = '/Users/samuelgeurts/Dropbox/20170327_3dpf_test/1_preprocessed';
            CompleteTemplate = LoadTemplateLink3(TemplatePath3dpf);
            
            temp = load([obj.SpotPath,'/stackinfo.mat']);
            stackinfo = temp.stackinfo;    
            
            for k = 1:numel(val)
               
                k1 = val(k);
                tform_1234 = stackinfo(k1).Registration(strcmp({stackinfo(k1).Registration.name},'tform_complete')).value; 

                CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',[obj.ImageInfo(k1).Name],'.tif']);           
                Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k1).ExtendedDeptOfField.IndexMatrix);
                Ialigned = imwarp(Icombined,tform_1234,'FillValues',255,'OutputView',CompleteTemplate.ref_temp);
                                
                cmbr = obj.BrainInfo(k1).ComputedMidBrain;
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
                
        function showgraph(obj)
            fsx = 8.5;
            fsy = 8.5;

            h1 = figure('PaperUnits','centimeters','Color',[1 1 1]);
            sng_figcm(fsx,fsy);
            plot(sort([obj.BrainInfo.Jaccard],'descend'))
            title('brain  validation')
            xlabel('fish number','FontSize',8,'FontName','arial');
            ylabel('Jaccard index','FontSize',8,'FontName','arial');
            set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
            set(gca,'XLim',[1 50])
            set(gca,'YLim',[0.48 1.02])
            sng_figcm(fsx,fsy,113.6)  
            set(0, 'currentfigure', h1)
            %export_fig(h1 ,['/Users/samuelgeurts/Desktop/','graph',num2str(1)], '-png', '-r600', '-nocrop');
            

            h2 = figure('PaperUnits','centimeters','Color',[1 1 1]);
            sng_figcm(fsx,fsy);
            plot(sort([obj.BrainInfo.Dice],'descend'))
            title('brain validation')                                    
            xlabel('fish number','FontSize',8,'FontName','arial');
            ylabel('Dice index','FontSize',8,'FontName','arial');
            set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
            set(gca,'XLim',[1 50])
            set(gca,'YLim',[-0.02 1.02])
            sng_figcm(fsx,fsy,113.6)  
            set(0, 'currentfigure', h2)
            %export_fig(h2 ,['/Users/samuelgeurts/Desktop/','graph',num2str(2)], '-png', '-r600', '-nocrop');
            
            h3 = figure('PaperUnits','centimeters','Color',[1 1 1]);
            sng_figcm(fsx,fsy);
            plot(sort([obj.SpotInfo.F1score],'descend'))
            title('spot validation')      
            xlabel('fish number','FontSize',8,'FontName','arial');
            ylabel('f1score','FontSize',8,'FontName','arial');
            set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
            set(gca,'XLim',[1 50])
            set(gca,'YLim',[-0.02 1.02])
            sng_figcm(fsx,fsy,113.6)  
            set(0, 'currentfigure', h3)
            %export_fig(h3 ,['/Users/samuelgeurts/Desktop/','graph',num2str(3)], '-png', '-r600', '-nocrop');
            
            h4 = figure('PaperUnits','centimeters','Color',[1 1 1]);
            sng_figcm(fsx,fsy);            
            plot(sort([obj.SpotInfo.Recall],'descend'))
            title('spot validation')     
            xlabel('fish number','FontSize',8,'FontName','arial');
            ylabel('Recall','FontSize',8,'FontName','arial');
            set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
            set(gca,'XLim',[1 50])
            set(gca,'YLim',[-0.02 1.02])
            sng_figcm(fsx,fsy,113.6)  
            set(0, 'currentfigure', h4)
            %export_fig(h4 ,['/Users/samuelgeurts/Desktop/','graph',num2str(4)], '-png', '-r600', '-nocrop');
            
            h5 = figure('PaperUnits','centimeters','Color',[1 1 1]);
            sng_figcm(fsx,fsy);           
            plot(sort([obj.SpotInfo.Precision],'descend'))
            title('spot validation')
            xlabel('fish number','FontSize',8,'FontName','arial');
            ylabel('Precision','FontSize',8,'FontName','arial');
            set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
            set(gca,'XLim',[1 50])
            set(gca,'YLim',[-0.02 1.02])
            sng_figcm(fsx,fsy,113.6)  
            set(0, 'currentfigure', h5)
            %export_fig(h5 ,['/Users/samuelgeurts/Desktop/','graph',num2str(5)], '-png', '-r600', '-nocrop');            
            
            h6 = figure('PaperUnits','centimeters','Color',[1 1 1]);
            sng_figcm(fsx,fsy)
            plot(sort([obj.SpotBrainInfo.F1score],'descend'))
            title('spot with computed brain')      
            xlabel('fish number','FontSize',8,'FontName','arial');
            ylabel('f1score','FontSize',8,'FontName','arial');
            set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
            set(gca,'XLim',[1 50])
            set(gca,'YLim',[-0.02 1.02])
            sng_figcm(fsx,fsy,113.6)  
            set(0, 'currentfigure', h6)
            %export_fig(h6 ,['/Users/samuelgeurts/Desktop/','graph',num2str(6)], '-png', '-r600', '-nocrop');
            
            h7 = figure('PaperUnits','centimeters','Color',[1 1 1]);
            sng_figcm(fsx,fsy);
            plot(sort([obj.SpotBrainInfo.Recall],'descend'))
            title('spot with computed brain')     
            xlabel('fish number','FontSize',8,'FontName','arial');
            ylabel('Recall','FontSize',8,'FontName','arial');
            set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
            set(gca,'XLim',[1 50])
            set(gca,'YLim',[-0.02 1.02])
            sng_figcm(fsx,fsy,113.6)  
            set(0, 'currentfigure', h7)
            %export_fig(h7 ,['/Users/samuelgeurts/Desktop/','graph',num2str(7)], '-png', '-r600', '-nocrop');
            
            h8 = figure('PaperUnits','centimeters','Color',[1 1 1]);
            sng_figcm(fsx,fsy);
            plot(sort([obj.SpotBrainInfo.Precision],'descend'))
            title('spot with computed brain') 
            xlabel('fish number','FontSize',8,'FontName','arial');
            ylabel('Precision','FontSize',8,'FontName','arial');
            set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
            set(gca,'XLim',[1 50])
            set(gca,'YLim',[-0.02 1.02])
            sng_figcm(fsx,fsy,113.6)  
            set(0, 'currentfigure', h8)
            %export_fig(h8,['/Users/samuelgeurts/Desktop/','graph',num2str(8)], '-png', '-r600', '-nocrop');
            
        end
                    
        function saveit(obj,savename)
            obj.date = datetime;
            if ~exist('savename','var')
                savename = inputname(1);
            end
            save([obj.SavePath,'/',savename,'.mat'],'obj');          
        end
                        
    end
        
end
%}   