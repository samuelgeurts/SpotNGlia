classdef SpotNGlia
    
    properties
        FishPath = []
        SourcePath = []
        SavePath = []
        
        AnnotatedMidBrainPath = []
        AnnotatedSpotPath = []
        
        SaveName = []
        InfoName = []
        
        Version = 'SpotNGlia 1.5.3'
        User = []
        OS = []
        Delimiter = ';'
        
        nfishes = []
        fishnumbers = []
        slicenumbers = []
        savedate = []
        
        ZFParameters = []
        
        StackInfo = []
        ImageInfo = []
        BatchInfo = []
        checkup = []
        
        InfoResult = []
        
        Computations = []
        Annotations = []
        
        BrainInfo = []
        BrainStats = []
        SpotInfo = []
        SpotStats = []
        SpotBrainInfo = []
        SpotBrainStats = []
        
        ImageInfoChecked_TF = false
        Sorting = []
    end
    
    properties(Transient = true)
        CompleteTemplate = []
    end
    
    methods
        %constructor function
        function obj = SpotNGlia(mode1)
            %obj.NewPath(0) default mode, UI for FishPath and SavePath, source is assumed to be in SpotNGlia folder
            %obj.NewPath(1) UI for FishPath, SavePath and SourcePath
            %mode 10,11,12 only works for S Geurts own pc's
            %obj.NewPath(10) this mode assumes the trainingset and set known paths
            %obj.NewPath(11) this mode assumes that Fishpath and SavePath are equal and user UI to select new path
            %obj.NewPath(12) switch between own pc's for MultiBatch on Seagate harddisk
            
            if ~exist('mode1', 'var')
                mode1 = 0;
            end
            
            if ~isempty(mode1)
                
                obj = NewPath(obj, mode1);
                
                [~, folder] = fileparts(obj.FishPath);
                obj.SaveName = strcat('SNG_', folder);
                obj.InfoName = strcat('INFO_', folder);
                
                load([obj.SourcePath, '/', 'zfinput.mat']); %#ok<LOAD>
                obj.ZFParameters = zfinput;
                obj.CompleteTemplate = load([obj.SourcePath, '/', 'Template3dpf', '.mat']);
                
                %file which contains all batch specific computation info
                matObj = matfile([obj.SavePath, '/', obj.InfoName, '.mat']);
                obj.InfoResult = whos(matObj);
                
            end
        end
        function obj = NewPath(obj, mode1)
            %obj.NewPath(0) default mode, UI for FishPath and SavePath, source is assumed to be in SpotNGlia folder
            %obj.NewPath(1) UI for FishPath, SavePath and SourcePath
            %mode 10,11,12 only works for S Geurts own pc's
            %obj.NewPath(10) this mode assumes the trainingset and set known paths
            %obj.NewPath(11) this mode assumes that Fishpath and SavePath are equal and user UI to select new path
            %obj.NewPath(12) switch between own pc's for MultiBatch on Seagate harddisk
            
            %msg1 = ['select image folder to process'];
            %msg2 = ['select destination folder to save'];
            %msg3 = ['select program source folder'];
            
            msg1 = ['select image folder to process:  ', obj.FishPath];
            msg2 = ['select destination folder to save:  ', obj.SavePath];
            msg3 = ['select program source folder:  ', obj.SourcePath];
            
            if mode1 == 0
                disp(msg1)
                obj.FishPath = uigetdir([], msg1);
                disp(msg2)
                obj.SavePath = uigetdir(obj.FishPath, msg2);
                
                % find path where SpotNGlia.m is saved and reconstruct source path
                SpotNGliaPath = mfilename('fullpath');
                %expression = ['\', filesep];
                splitStr = regexp(SpotNGliaPath, filesep, 'split');
                path = [];
                for k1 = 1:numel(splitStr) - 3
                    path = [path, splitStr{k1}, filesep]; %#ok<AGROW>
                end
                obj.SourcePath = [path, 'Source'];
                
                %obj.OS = getenv('OS');
                if isempty(obj.OS)
                    obj.User = getenv('USER');
                else
                    obj.User = getenv('username');
                end
            elseif mode1 == 1
                disp(msg1)
                obj.FishPath = uigetdir([], msg1);
                disp(msg2)
                obj.SavePath = uigetdir([], msg2);
                disp(msg3)
                obj.SourcePath = uigetdir([], msg3);
                
                %obj.OS = getenv('OS');
                if isempty(obj.OS)
                    obj.User = getenv('USER');
                else
                    obj.User = getenv('username');
                end
            elseif (mode1 == 10) || (mode1 == 11) || (mode1 == 12)
                if strcmp(getenv('OS'), 'Windows_NT') && strcmp(getenv('username'), '260018')
                    BasePath = 'C:\Users\260018\Dropbox';
                    User = '260018';
                elseif strcmp(getenv('OS'), 'Windows_NT') && strcmp(getenv('username'), 'SNGeu')
                    BasePath = 'C:\Users\SNGeu\Dropbox';
                    User = 'SNGeu';
                elseif strcmp(getenv('USER'), 'samuelgeurts')
                    BasePath = '/Users/samuelgeurts/Dropbox';
                    User = 'samuelgeurts';
                elseif strcmp(getenv('USER'), 'Anouk')
                    BasePath = '/Users/Anouk/Dropbox';
                    User = 'Anouk';
                else
                    error('unknown platform and user');
                end
                
                %this modes 10,11,12 only work for my own pc and are used for development (samuel geurts)
                if mode1 == 10
                    %this mode assumes the trainingset and set known paths
                    obj.FishPath = [BasePath, '/', '20170327_3dpf'];
                    obj.SavePath = [BasePath, '/', 'SpotNGlia Destination'];
                    obj.SourcePath = [BasePath, '/', 'SpotNGlia Source'];
                    obj.AnnotatedMidBrainPath = [BasePath, '/', 'SpotNGlia Source', '/', 'Roi brain'];
                    obj.AnnotatedSpotPath = [BasePath, '/', 'SpotNGlia Source', '/', 'Roi microglia'];
                elseif mode1 == 11
                    %this mode assumes that Fishpath and SavePath are equal and user UI to select new path
                    %
                    disp(msg1)
                    obj.FishPath = uigetdir([], msg1);
                    obj.SavePath = obj.FishPath;
                    obj.SourcePath = [BasePath, '/', 'SpotNGlia Source'];
                elseif mode1 == 12
                    %this mode automatically changes paths between users '260018' and 'samuelgeurts' an SNGeu
                    %i.e. between macbookpro and erasmus pc and medion pc.
                    %Only when fishes are on harddisk
                    
                    if strcmp(User, 'SNGeu') && strcmp(obj.User, 'samuelgeurts')
                        splitStr = regexp(obj.FishPath, '\/', 'split');
                        if strcmp(splitStr{2}, 'Volumes')
                            filename = 'D:';
                            for k2 = 4:numel(splitStr)
                                filename = [filename, '/', splitStr{k2}]; %#ok<AGROW>
                            end
                            obj.FishPath = filename;
                        end
                    end
                    
                    if strcmp(User, 'samuelgeurts') && strcmp(obj.User, 'SNGeu')
                        expression = '\';
                        splitStr = regexp(obj.FishPath, expression, 'split');
                        if strcmp(splitStr{1}, 'D:')
                            filename = '/Volumes/Seagate Expansion Drive';
                            for k2 = 2:numel(splitStr)
                                filename = [filename, '/', splitStr{k2}]; %#ok<AGROW>
                            end
                            obj.FishPath = filename;
                        end
                    end
                    
                    if strcmp(User, '260018') && strcmp(obj.User, 'samuelgeurts')
                        expression = '\/';
                        splitStr = regexp(obj.FishPath, expression, 'split');
                        
                        if strcmp(splitStr{1}, 'Volumes')
                            filename = 'H:';
                            for k2 = 3:numel(splitStr)
                                filename = [filename, '/', splitStr{k2}]; %#ok<AGROW>
                            end
                            obj.FishPath = filename;
                        end
                    end
                    
                    if strcmp(User, 'samuelgeurts') && strcmp(obj.User, '260018')
                        expression = '\\';
                        splitStr = regexp(obj.FishPath, expression, 'split');
                        
                        if strcmp(splitStr{1}, 'H:')
                            filename = '/Volumes/Seagate Expansion Drive';
                            for k2 = 2:numel(splitStr)
                                filename = [filename, '/', splitStr{k2}]; %#ok<AGROW>
                            end
                            obj.FishPath = filename;
                        end
                    end
                    obj.SavePath = obj.FishPath;
                    obj.SourcePath = [BasePath, '/', 'SpotNGlia Source'];
                end
                obj.User = User;
            else
                error('unknown mode')
            end
            
            obj.OS = getenv('OS');
            
        end
        function obj = SliceCombination(obj, slicenumbers)
            %computes imageinfo and stackinfo
            %       sorting = ['date'/'name']
            %       CheckImageInfo_TF [true/false]
            %Example:   obj = SliceCombination(obj)
            %Example:   obj = SliceCombination(obj,'date',true)
            %Example:   obj = SliceCombination(obj,'name',true)
            
            if ~isempty(obj.ImageInfo)
                answer = questdlg('Are you sure to overwrite ImageInfo and StackInfo', '', 'Ok', 'Cancel', 'Cancel');
                if strcmp(answer, 'Ok')
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
                obj.Sorting = questdlg('Which sorting algorithm do you want?', 'Choose sorting method', 'Date', 'Name', 'Date');
                if strcmp(obj.Sorting, 'Date')
                    obj.ZFParameters = sng_zfinput(obj.ZFParameters, 0, 'imageinfo', 'stackselection', 'sorting', 'Date', 'high'); %date or name
                elseif strcmp(obj.Sorting, 'Name')
                    obj.ZFParameters = sng_zfinput(obj.ZFParameters, 0, 'imageinfo', 'stackselection', 'sorting', 'Name', 'high'); %date or name
                end
            end
            
            if isempty(obj.ImageInfo)
                % compute ImageInfo and StackInfo
                dirinfo = dir([obj.FishPath, '/*.', 'tif']);
                
                %removes filenames from folder which starts with '._'
                %an issue induces sometimes by extern harddisk and macs
                ind = regexp({dirinfo.name}, '._');
                for k1 = 1:numel(ind)
                    if ~isempty(ind{k1}) && ind{k1}(1) == 1
                        TF(k1) = false;
                    else
                        TF(k1) = true;
                    end
                end
                dirinfo = dirinfo(TF);
                
                %removes unimportant fields
                imageinfotemp = rmfield(dirinfo, {'isdir', 'datenum', 'bytes'});
                
                if ~isfield(imageinfotemp, 'folder')
                    [imageinfotemp.folder] = deal(obj.FishPath);
                end
                
                % make a selection of fishslices if fishslices is given as input
                if ~exist('slicenumbers', 'var')
                    obj.ImageInfo = imageinfotemp;
                    slicenumbers = 1:numel(imageinfotemp);
                elseif (max(slicenumbers) <= numel(imageinfotemp))
                    obj.ImageInfo = imageinfotemp(slicenumbers);
                else
                    error('at least one fish does not exist in StackInfo')
                end
                obj.slicenumbers.slicecombination = slicenumbers;
                
                [obj.ImageInfo] = ImageInfoSNG(obj.ImageInfo, obj.ZFParameters);
                [obj.StackInfo] = StackInfoSNG(obj.ImageInfo);
                
                ImageInfo = obj.ImageInfo; %#ok<NASGU>
                StackInfo = obj.StackInfo; %#ok<NASGU>
                
                save([obj.SavePath, '/', obj.InfoName, '.mat'], 'ImageInfo')
                save([obj.SavePath, '/', obj.InfoName, '.mat'], 'StackInfo', '-append')
            end
            
            if ~obj.ImageInfoChecked_TF
                
                %openvar('ImageInfo')
                %openvar('StackInfo')
                
                str = sprintf(['Check slice combination in ImageInfo and StackInfo.', ...
                    '\nApply corrections in "imageinfo.CorNextStack".', ...
                    '\n   Set value "1" for new fish', ...
                    '\n   Set value "2" for removing image from stack.', ...
                    '\n   Set value "0" if slice belongs to previous slice']);
                msgbox(str);
            end
            obj.saveit
            obj.nfishes = numel(obj.StackInfo);
            
        end
        function obj = PreProcession(obj, fishnumbers)
            
            h = waitbar(0, 'Preprocession', 'Name', 'SpotNGlia');
            [obj.StackInfo] = StackInfoSNG(obj.ImageInfo);
            obj.nfishes = numel(obj.StackInfo);
            
            %if ~isempty(obj.ImageInfo)
            %    [obj.StackInfo] = StackInfoSNG(obj.ImageInfo);
            %else
            %    error('First run SliceCombination');
            %end
            
            if ~obj.ImageInfoChecked_TF
                answer = questdlg('Do you confirm ImageInfo?', '', 'Yes', 'No', 'No');
                if strcmp(answer, 'No')
                    msgbox('First run SliceCombination');
                elseif strcmp(answer, 'Yes')
                    obj.ImageInfoChecked_TF = true;
                end
            end
            
            if obj.ImageInfoChecked_TF
                obj.slicenumbers.stackinfo = find([obj.ImageInfo.CorNextStack] ~= 2);
                
                if ~exist('fishnumbers', 'var')
                    fishnumbers = 1:numel(obj.StackInfo);
                elseif (max(fishnumbers) > numel(obj.StackInfo))
                    error('at least one fish does not exist in StackInfo')
                end
                nfishes = numel(fishnumbers); %#ok<*PROPLC>
                obj.fishnumbers.preprocession = fishnumbers;
                
                %only preallocate if more than 1 fish has to be processed
                if nfishes > 1
                    PreprocessionInfo(nfishes) = struct('ColorWarp', [], ...
                        'SliceWarp', [], ...
                        'NNC_Img_before', [], ...
                        'NNC_Img_after', [], ...
                        'NNC_BP_before', [], ...
                        'NNC_BP_after', []);
                end
                
                %folder to save CorrectedFish
                TempFolderName = ([obj.SavePath, '/', 'CorrectedFish']);
                if ~exist(TempFolderName, 'dir')
                    mkdir(TempFolderName)
                end
                
                for k1 = 1:nfishes
                    
                    waitbar(k1/nfishes, h, ['Preprocession ', num2str(k1), ' / ', num2str(nfishes)])
                    
                    %select slices
                    fn = fishnumbers(k1);
                    
                    ImageSlice = cell(1, obj.StackInfo(fn).stacksize); %preallocate for every new slice
                    for k2 = 1:obj.StackInfo(fn).stacksize
                        ImageSlice{k2} = imread([obj.FishPath, '/', obj.StackInfo(fn).imagenames{k2}]);
                        ImageSlice{k2} = im2uint8(ImageSlice{k2});
                    end
                    
                    %                    if savefig1_TF
                    %                        sng_SaveCell2TiffStack(ImageSlice,[SavePath,'/stack/',stackinfo(k1).stackname,'.tif'])
                    %                    end
                    
                    %[CorrectedSlice,ppoutput(k1)] = PreprocessionSNG(ImageSlice,obj.ZFParameters);
                    if k1 == 1
                        [CorrectedFishTemp, PreprocessionInfo] = PreprocessionSNG(ImageSlice, obj.ZFParameters);
                    else
                        [CorrectedFishTemp, PreprocessionInfo(k1)] = PreprocessionSNG(ImageSlice, obj.ZFParameters);
                    end
                    
                    sng_SaveCell2TiffStack(CorrectedFishTemp, [TempFolderName, '/', obj.StackInfo(k1).stackname, '.tif'])
                    
                    %obj.CorrectedFish(k1).image = CorrectedFishTemp; <- %save to object takes to much space
                    
                    %TODO for imaging   [CorrectedSlice,ppoutput,ImageSliceCor,FiltIm] = PreprocessionSNG(ImageSlice,zfinput)
                    %                    if savefig2_TF
                    %                        sng_SaveCell2TiffStack(CorrectedSlice,[obj.SavePath,'/',StackInfo(k1).stackname,'.tif'])
                    %                    end
                end
            end
            obj.saveit
            delete(h)
            save([obj.SavePath, '/', obj.InfoName, '.mat'], 'PreprocessionInfo', '-append')
            
        end
        function obj = ExtendedDeptOfField(obj, fishnumbers)
            
            h = waitbar(0, 'Extended Dept of Field', 'Name', 'SpotNGlia');
            
            if ~exist('fishnumbers', 'var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif max(fishnumbers) > numel(obj.StackInfo)
                error('at least one fish does not exist in PreProcessionInfo')
            end
            nfishes = numel(fishnumbers);
            
            obj.fishnumbers.extendeddeptoffield = fishnumbers;
            
            
            %{
                              %Stack has to be added but than input parameters are needed,
                              so the function PreprocessionSNG has to be separated in that
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
                ExtendedDeptOfFieldInfo(nfishes) = struct('IndexMatrix', [], ...
                    'variance_sq', []);
            end
            
            %folder to save CombinedFish
            TempFolderName = ([obj.SavePath, '/', 'CombinedFish']);
            if ~exist(TempFolderName, 'dir')
                mkdir(TempFolderName)
            end
            
            for k1 = 1:nfishes
                waitbar(k1/nfishes, h, ['Extended Dept of Field ', num2str(k1), ' / ', num2str(nfishes)])
                fn = fishnumbers(k1);
                
                CorrectedFish = sng_openimstack2([obj.SavePath, '/', 'CorrectedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
                [CombinedFishTemp, ExtendedDeptOfFieldInfo(k1)] = ExtendedDeptofFieldSNG(CorrectedFish, obj.ZFParameters);
                %obj.CombinedFish(k1).image = CombinedFishTemp; because it takes to much space
                imwrite(uint8(CombinedFishTemp), [TempFolderName, '/', obj.StackInfo(k1).stackname, '.tif'], ...
                    'WriteMode', 'overwrite', 'Compression', 'none');
            end
            
            obj.saveit
            save([obj.SavePath, '/', obj.InfoName, '.mat'], 'ExtendedDeptOfFieldInfo', '-append')
            delete(h)
        end
        function obj = Registration(obj, fishnumbers)
            
            h = waitbar(0, 'Registration', 'Name', 'SpotNGlia');
            
            if ~exist('fishnumbers', 'var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif max(fishnumbers) > numel(obj.StackInfo)
                error('at least one fish does not exist in ExtendedDeptOfFieldInfo')
            end
            nfishes = numel(fishnumbers);
            obj.fishnumbers.registration = fishnumbers;
            
            %load complete template
            obj = obj.LoadTemplate;
            
            
            %folder to save AlignedFish
            TempFolderName = ([obj.SavePath, '/', 'AlignedFish']);
            if ~exist(TempFolderName, 'dir')
                mkdir(TempFolderName)
            end
            %preallocation
            MaxFishColor = zeros(nfishes, 3);
            MeanFishColor = zeros(nfishes, 3);
            RegistrationInfo = cell(nfishes, 1);
            
            for k1 = 1:nfishes
                waitbar(k1/nfishes, h, ['Registration ', num2str(k1), ' / ', num2str(nfishes)])
                fn = fishnumbers(k1);
                CombinedFish = imread([obj.SavePath, '/', 'CombinedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
                [AlignedFishTemp, RegistrationInfo{k1, 1}] = AlignmentSNG(CombinedFish, obj.CompleteTemplate, obj.ZFParameters);
                %obj.AlignedFish(k1).image = AlignedFishTemp; %because it
                %takes to much space
                imwrite(uint8(AlignedFishTemp), [TempFolderName, '/', obj.StackInfo(k1).stackname, '.tif'], ...
                    'WriteMode', 'overwrite', 'Compression', 'none');
                %generates matrix with MeanFishColors
                MeanFishColor(k1, 1:3) = RegistrationInfo{k1}(strcmp({RegistrationInfo{k1}.name}, 'MeanFishColor')).value;
                MaxFishColor(k1, 1:3) = RegistrationInfo{k1}(strcmp({RegistrationInfo{k1}.name}, 'MaxFishColor')).value;
                
            end
            obj.BatchInfo.MeanMaxFishColor = mean(MaxFishColor(:));
            obj.BatchInfo.StdMaxFishColor = std(MaxFishColor(:));
            obj.BatchInfo.MeanFishColor = mean(MeanFishColor(:));
            obj.BatchInfo.StdFishColor = std(MeanFishColor(:));
            
            
            obj.saveit
            save([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo', '-append')
            delete(h)
        end
        function obj = BrainSegmentation(obj, fishnumbers)
            h = waitbar(0, 'Brain Segmentation', 'Name', 'SpotNGlia');
            
            if ~exist('fishnumbers', 'var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif max(fishnumbers) > numel(obj.StackInfo)
                error('at least one fish does not exist in RegistrationInfo')
            end
            nfishes = numel(fishnumbers);
            obj.fishnumbers.brainsegmentation = fishnumbers;
            
            %load complete template
            obj = obj.LoadTemplate;
            
            
            if nfishes > 1
                BrainSegmentationInfo(nfishes) = struct('EdgeFilterWidth', [], ...
                    'ShortestPath', [], ...
                    'ShortestPathValue', [], ...
                    'BrainEdge', [], ...
                    'PolarTransform', []);
            end
            
            for k1 = 1:nfishes
                fn = fishnumbers(k1);
                waitbar(k1/nfishes, h, ['Brain Segmentation ', num2str(k1), ' / ', num2str(nfishes)])
                AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
                [~, BrainSegmentationInfo(k1)] = MidBrainDetectionSNG(AlignedFish, obj.CompleteTemplate, obj.ZFParameters);
                %obj.BrainFish(k1).image = BrainFishTemp;
            end
            obj.saveit
            save([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo', '-append')
            
            delete(h)
        end
        function obj = SpotDetection(obj, fishnumbers)
            
            obj = obj.LoadTemplate;
            load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo')
            
            %load checkup if exist from pre SpotNGlia1.4.0
            obj = FillCheckup(obj);
            
            h = waitbar(0, 'SpotDetection', 'Name', 'SpotNGlia');
            
            if ~exist('fishnumbers', 'var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif max(fishnumbers) > numel(obj.StackInfo)
                error('at least one fish does not exist in BrainInfo')
            end
            nfishes = numel(fishnumbers);
            obj.fishnumbers.spotdetection = fishnumbers;
            
            %load complete template
            
            %if nfishes > 1
            %	spoutput(nfishes) = struct('EdgeFilterWidth',[],...
            %        'ShortestPath',[],...
            %        'ShortestPathValue',[],...
            %        'BrainEdge',[]);
            %end
            SpotsDetected = cell(nfishes, 1);
            SpotParameters = cell(nfishes, 1);
            
            
            if nfishes > 1
                SpotDetectionInfo(nfishes) = struct('Multiproduct', [], ...
                    'MultiproductThreshold', []);
            end
            
            for k1 = 1:nfishes
                fn = fishnumbers(k1);
                waitbar(k1/nfishes, h, ['SpotDetection ', num2str(k1), ' / ', num2str(nfishes)])
                AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
                [SpotsDetected{k1}, SpotParameters{k1}, SpotDetectionInfo(k1)] = SpotDetectionSNG(AlignedFish, obj.CompleteTemplate, BrainSegmentationInfo(k1).BrainEdge, obj.ZFParameters);
                
                %update Computations
                obj.Computations(k1).Spots = reshape([SpotsDetected{k1}.Centroid], 2, numel(SpotsDetected{k1}))';
                obj.Computations(k1).Midbrain = fliplr(BrainSegmentationInfo(k1).BrainEdge);
                obj.Computations(k1).Counts = numel(SpotsDetected{k1});
                
                %update checkup after new spot computations
                if ~isempty(obj.checkup)
                    obj.checkup(k1).Spots = obj.SpotsInsiteArea(SpotParameters{k1}, obj.checkup(k1).Midbrain);
                    obj.checkup(k1).Counts = size(obj.checkup(k1).Spots, 1);
                end
            end
            
            obj.saveit
            save([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotDetectionInfo', '-append')
            save([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected', '-append')
            save([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters', '-append')
            
            obj.buildsheet
            
            delete(h)
        end
        function obj = SpotDetection2(obj, fishnumbers)
            
            INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo', 'BrainSegmentationInfo');
            TEMPLATE = load([obj.SourcePath, '/', 'Template3dpf.mat'], 'ref_temp', 'SVAP_index', 'SpotVectorArrayProbability');
            
            %load checkup if exist from pre SpotNGlia1.4.0
            obj = FillCheckup(obj);
            
            h = waitbar(0, 'SpotDetection', 'Name', 'SpotNGlia');
            
            if ~exist('fishnumbers', 'var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif max(fishnumbers) > numel(obj.StackInfo)
                error('at least one fish does not exist in BrainInfo')
            end
            
            nfishes = numel(fishnumbers);
            
            SpotsDetected = cell(nfishes, 1);
            SpotParameters = cell(nfishes, 1);
            if nfishes > 1
                SpotDetectionInfo(nfishes) = struct('Multiproduct', [], ...
                    'MultiproductThreshold', []);
            end
            
            for k1 = 1:nfishes
                fn = fishnumbers(k1);
                waitbar(k1/nfishes, h, ['SpotDetection on slice ', num2str(k1), ' / ', num2str(nfishes)])
                tform_complete = INFO.RegistrationInfo{fn}(strcmp({INFO.RegistrationInfo{fn}.name}, 'tform_complete')).value;
                CorrectedFish = sng_openimstack2([obj.SavePath, '/', 'CorrectedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
                AlignedSlice = cell(1, numel(CorrectedFish));
                for k2 = 1:numel(CorrectedFish)
                    AlignedSlice{k2} = imwarp(CorrectedFish{k2}, tform_complete, 'FillValues', 255, 'OutputView', TEMPLATE.ref_temp);
                end
                [SpotsDetected{fn}, SpotParameters{fn}, SpotDetectionInfo(fn)] = SpotDetectionSliceSNG(AlignedSlice, TEMPLATE, INFO.BrainSegmentationInfo(fn).BrainEdge, obj.ZFParameters);
                
                %update Computations
                obj.Computations(fn).Spots = reshape([SpotsDetected{fn}.Centroid], 2, numel(SpotsDetected{fn}))';
                obj.Computations(fn).Midbrain = fliplr(INFO.BrainSegmentationInfo(fn).BrainEdge);
                obj.Computations(fn).Counts = numel(SpotsDetected{fn});
                
                %update checkup after new spot computations
                if ~isempty(obj.checkup)
                    obj.checkup(fn).Spots = obj.SpotsInsiteArea(SpotParameters{fn}, obj.checkup(fn).Midbrain);
                    obj.checkup(fn).Counts = size(obj.checkup(fn).Spots, 1);
                end
            end
            
            obj.saveit
            save([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotDetectionInfo', '-append')
            save([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected', '-append')
            save([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters', '-append')
            
            obj.buildsheet
            
            delete(h)
        end
        function obj = FillComputations(obj, BrainSegmentationInfo, SpotsDetected)
            %temporary function if an obj is opened that is created with SpotNGlia 1.4.0 or before
            if isempty(obj.Computations)
                if ~exist('BrainSegmentationInfo', 'var')
                    load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo')
                end
                if ~exist('SpotsDetected', 'var')
                    load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected')
                end
                for k = 1:obj.nfishes
                    obj.Computations(k).Spots = reshape([SpotsDetected{k}.Centroid], 2, numel(SpotsDetected{k}))'; %#ok<IDISVAR,USENS>
                    obj.Computations(k).Midbrain = fliplr(BrainSegmentationInfo(k).BrainEdge);
                    obj.Computations(k).Counts = size(obj.Computations(k).Spots, 1);
                end
            end
        end
        function obj = FillCheckup(obj)
            %temporary function if an obj is opened that is created with SpotNGlia 1.4.0 or before
            if isempty(obj.checkup)
                load([obj.SavePath, '/', obj.InfoName, '.mat'], 'checkup')
                if exist('checkup', 'var')
                    obj.checkup = checkup;
                end
            end
        end
        function obj = CorrectCheckBrain(obj)
        %this function updates the spots according to previous added or removed spots
        
%         obj.checkup(k1).Spots
%         obj.checkup(k1).SpotAdditions
%         obj.checkup(k1).SpotRemovals
%         
%         
%         [obj.checkup(k1).Spots]
%         
%         
%         [TF, ~] = ismember([obj.checkup(k1).Spots], [obj.checkup(k1).SpotAdditions], 'rows');
%         sc.XData(TF) = [];
%         sc.YData(TF) = [];
%         %add previous added spots
%         sc.XData = [sc.XData, ph2.XData];
%         sc.YData = [sc.YData, ph2.YData];
        end
        function obj = CompleteProgram(obj, fishnumbers)
            
            [obj.StackInfo] = StackInfoSNG(obj.ImageInfo);
            
            if ~exist('fishnumbers', 'var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif (max(fishnumbers) > numel(obj.StackInfo))
                error('at least one fish does not exist in StackInfo')
            end
            
            obj = obj.PreProcession(fishnumbers);
            obj = obj.ExtendedDeptOfField(fishnumbers);
            obj = obj.Registration(fishnumbers);
            obj = obj.BrainSegmentation(fishnumbers);
            obj = obj.SpotDetection(fishnumbers);
            
            % obj.ShowFish([], true)
        end
        function obj = HistPar(obj, fishnumbers)
            %function to add a mean fish color parameter to the RegisterationInfo variable saved in INFO_...
            %its function can be removed later on as it is also added to Registration
            
            load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo')
            
            if ~exist('fishnumbers', 'var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif max(fishnumbers) > numel(obj.StackInfo)
                error('at least one fish does not exist in RegistrationInfo')
            end
            
            nfishes = numel(fishnumbers);
            for k1 = 1:nfishes
                fn = fishnumbers(k1);
                AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
                threshold = RegistrationInfo{k1}(strcmp({RegistrationInfo{k1}.name}, 'BackgroundThreshold')).value;
                
                %figure;imagesc(AlignedFish(:,700:end,:));
                
                for k3 = 1:3
                    Img = AlignedFish(:, 700:end, k3);
                    MeanFishColor(k1, k3) = mean(Img(:) < threshold(k3));
                    h = hist(Img(:), 0:1:255);
                    h2 = smoothdata(h, 'gaussian', 6);
                    [~, MaxFishColor(k1, k3)] = max(h2(1:floor(threshold(k3))));
                end
                
                row = find(strcmp({RegistrationInfo{k1}.name}, 'MeanFishColor'));
                RegistrationInfo{k1} = sng_StructFill(RegistrationInfo{k1}, {'Registration', 'Background Removal', 'MeanFishColor', MeanFishColor(k1, :)}, row);
                row = find(strcmp({RegistrationInfo{k1}.name}, 'MaxFishColor'));
                RegistrationInfo{k1} = sng_StructFill(RegistrationInfo{k1}, {'Registration', 'Background Removal', 'MaxFishColor', MaxFishColor(k1, :)}, row);
                
            end
            save([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo', '-append')
            
            obj.BatchInfo.MeanMaxFishColor = mean(MaxFishColor(:));
            obj.BatchInfo.StdMaxFishColor = std(MaxFishColor(:));
            obj.BatchInfo.MeanFishColor = mean(MeanFishColor(:));
            obj.BatchInfo.StdMeanFishColor = std(MeanFishColor(:));
        end
        
        ShowFishHeadHist(obj, fishnumber)
        ShowMaxFishHist(obj, fishnumber)
        obj = CheckFish(obj, ifish, INFO)
    end
    
    methods(Hidden = true)       
        Mask = BrainMask(obj, fishnumbers)
        %BrainOptimization(obj, fishnumbers)
        %obj2 = SpotOptimization2(obj, fishnumbers, zfinputlist)
        
        function obj = buildsheet(obj)
            
            obj = FillComputations(obj); %for pre SNG1.4.0
            obj = FillCheckup(obj);
            
            if ~isempty(obj.checkup)
                nspots2 = [obj.checkup.Counts];
                nspots2([obj.checkup.Include] ~= 1) = NaN;
                
                Sheet = [{obj.StackInfo.stackname}', ...
                    {obj.StackInfo.stacksize}', ...
                    num2cell([obj.Computations.Counts])', ...
                    num2cell(nspots2)'];
                title = {obj.InfoName, 'Images', 'Computed', 'Corrected'};
            else
                Sheet = [{obj.StackInfo.stackname}', ...
                    {obj.StackInfo.stacksize}', ...
                    num2cell([obj.Computations.Counts])'];
                title = {obj.InfoName, 'Images', 'Computed'};
            end
            ds = cell2dataset([title; Sheet]);
            export(ds, 'file', [obj.SavePath, '/', obj.InfoName, '.csv'], 'delimiter', obj.Delimiter)
        end
        function saveit(obj)
            %obj.savedate = datetime;
            %dt = strcat(string(year(date)),string(month(date)),string(day(date)));
            %firstim = obj.ImageInfo(1).name;
            %firstn = char(regexp(firstim,'\d+.tif','match'));
            %firstn = strrep(firstn,'.tif','');
            %name = strrep(firstim, [firstn,'.tif'],'');
            
            save(strcat(obj.SavePath, '/', obj.SaveName, '.mat'), 'obj');
            
            
            %{
                        load([obj.SavePath,'/',obj.InfoName,'.mat'], 'SpotsDetected')
                        if exist('SpotsDetected')
                            for k1 = 1:numel(SpotsDetected)
                                nspots(k1,1) = numel(SpotsDetected{k1});
                            end
 
                            Sheet = [{obj.StackInfo.stackname}',{obj.StackInfo.stacksize}',num2cell(nspots)]
                            title = {obj.InfoName,'images','nspots'}
 
                            ds = cell2dataset([title;Sheet]);
                            export(ds,'file',[P{1},'/',I{1},'.csv'],'delimiter',',')
                        end
            %}
        end
        function obj = LoadTemplate(obj)
            if isempty(obj.CompleteTemplate)
                obj.CompleteTemplate = load([obj.SourcePath, '/', 'Template3dpf', '.mat']);
            end
        end       
        function obj = LoadAnnotations(obj)
            
            %annotated midbrain roi
            if isempty(obj.AnnotatedMidBrainPath) || (obj.AnnotatedMidBrainPath(1) == 0)
                disp(['Select Annotated MidBrain Path for ',obj.SaveName])
                obj.AnnotatedMidBrainPath = uigetdir([], 'Select Annotated MidBrain Path');
            end
            
            %annotated spots
            if isempty(obj.AnnotatedSpotPath) || (obj.AnnotatedSpotPath(1) == 0)
                disp(['Select Annotated Spot Path for ',obj.SaveName])
                obj.AnnotatedSpotPath = uigetdir([], 'Select Annotated Spot Path');
            end
            
            %if one of the two paths is not empty and a char (dir) then load RegistrationInfo
            if (~isempty(obj.AnnotatedMidBrainPath) && ischar(obj.AnnotatedMidBrainPath)) || ...
                    (~isempty(obj.AnnotatedSpotPath) && ischar(obj.AnnotatedSpotPath))
                nfishes = numel(obj.StackInfo);
                tform_1234 = cell(nfishes, 1);
                %annotated brain and spots
                load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo')
                for k1 = 1:nfishes
                    tform_1234{k1} = RegistrationInfo{k1}(strcmp({RegistrationInfo{k1}.name}, 'tform_complete')).value;
                end
            end
            
            %if midbrain path exists in object than load, transform, save midbrain
            if (~isempty(obj.AnnotatedMidBrainPath) && ischar(obj.AnnotatedMidBrainPath))
                ambr = cell(nfishes, 1);
                for k1 = 1:nfishes
                    if exist([obj.AnnotatedMidBrainPath, '/', obj.StackInfo(k1).stackname, '.zip'], 'file')
                        RoiBrain = ReadImageJROI([obj.AnnotatedMidBrainPath, '/', obj.StackInfo(k1).stackname, '.zip']);
                        amb = sng_roicell2poly(RoiBrain, 1);
                        [ambr{k1}(:, 1), ambr{k1}(:, 2)] = transformPointsForward(tform_1234{k1}, amb(:, 1), amb(:, 2));
                        ambr{k1} = double(ambr{k1});
                    end
                end
                [obj.Annotations(1:nfishes).MidBrain] = ambr{:};
            end
            
            %if spot path exists in object than load, transform, save spot
            if (~isempty(obj.AnnotatedSpotPath) && ischar(obj.AnnotatedSpotPath))
                SpotAnn = cell(nfishes, 1);
                for k1 = 1:nfishes
                    
                    if exist([obj.AnnotatedSpotPath, '/', obj.StackInfo(k1).stackname, '.roi'], 'file')
                        RoiMicroglia = ReadImageJROI([obj.AnnotatedSpotPath, '/', obj.StackInfo(k1).stackname, '.roi']);
                        SpotAnn{k1} = zeros(size(RoiMicroglia.mfCoordinates, 1), 2);
                        [SpotAnn{k1}(:, 1), SpotAnn{k1}(:, 2)] = transformPointsForward(tform_1234{k1}, ...
                            RoiMicroglia.mfCoordinates(:, 1), RoiMicroglia.mfCoordinates(:, 2));
                    end
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
        function obj = BrainVal(obj, BrainSegmentationInfo, fishnumbers)
            %fishnumbers has to correspondent with BrainSegmentationInfo if given
            if ~exist('BrainSegmentationInfo', 'var')
                load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo')
            end
            
            if ~exist('fishnumbers', 'var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif (max(fishnumbers) > numel(obj.StackInfo))
                error('at least one fish does not exist in StackInfo')
            end
            nfishes = numel(fishnumbers);
            
            
            Jaccard = zeros(nfishes, 1);
            Dice = zeros(nfishes, 1);
            
            cmbr = {BrainSegmentationInfo.BrainEdge};
            ambr = {obj.Annotations.MidBrain};
            
            for k1 = 1:nfishes
                fn = fishnumbers(k1);
                
                mx = round(max([cmbr{k1}; fliplr(ambr{fn})]));
                
                a = poly2mask(cmbr{k1}(:, 2), cmbr{k1}(:, 1), mx(1), mx(2));
                b = poly2mask(ambr{fn}(:, 1), ambr{fn}(:, 2), mx(1), mx(2));
                
                intersection = (a & b);
                union = (a | b);
                
                Jaccard(k1) = sum(intersection(:)) / sum(union(:));
                %Overlap(k1) = sum(intersection(:)) / min(sum(a(:)),sum(b(:)));
                Dice(k1) = (2 * sum(intersection(:))) / (sum(union(:)) + sum(intersection(:)));
                
            end
            
            %fill the BrainInfo Structure (per fish)
            temp = num2cell(Jaccard); [obj.BrainInfo(1:nfishes).Jaccard] = temp{:};
            temp = num2cell(Dice); [obj.BrainInfo(1:nfishes).Dice] = temp{:};
            
            %fill BrainStats Structure, statistical values
            obj.BrainStats.FiveNumberSummaryJaccard = sng_FiveNumberSum(Jaccard);
            obj.BrainStats.MeanJaccard = mean(Jaccard);
            
            obj.BrainStats.FiveNumberSummaryDice = sng_FiveNumberSum(Dice);
            obj.BrainStats.MeanDice = mean(Dice);
            
        end
        function obj = SpotVal(obj, SpotParameters)
            
            if ~exist('SpotParameters', 'var')
                load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters');
            end
            %checks first if Annoations exist, than of Midbrain exist, than of all fiels are empty
            if isempty(obj.Annotations) || ~isfield(obj.Annotations, 'MidBrain') || all(cellfun(@isempty, {obj.Annotations.MidBrain}))
                AMBR = {obj.checkup.Midbrain};
            else
                AMBR = {obj.Annotations.MidBrain};
            end
            
            
            if sum(~cellfun('isempty', {obj.Annotations.Spots})) ~= numel(SpotParameters)
                warning('not all fishes are annotated by hand')
            end
            fishnumbers = 1:numel(SpotParameters);
            %removes fishnumbers from the list if not all annotations are known
            fishnumbers = fishnumbers(~cellfun('isempty', {obj.Annotations.Spots}));
            nfishes = numel(fishnumbers);
            
            LinkDistance = cell(nfishes, 1);
            CorrectSpots = cell(nfishes, 1);
            nCorrect = zeros(nfishes, 1);
            FalsePosSpots = cell(nfishes, 1);
            nFalsePos = zeros(nfishes, 1);
            FalseNegSpots = cell(nfishes, 1);
            nFalseNeg = zeros(nfishes, 1);
            Precision = zeros(nfishes, 1);
            Recall = zeros(nfishes, 1);
            F1score = zeros(nfishes, 1);
            AbsDifference = zeros(nfishes, 1);
            RelDifference = zeros(nfishes, 1);
            
            SpotCom = cell(nfishes, 1);
            ambsS = cell(nfishes, 1);
            
            %value added to prevent for substraction by zero
            a = 0.001;
            
            
            for k1 = 1:nfishes
                Spotpar = SpotParameters{k1};
                %ambr = obj.Annotations(k1).MidBrain;
                ambs = obj.Annotations(k1).Spots;
                ambr = AMBR{k1};
                
                if ~isempty(Spotpar)
                    
                    %select spot insite annotated brainregion
                    [spotcentroids] = reshape([Spotpar.Centroid], 2, numel(Spotpar))';
                    [ins, ~] = inpolygon(spotcentroids(:, 1), spotcentroids(:, 2), ambr(:, 1), ambr(:, 2));
                    
                    SpotparS = Spotpar( ...
                        ins' & ...
                        [Spotpar.MinProbability] & ... %selection of spotinfo based on conditions
                        [Spotpar.LargerThan] & ...
                        [Spotpar.SmallerThan]);
                    
                    SpotCom{k1} = reshape([SpotparS.Centroid], 2, numel(SpotparS))';
                else
                    SpotCom{k1} = [];
                end
                
                %because the annotated midbrain does not corresponds
                %exactly with the annotated spots, de annotated spots are
                %also filtered to compute the performance of the
                %spotdetection only
                [ins2, ~] = inpolygon(ambs(:, 1), ambs(:, 2), ambr(:, 1), ambr(:, 2));
                ambsx = ambs(:, 1); ambsx = ambsx(ins2);
                ambsy = ambs(:, 2); ambsy = ambsy(ins2);
                ambsS{k1} = [ambsx, ambsy];
                
                [Correct, FalsePos, FalseNeg, link] = sng_CoordinateMatching ...
                    (SpotCom{k1}, ambsS{k1}, 10);
                
                LinkDistance{k1} = link;
                
                if exist('Correct', 'var')
                    CorrectSpots{k1} = Correct;
                    nCorrect(k1) = size(Correct, 1);
                else
                    CorrectSpots{k1} = [];
                    nCorrect(k1) = 0;
                end
                if exist('FalsePos', 'var')
                    FalsePosSpots{k1} = FalsePos;
                    nFalsePos(k1) = size(FalsePos, 1);
                else
                    FalsePosSpots{k1} = [];
                    nFalsePos(k1) = 0;
                end
                if exist('FalseNeg', 'var')
                    FalseNegSpots{k1} = FalseNeg;
                    nFalseNeg(k1) = size(FalseNeg, 1);
                else
                    FalseNegSpots{k1} = [];
                    nFalseNeg(k1) = 0;
                end
                
                Precision(k1) = nCorrect(k1) / (nCorrect(k1) + nFalsePos(k1) + a);
                Recall(k1) = nCorrect(k1) / (nCorrect(k1) + nFalseNeg(k1) + a);
                F1score(k1) = 2 * (Precision(k1) * Recall(k1)) / (Precision(k1) + Recall(k1) + a);
                AbsDifference(k1) = nFalsePos(k1) - nFalseNeg(k1);
                RelDifference(k1) = (nFalsePos(k1) - nFalseNeg(k1)) / nCorrect(k1);
                
            end
            obj.SpotStats.FiveNumberSummaryPrecision = sng_FiveNumberSum(Precision);
            obj.SpotStats.MeanPrecision = mean(Precision);
            
            obj.SpotStats.FiveNumberSummaryRecall = sng_FiveNumberSum(Recall);
            obj.SpotStats.MeanRecall = mean(Recall);
            
            obj.SpotStats.FiveNumberSummaryF1score = sng_FiveNumberSum(F1score);
            obj.SpotStats.MeanF1score = mean(F1score);
            
            obj.SpotStats.FiveNumberSummaryAbsDifference = sng_FiveNumberSum(AbsDifference);
            obj.SpotStats.MeanAbsDifference = mean(AbsDifference);
            
            obj.SpotStats.FiveNumberSummaryRelDifference = sng_FiveNumberSum(RelDifference);
            obj.SpotStats.MeanRelDifference = mean(RelDifference);
            
            
            %store variables in obj.SpotInfo
            [obj.SpotInfo(1:nfishes).LinkDistance] = LinkDistance{:};
            [obj.SpotInfo(1:nfishes).CorrectSpots] = CorrectSpots{:};
            temp = num2cell(nCorrect); [obj.SpotInfo(1:nfishes).nCorrect] = temp{:};
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
        function obj = TtestVal(obj, SpotN)
            %computes t-test if the number of spots are known and saved ouder
            %if checkup is generated, TtestVal will compute on those and store in
            %obj.SpotBrainStats. This is not very pretty and should be changed later
            
            obj = FillCheckup(obj); %pre1.4.0
            obj = FillComputations(obj);
            
            if ~exist('SpotN', 'var')
                SpotN = [obj.Computations.Counts];
                
                %would be easier if SpotN containes a counts column
                for k1 = 1:numel(obj.checkup)
                    SpotN(k1) = size(obj.checkup(k1).Spots, 1);
                    if obj.checkup(k1).Include == 0
                        SpotN(k1) = NaN;
                    end
                end
            end
            
            if ~isfield((obj.Annotations), 'Counts')
                for l2 = 1:numel(obj.Annotations)
                    obj.Annotations(l2).Counts = size(obj.Annotations(l2).Spots, 1);
                end
            end
            
            SpotAnnN = [obj.Annotations.Counts];
            
            
            if numel(SpotN) ~= numel(SpotAnnN)
                warning('number of computed fishes is different from annotated fishes');
                SpotN = SpotN(1:numel(SpotAnnN));
            end
            
            ind = logical(~isnan(SpotAnnN).* ~isnan(SpotN));
            
            %removes emptys and nans
            SAN = SpotAnnN(ind);
            SN = SpotN(ind);
            
            [h, p] = ttest(SAN, SN);
            mn = mean(abs(SN-SAN));
            st = std(abs(SN-SAN));
            RMSD = sqrt(mean((SAN - SN).^2));
            
            obj.SpotBrainStats.ttest = h;
            obj.SpotBrainStats.ttestval = p;
            obj.SpotBrainStats.MeanAbsDifference = mn;
            obj.SpotBrainStats.StdAbsDifference = st;
            obj.SpotBrainStats.RMSD = RMSD;
            
        end
        function obj = SpotBrainVal(obj)
            
            %load([obj.SavePath,'/',obj.InfoName,'.mat'],'SpotParameters')
            load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected')
            
            nfishes = numel(obj.StackInfo);
            
            LinkDistance = cell(nfishes, 1);
            CorrectSpots = cell(nfishes, 1);
            nCorrect = zeros(nfishes, 1);
            FalsePosSpots = cell(nfishes, 1);
            nFalsePos = zeros(nfishes, 1);
            FalseNegSpots = cell(nfishes, 1);
            nFalseNeg = zeros(nfishes, 1);
            Precision = zeros(nfishes, 1);
            Recall = zeros(nfishes, 1);
            F1score = zeros(nfishes, 1);
            AbsDifference = zeros(nfishes, 1);
            RelDifference = zeros(nfishes, 1);
            
            %SpotCom = cell(nfishes,1);
            
            %value added to prevent for substraction by zero
            a = 0.001;
            
            for k1 = 1:nfishes
                
                Spotpar = SpotsDetected{k1};
                %ambr = obj.BrainInfo(k1).AnnotatedMidBrainRegistrated;
                
                SpotCom = reshape([Spotpar.Centroid], 2, numel(Spotpar))';
                ambs = obj.Annotations(k1).Spots;
                
                
                [Correct, FalsePos, FalseNeg, link] = sng_CoordinateMatching ...
                    (SpotCom, ambs, 10);
                
                LinkDistance{k1} = link;
                
                if exist('Correct', 'var')
                    CorrectSpots{k1} = Correct;
                    nCorrect(k1) = size(Correct, 1);
                else
                    CorrectSpots{k1} = [];
                    nCorrect(k1) = 0;
                end
                if exist('FalsePos', 'var')
                    FalsePosSpots{k1} = FalsePos;
                    nFalsePos(k1) = size(FalsePos, 1);
                else
                    FalsePosSpots{k1} = [];
                    nFalsePos(k1) = 0;
                end
                if exist('FalseNeg', 'var')
                    FalseNegSpots{k1} = FalseNeg;
                    nFalseNeg(k1) = size(FalseNeg, 1);
                else
                    FalseNegSpots{k1} = [];
                    nFalseNeg(k1) = 0;
                end
                
                Precision(k1) = nCorrect(k1) / (nCorrect(k1) + nFalsePos(k1) + a);
                Recall(k1) = nCorrect(k1) / (nCorrect(k1) + nFalseNeg(k1) + a);
                F1score(k1) = 2 * (Precision(k1) * Recall(k1)) / (Precision(k1) + Recall(k1) + a);
                AbsDifference(k1) = abs(nFalsePos(k1)-nFalseNeg(k1));
                RelDifference(k1) = abs(nFalsePos(k1)-nFalseNeg(k1)) / nCorrect(k1);
            end
            obj.SpotBrainStats.FiveNumberSummaryPrecision = sng_FiveNumberSum(Precision);
            obj.SpotBrainStats.MeanPrecision = mean(Precision);
            obj.SpotBrainStats.FiveNumberSummaryRecall = sng_FiveNumberSum(Recall);
            obj.SpotBrainStats.MeanRecall = mean(Recall);
            obj.SpotBrainStats.FiveNumberSummaryF1score = sng_FiveNumberSum(F1score);
            obj.SpotBrainStats.MeanF1score = mean(F1score);
            obj.SpotBrainStats.FiveNumberSummaryAbsDifference = sng_FiveNumberSum(AbsDifference);
            obj.SpotBrainStats.MeanAbsDifference = mean(abs(AbsDifference));
            obj.SpotBrainStats.FiveNumberSummaryRelDifference = sng_FiveNumberSum(RelDifference);
            obj.SpotBrainStats.MeanRelDifference = mean(RelDifference);
            
            %store variables in obj.SpotInfo
            [obj.SpotBrainInfo(1:nfishes).LinkDistance] = LinkDistance{:};
            [obj.SpotBrainInfo(1:nfishes).CorrectSpots] = CorrectSpots{:};
            temp = num2cell(nCorrect); [obj.SpotBrainInfo(1:nfishes).nCorrect] = temp{:};
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
        
        function ShowBoxPlot(obj, exportit)
            %Example show and save
            %   ShowBoxPlot(obj,1)
            %Example only show
            %   obj.ShowBoxPlot
            
            %%% brainval spotval
            fsx = 6; fsy = 10;
            
            if isfield(obj.BrainInfo, 'Jaccard')
                bj = obj.BrainStats.FiveNumberSummaryJaccard(3);
                bd = obj.BrainStats.FiveNumberSummaryDice(3);
                
                [h1, g1] = setfigax1;
                boxplot(g1, [ ...
                    obj.BrainInfo.Jaccard; ...
                    obj.BrainInfo.Dice]', {'Jaccard', 'Dice'});
                %title('Brain Validation')
                
                setfigax2(g1)
                text(1.2, bj, sprintf('%.2f', bj), 'FontName', 'arial', 'FontSize', 8)
                text(2.2, bd, sprintf('%.2f', bd), 'FontName', 'arial', 'FontSize', 8)
                realsizeandsave(h1)
                
                
            end
            
            if isfield(obj.SpotInfo, 'Precision')
                sp = obj.SpotStats.FiveNumberSummaryPrecision(3);
                sr = obj.SpotStats.FiveNumberSummaryRecall(3);
                sf = obj.SpotStats.FiveNumberSummaryF1score(3);
                [h2, g2] = setfigax1;
                boxplot(g2, [ ...
                    obj.SpotInfo.Precision; ...
                    obj.SpotInfo.Recall; ...
                    obj.SpotInfo.F1score]', {'Precision', 'Recall', 'F1score'});
                %title('Spot Validation')
                setfigax2(g2)
                text(1.3, sp, sprintf('%.2f', sp), 'FontName', 'arial', 'FontSize', 8)
                text(2.3, sr, sprintf('%.2f', sr), 'FontName', 'arial', 'FontSize', 8)
                text(3.3, sf, sprintf('%.2f', sf), 'FontName', 'arial', 'FontSize', 8)
                set(gca, 'XLim', [0.5, 3.8])
                realsizeandsave(h2)
            end
            
            if isfield(obj.SpotBrainInfo, 'Precision')
                bsp = obj.SpotBrainStats.FiveNumberSummaryPrecision(3);
                bsr = obj.SpotBrainStats.FiveNumberSummaryRecall(3);
                bsf = obj.SpotBrainStats.FiveNumberSummaryF1score(3);
                
                [h3, g3] = setfigax1;
                boxplot(g3, [ ...
                    obj.SpotBrainInfo.Precision; ...
                    obj.SpotBrainInfo.Recall; ...
                    obj.SpotBrainInfo.F1score]', {'Precision', 'Recall', 'F1score'});
                %title('Spot Validation on Computed Brain')
                
                setfigax2(g3)
                
                text(1.3, bsp, sprintf('%.2f', bsp), 'FontName', 'arial', 'FontSize', 8)
                text(2.3, bsr, sprintf('%.2f', bsr), 'FontName', 'arial', 'FontSize', 8)
                text(3.3, bsf, sprintf('%.2f', bsf), 'FontName', 'arial', 'FontSize', 8)
                set(gca, 'XLim', [0.5, 3.8])
                realsizeandsave(h3)
            end
            
            function [figurehandle, axishandle] = setfigax1
                %fsx = 5;fsy = 10;
                figurehandle = figure('PaperUnits', 'centimeters', 'Color', [1, 1, 1]);
                axishandle = gca;
                sng_figcm(fsx, fsy);
            end
            function setfigax2(axishandle)
                set(axishandle, 'FontName', 'arial', 'FontSize', 8, 'XGrid', 'on', 'YGrid', 'on');
                set(axishandle, 'Units', 'centimeters', 'Position', [1.2, 1.2, fsx - 1.7, fsy - 1.7]);
                set(axishandle, 'YLim', [-0.02, 1.02]);
            end
            function realsizeandsave(figurehandle)
                if exist('exportit', 'var') && exportit
                    export_fig(figurehandle, ['/Users/samuelgeurts/Desktop/', 'boxplot', num2str(figurehandle.Number)], '-png', '-r600', '-nocrop');
                end
                ScaledFigure.calibrateDisplay(113.6); %113.6 for this screen
                ScaledFigure(figurehandle, 'reuse');
                set(figurehandle, 'Units', 'Centimeters');
                set(figurehandle, 'Position', (get(figurehandle, 'Position') + [fsx, 0, 0, 0]));
            end
            
            
        end
        function ShowFishVal(obj, fishnumbers, ComOrAnn)
            %if ComOrAnn is not given, the computed brain is used
            
            if ~exist('fishnumbers', 'var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif (max(fishnumbers) > numel(obj.StackInfo))
                error('at least one fish does not exist in StackInfo')
            end
            
            %CompOrAnn determines if the computed or the annotated brain is
            %displayed inlcudin ghte Correct/FalsePos/FalseNeg
            
            load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo')
            
            %if ~isempty({obj.Annotations.MidBrain})
            %    ambr = {obj.Annotations.MidBrain};
            %end
            
            if isempty(obj.Annotations) || ~isfield(obj.Annotations, 'MidBrain')
                AMBR = {obj.checkup.Midbrain};
            else
                AMBR = {obj.Annotations.MidBrain};
            end
            
            for k1 = 1:numel(fishnumbers)
                fn = fishnumbers(k1);
                
                cmbr = BrainSegmentationInfo(fn).BrainEdge;
                AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
                
                %%
                figure; imshow(uint8(AlignedFish))
                hold on;
                plot(cmbr(:, 2), cmbr(:, 1), 'Color', [255, 75, 75]/255, 'LineWidth', 2);
                
                if exist('ambr', 'var')
                    plot(AMBR{fn}(:, 1), AMBR{fn}(:, 2), 'Color', [75, 75, 255]/255, 'LineWidth', 2);
                end
                
                if ~exist('ComOrAnn', 'var') || strcmp(ComOrAnn, 'Com')
                    Correct = obj.SpotBrainInfo(fn).CorrectSpots;
                    FalsePos = obj.SpotBrainInfo(fn).FalsePosSpots;
                    FalseNeg = obj.SpotBrainInfo(fn).FalseNegSpots;
                    
                    ps = obj.SpotBrainInfo(fn).Precision; %used for textbox
                    rs = obj.SpotBrainInfo(fn).Recall;
                    fs = obj.SpotBrainInfo(fn).F1score;
                    
                elseif strcmp(ComOrAnn, 'Ann')
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
                
                scatter([Correct(:, 1); FalsePos(:, 1)], [Correct(:, 2); FalsePos(:, 2)], 400, 'LineWidth', 2, 'MarkerEdgeColor', 1/255*[255, 75, 75]);
                scatter([Correct(:, 1); FalseNeg(:, 1)], [Correct(:, 2); FalseNeg(:, 2)], 270, 'LineWidth', 2, 'MarkerEdgeColor', 1/255*[75, 75, 255]);
                legend('Computed brain', 'Annotated brain', 'Computed spots', 'Annotated spots')
                
                
                nc = size(Correct, 1) + size(FalsePos, 1); %number of computed spots
                na = size(Correct, 1) + size(FalseNeg, 1); %number of annotated spots
                jb = obj.BrainInfo(fn(k1)).Jaccard; %Jaccard number
                db = obj.BrainInfo(fn(k1)).Dice; % Dice number
                
                
                str = sprintf(['BRAIN PARAMETERS', ...
                    '\n', 'Jaccard: %.2f', ...
                    '\n', 'dice: %.2f', ...
                    '\n\n', 'SPOT PARAMETERS', ...
                    '\n', 'Precision: %.2f', ...
                    '\n', 'Recall: %.2f', ...
                    '\n', 'F1score: %.2f', ...
                    '\n\n', 'Computed spots: %.0f', ...
                    '\n', 'Annotated spots: %.0f', ...
                    '\n\n', 'Correct spots: %.0f', ...
                    '\n', 'False Positives: %.0f', ...
                    '\n', 'False Negatives: %.0f'], ...
                    jb, db, ps, rs, fs, nc, na, ...
                    size(Correct, 1), ...
                    size(FalsePos, 1), ...
                    size(FalseNeg, 1));
                
                annotation('textbox', [.1, .63, .7, .3], 'String', str, 'FitBoxToText', 'on');
                
            end
        end
        function ShowFish(obj, fishnumbers, saveimage)
            
            if ~exist('fishnumbers', 'var') || isempty(fishnumbers)
                fishnumbers = 1:numel(obj.StackInfo);
            elseif (max(fishnumbers) > numel(obj.StackInfo))
                error('at least one fish does not exist in StackInfo')
            end
            
            if exist('saveimage', 'var') && saveimage
                TempFolderName = ([obj.SavePath, '/', 'FinalResult']);
                if ~exist(TempFolderName, 'dir')
                    mkdir(TempFolderName)
                end
            end
            
            load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo')
            load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected')
            
            for k1 = 1:numel(fishnumbers)
                fn = fishnumbers(k1);
                
                cmbr = BrainSegmentationInfo(fn).BrainEdge;
                AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
                cmbs = reshape([SpotsDetected{fn}.Centroid], 2, numel(SpotsDetected{fn}))';
                
                if exist('saveimage', 'var') && saveimage
                    figure('visible', 'off')
                else
                    figure
                end
                imshow(uint8(AlignedFish), 'Border', 'tight', 'InitialMagnification', 50)
                %imshow(uint8(AlignedFish))
                hold on;
                plot(cmbr(:, 2), cmbr(:, 1), 'Color', [255, 75, 75]/255, 'LineWidth', 2);
                
                scatter(cmbs(:, 1), cmbs(:, 2), 400, 'LineWidth', 2, 'MarkerEdgeColor', 1/255*[255, 75, 75]);
                
                nc = size(cmbs, 1); %number of computed spots
                str = sprintf('Computed Spots: %.0f', nc);
                annotation('textbox', [.1, .63, .7, .3], 'String', str, 'FitBoxToText', 'on', 'FontSize', 15);
                
                if exist('saveimage', 'var') && saveimage
                    saveas(gcf, [TempFolderName, '/', obj.StackInfo(fn).stackname], 'tif');
                    close(gcf)
                end
            end
        end
        
        function obj = Dummy(obj,input)
           disp('Dummy function') 
           if exist('input','var')
               disp(num2str(input))
           end
           disp([obj.SaveName])
        end
    end
    
    methods(Static)
        function SpotCoordsFiltered = SpotsInsiteArea(SpotParametersSingle, SpotCoords)
            %selects only Spots instite area 'SpotParametersSingle'
            
            %{
            SpotParametersSingle = SpotParameters{ifish};
            SpotCoords = [X3',Y3']
            %}
            
            Xrow = SpotCoords(:, 1);
            Yrow = SpotCoords(:, 2);
            
            %% compute new spots insite area
            rc = reshape([SpotParametersSingle.Centroid], 2, numel(SpotParametersSingle))';
            [in, ~] = inpolygon(rc(:, 1), rc(:, 2), Xrow, Yrow);
            
            SpotsDetec = SpotParametersSingle(in' == 1 & ...
                [SpotParametersSingle.LargerThan] == 1 & ...
                [SpotParametersSingle.SmallerThan] == 1 & ...
                [SpotParametersSingle.MinProbability] == 1);
            
            [SpotCoordsFiltered] = reshape([SpotsDetec.Centroid], 2, numel(SpotsDetec))';
        end
    end
    
    
end

