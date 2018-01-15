function ReleaseSpotNGlia

SavePath = uigetdir

[fList, pList] = matlab.codetools.requiredFilesAndProducts('SpotNGlia.m')


%% save dependencies to specific folder


List = fList'
expression = ['\', filesep];

for k1 = 1:numel(List)
    
    splitStr = regexp(fList{k1}, expression, 'split');
    
    if strcmp(splitStr{6}, '3p_functions')
        %seperate folder for third party functions
        
        path = [SavePath, filesep, splitStr{6}, filesep, splitStr{7}];
        if ~exist(path, 'dir')
            mkdir(path)
        end
        copyfile(fList{k1}, path)
    elseif strcmp(splitStr{6}, 'SpotNGlia')
        %separate folder for SpotNGlia functions
        path = SavePath;
        for k2 = 6:numel(splitStr) - 1
            path = [path, filesep, splitStr{k2}]; %#ok<AGROW>
        end
        if ~exist(path, 'dir')
            mkdir(path)
        end
        copyfile(fList{k1}, path)
    elseif strcmp(splitStr{6}, 'sng_functions')
        %separate folder for supporting sng functions
        path = [SavePath, filesep, 'sng_functions'];
        if ~exist(path, 'dir')
            mkdir(path)
        end
    else
        %functions not on paths above are saved directly to the SavePath
        path = SavePath;
    end
    
    copyfile(fList{k1}, path)
end

%save source file
BasePath = '/Users/samuelgeurts/Dropbox';
SourcePath = [BasePath, '/', 'SpotNGlia Source'];
SourcePathNew = [SavePath,filesep,'Source']


mkdir(SourcePathNew)
copyfile([SourcePath,filesep,'Template3dpf.mat'], SourcePathNew)
copyfile([SourcePath,filesep,'zfinput.mat'], SourcePathNew)

%save script files
copyfile('/Users/samuelgeurts/Dropbox/My Matlab/SpotNGlia/User Methods SpotNGlia.m',SavePath);


end