if ~exist('FolderPath10')
    disp('select loadpath')
    FolderPath10 = uigetdir('/Volumes/Macintosh SSD/Users/samuelgeurts/Dropbox','select folder for loading')
end

load([FolderPath10,'/imageinfo.mat'],'imageinfo')
load([FolderPath10,'/stackinfo.mat'],'stackinfo')
load([FolderPath10,'/input.mat'],'zfinput')


for k1 = 1:numel(stackinfo)
    
    
        qqq = imread([FolderPath10,'/',stackinfo(k1).stackname,'-BrainSegmentation.tif']); 
        bwbb = bwboundaries(qqq);


        stackinfo(k1).BrainSegmentation.BrainEdge = bwbb{1}
        
        
        figure;imagesc(qqq)
end
        
    save([FolderPath,'/stackinfo.mat'],'stackinfo')
        
        
        a/(a + b + c)