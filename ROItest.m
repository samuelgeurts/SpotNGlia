[file path] = uigetfile('\\bioinf-filesrv2\cluster15\Ham\Research_group\Data voor Samuel\Counts')
%[file path] = uigetfile('/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Data Nynke/')
ROI = ReadImageJROI([path file])
s = ROI.vnSlices
cc = ROI.mfCoordinates

%[file path] = uigetfile('/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Data Nynke/')
[file path] = uigetfile('\\bioinf-filesrv2\cluster15\Ham\Research_group\Data voor Samuel\outbox')
info = imfinfo([path file])
for k = 1:numel(info)
    
    sls = find(s == k)
    
    Img{k} = imread([path file],k)
    figure;imagesc(Img{k})
    hold on
    scatter(cc(sls,1),cc(sls,2))
end

s = ROI.vnSlices
cc = ROI.mfCoordinates

