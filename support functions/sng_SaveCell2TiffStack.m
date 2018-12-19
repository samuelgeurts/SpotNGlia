function sng_SaveCell2TiffStack(CorrectedSlice,FullSavePath)
%save cell structure with image to single tiff file
%overwrites if already exist


%FullSavePath = [SavePath,'/',stackinfo(k1).stackname,'.tif']

SavePath = fileparts(FullSavePath);

if ~exist(SavePath,'dir');
    mkdir(SavePath);
end

imwrite(uint8(CorrectedSlice{1}),FullSavePath, 'WriteMode', 'overwrite',  'Compression','none'); 
for k3 = 2:numel(CorrectedSlice)
    imwrite(uint8(CorrectedSlice{k3}),FullSavePath,'WriteMode','append','Compression','none');
end


end


