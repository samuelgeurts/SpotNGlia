function sng_SaveCell2TiffStack(CorrectedSlice,FullSavePath)
%save cell structure with image to single tiff file
%overwrites if already exist

%FullSavePath = [SavePath,'/',imageinfo(1).batchname,'.tif']

imwrite(uint8(CorrectedSlice{1}),FullSavePath, 'WriteMode', 'overwrite',  'Compression','none'); 
for k3 = 2:numel(CorrectedSlice)
    imwrite(uint8(CorrectedSlice{k3}),[FullSavePath,'.tif'],'WriteMode','append','Compression','none');
end


end