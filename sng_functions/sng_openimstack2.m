function Imstack = sng_openimstack2(Filename)
% opens stacked tif images and saves in cellstructure

% Filename = 'STACKcorrect0026-0027-0028.tif'
% imread('Filename');           %only opens first image of stack

info = imfinfo(Filename);
t=Tiff(Filename,'r');

% %image in tiff file is saved in strips
% %if image contains 3 stacks, there are 3 offsets
% info.StripOffsets;           %show stripoffsets
% offset1 = info.StripOffsets; %first offset
% sbc = info.StripByteCounts;  %

%does not work because stacks are not saved as subfiles
%offsets = getTag('t.TagID,SubIFD')     

Imstack = cell(1,numel(info));

for i = 1:numel(info)
    setDirectory(t,i)
    %Imstack{i} = t.read();
    Imstack{i} = t.readRGBAImage();
    Imstack{i} = im2uint8(Imstack{i});
end

end