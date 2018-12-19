function Icombined = sng_SliceCombine(ImgCell,IndexMatrix)
%creates combined image based on Indexmatrix called in sng_StackDof

Icombined = zeros(size(ImgCell{1}),'uint8');
for k = 1:numel(ImgCell)
    ImgCell{k}(~repmat(IndexMatrix(:,:,k),1,1,3)) = 0;
    Icombined = Icombined + ImgCell{k};
end