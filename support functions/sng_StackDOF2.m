function [IndexMatrix, variance_sq, Icombined] = sng_StackDOF2(ImgCell,win_size,filt)
%Version sng_StackDOF2
%   change the slicen ouput to a logic variable IndexMatrix which costs
%   less memory and can be applied faster to acchieve the Icombined image

%{
Img = CorrectedSlice;
win_size = variancedisksize
%}

if ~exist('win_size','var')
    win_size = 3;
end
if ~exist('filt','var')
    filt = 'max';
end

n = numel(ImgCell);

[M,N,~] = size(ImgCell{1});
variance_sq = zeros(M,N,size(ImgCell,1));


% Making spherical structuring element for determing the local variance
% with radisu win_size
 struc_el = fspecial('disk', win_size)>0;
 
 % Calculating the local variance for each pixel and each slice
 for i = 1:n 
    variance = stdfilt(ImgCell{i}, struc_el);
    variance_sq(:,:,i) = sum(variance,3);
 
    %variance_sq(:,:,i) = sqrt(sum(variance.^2,3));

 end
 
% Making an MxN image with the slice indexes of the maximum/minimum local variance values
if strcmp(filt,'min')
    [~,slicen] = min(variance_sq,[],3);    
elseif strcmp(filt,'max')   
    [~,slicen] = max(variance_sq,[],3);
else
    error('wronginput')
end

%create a logic image which selects pixels form the different slices
IndexMatrix = false(size(variance_sq));
Btemp = (slicen(:) - 1) * M*N; 
Ctemp = 1:M*N;
IndexMatrix(Btemp+Ctemp') = true;

if nargout >= 3
Icombined = sng_SliceCombine(ImgCell,IndexMatrix);
end



%{
figure;imagesc(Img{1})
figure;imagesc(Img{2})
figure;imagesc(Img{3})
figure;imagesc(Img{4})


figure;imagesc(variance_sq(:,:,1))
figure;imagesc(variance_sq(:,:,2))
figure;imagesc(variance_sq(:,:,3))
figure;imagesc(variance_sq(:,:,4))

figure;imagesc(Icombined)

%}
%stacked double used
%value win_size is char, should be int

