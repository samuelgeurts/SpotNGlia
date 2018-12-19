function [Img2,crop] = sng_RGB_IATwarp2(Img1,ECCWarp)
%warp image according to warp info and remove not overlapping boundary
%Version sng_RGB_IATwarp2
%   change cropping method and add output crop values

%{
%Example
Img1 = ImageSlice{1}(:,:,1:3);
a = ECCWarp{1}{1}
b = ECCWarp{1}{2}
clear ECCWarp
ECCWarp{1} = a
ECCWarp{2} = b
%}

%preallocation
[M,N,Oh] = size(Img1);

[Img1(:,:,2), supportECC{1}] = iat_inverse_warping(Img1(:,:,2), ECCWarp{1},'translation' , 1:N, 1:M);
[Img1(:,:,3), supportECC{2}] = iat_inverse_warping(Img1(:,:,3), ECCWarp{2},'translation' , 1:N, 1:M);
support = supportECC{1} .* supportECC{2};

%{
figure;imagesc(Img1)
figure;imagesc(support)
%}

%old crop method
%{
%crop image, find the largest rectangle consist of ones
A = numel(find(support(round(M/2),1:end)));
B = numel(find(support(1:end,round(N/2))));
Img2 = reshape(Img1(repmat(logical(support),1,1,3)),B,A,3);        
%figure;imshow(uint8(Img2{j}));
%}

xr1 = find(support(round(M/2),1:end),1,'first');
xr2 = find(support(round(M/2),1:end),1,'last');
yr1 = find(support(1:end,round(N/2)),1,'first');
yr2 = find(support(1:end,round(N/2)),1,'last');

Img2 = Img1(yr1:yr2,xr1:xr2,:);
crop = [yr1,yr2;xr1,xr2];


end

