function [ECCWarp,Img2] = sng_RGB2(Img1,scale,transform,levels,iterations)
%RGB2 does te same als RGB but the output functions are turned for speed
%increase as you only want to now ECCWarp
    

%set parameters
%{

scale = 1/4;
par.transform = 'translation'; 
par.levels = 2;
par.iterations = 10; %iterations per level
Img1 = tempImage;
%}

par.transform = transform;
par.levels = levels;
par.iterations = iterations; %iterations per level

Img1sc = imresize(Img1,scale);

ECCWarp{1} = iat_ecc(Img1sc(:,:,2), Img1sc(:,:,1), par)*(1/scale);
ECCWarp{2} = iat_ecc(Img1sc(:,:,3), Img1sc(:,:,1), par)*(1/scale);
% ECCWarp = iat_LucasKanade(Img1sc(:,:,j), Img1sc(:,:,1), par)*(1/scale);



if nargout == 2
    [M,N,Oh] = size(Img1);
    [Img1(:,:,2), supportECC{1}] = iat_inverse_warping(Img1(:,:,2), ECCWarp{1},'translation' , 1:N, 1:M);
    [Img1(:,:,3), supportECC{2}] = iat_inverse_warping(Img1(:,:,3), ECCWarp{2},'translation' , 1:N, 1:M);
    support = supportECC{1} .* supportECC{2};
    %crop image, find the largest rectangle consist of ones
    A = numel(find(support(round(M/2),1:end)));
    B = numel(find(support(1:end,round(N/2))));
    Img2 = reshape(Img1(repmat(logical(support),1,1,3)),B,A,3);        
    %figure;imagesc(Img2);
end

end