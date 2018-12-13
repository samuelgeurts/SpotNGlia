function [ECCWarp,cellImg1] = sng_AllignFishBatch3(cellImg1,scale,transform,levels,iterations,varargin)
%allign all fishes with different dept of field from location 
%the first fish acts as a template
%imput is at least 2 images in a cell structure
%iatool.net/ for more info
%by the translation method, information has been lost, consider to use only
%the translation values 
%examples
%sng_AllignFishBatch3(cellImg1,scale,transform,levels,iterations,Fullpath)   


% scale = 1/4;
% transform = 'affine'; %'translation'; 
% levels = 3;
% iterations = 20; %iterations per level
% cellImg1 = Meanfish;

%algorithm parameters
par.transform = transform;
par.levels = levels;
par.iterations = iterations; %iterations per level

Img1sc = imresize(cellImg1{1},scale); %resize to speed up
Img1sc = Img1sc(1:round(size(Img1sc,1)*(11/12)),:,1:3); %to exclude scalebar

[M,N,Oh] = size(cellImg1{1});
support = true(M,N);
ECCWarp = cell(1,numel(cellImg1));

for j = 2:numel(cellImg1)
    Img2sc = imresize(cellImg1{j},scale);
    %to remove the scale bar
    Img2sc = Img2sc(1:round(size(Img2sc,1)*(11/12)),:,1:3);
    ECCWarp{j} = iat_ecc(Img2sc, Img1sc, par)*(1/scale);
end

if nargin >= 6 || nargout >= 2
    for j = 2:numel(cellImg1)
        [cellImg1{j}, supportECC] = iat_inverse_warping(cellImg1{j}, ECCWarp{j}, par.transform, 1:N, 1:M);
        cellImg1{j} = uint8(cellImg1{j});
        support = support & logical(supportECC);
    end
    %crop image, find the largest rectangle consist of ones
    A = numel(find(support(round(M/2),1:end)));
    B = numel(find(support(1:end,round(N/2))));
    
    for j = 1:numel(cellImg1)
        cellImg1{j} = reshape(cellImg1{j}(repmat(support,1,1,Oh)),B,A,Oh);
        if nargin >= 6
            Fullpath = varargin{1};
            imwrite(uint8(cellImg1{j}),Fullpath,...
                'WriteMode', 'append', 'Compression','none');
        end
        %figure;imshow(uint8(cellImg1{j}));
    end
end
    
    
    
    
    
    
end


% for k=1:numel(cellImg1)
%     figure;imagesc(uint8(Meanfish{k}))        
%     figure;imagesc(uint8(cellImg1{k}))
%     figure;imagesc(support);
% end
% sng_figureslide
