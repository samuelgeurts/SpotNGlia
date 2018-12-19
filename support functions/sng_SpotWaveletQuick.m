function [MultiProductTh,MultiProduct,A,W,Wth] = sng_SpotWavelet(GI,ScaleLevels,ScaleBase,MPlevels,MPthreshold,image_TF)
% this function computes the multiproduct of the wavelet based
% spotdetection. Scale dependent thresholding is not applied in the
% version.
% sng_SpotWaveletQuick comes from sng_SpotWavelet but it computes only the 
% wavelet planes used for the multiproduct
% example:
%{
[MultiProductTh,MultiProduct,A,W,Wth] = sng_SpotWavelet(GI,8,0.5,3:8,3000)

%}
% % % % % ScaleLevels = max(MPlevels)
%Scalelevels became useless 
nplanes = numel(MPlevels); %number of waveletplanes
[sx,sy] = size(GI);
GI = double(GI);
%preallocation
A = zeros(sx,sy,nplanes+1);

gsigma = [2.^(ScaleBase*[MPlevels(1)-1,MPlevels])];
for k = 1:nplanes+1
    %gaussian filter from A0   
    H = fspecial('Gaussian',2*ceil(2*gsigma(k))+1,gsigma(k));
    A(:,:,k) = imfilter(GI,H,'symmetric');
end
    % wavelet coefficients
    W = A(:,:,1:end-1)-A(:,:,2:end);
    

% Hard Thresholding (scale dependent) <-doesnt work well,fix!

    %sigma{i} = 1.4826*mad((W{i}(:))',1);       %median absolute deviation
    %t{i} = Kthreshold*sigma{i};
    %Wth{i} = ((W{i})>=t{i}).*W{i};               %bright spots of dark background
    %Wth{i} = ((W{i})<=t{i}).*W{i};                %dark spots of bright background
    %Wth{i} = (abs(W{i})<=t{i}).*W{i};                %both %does not work well
    %sng_lim(Wth{i})
                
   Wth = abs(((W)<= 0).*W);   %dark spots of bright background

%% Multiproduct

MultiProduct = prod(Wth(:,:,1:nplanes),3);
MultiProductTh = MultiProduct >= MPthreshold;

%{
    figure;imagesc(MultiProduct);axis off tight equal
    figure;imagesc(abs(log(MultiProduct+1)));axis off tight equal   

    figure;imagesc(MultiProductTh);axis off tight equal
    colormap gray
    set(gca,'position',[0 0 1 1],'units','normalized')
    truesize

    bwb = bwboundaries(Ibrain);

    figure;imagesc(MultiProductTh);axis off tight equal
    hold on;plot(bwb{1}(:,2),bwb{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',4)
    colormap gray
    set(gca,'position',[0 0 1 1],'units','normalized')
    truesize    

    figure;imagesc(abs(log(Multiproduct+1)));    
%}

%{
    figure;imagesc(MaskedIm);axis off tight equal
    hold on;plot(bwb{1}(:,2),bwb{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',4)
    colormap gray
    set(gca,'position',[0 0 1 1],'units','normalized')
    truesize
%}    
end