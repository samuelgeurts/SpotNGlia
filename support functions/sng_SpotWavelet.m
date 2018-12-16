function [MultiProductTh,MultiProduct,A,W,Wth] = sng_SpotWavelet(GI,ScaleLevels,ScaleBase,MPlevels,MPthreshold,image_TF)
% this function computes the multiproduct of the wavelet based
% spotdetection. Scale dependent thresholding is not applied in the
% version.
% example:
%{
[MultiProductTh,MultiProduct,A,W,Wth] = sng_SpotWavelet(GI,8,0.5,3:8,3000)

%}
% % % % % ScaleLevels = max(MPlevels)

[sx,sy] = size(GI);

%preallocation
A = zeros(sx,sy,ScaleLevels);


A(:,:,1) = double(GI);
gsigma = [1,2.^(ScaleBase*(1:ScaleLevels))];
for k = 2:ScaleLevels+1
    %gaussian filter from A0
    H = fspecial('Gaussian',2*ceil(2*gsigma(k))+1,gsigma(k));
    A(:,:,k) = imfilter(A(:,:,1),H,'symmetric');
end
    % wavelet coefficients
    W = A(:,:,1:ScaleLevels)-A(:,:,2:ScaleLevels+1);

% Hard Thresholding (scale dependent) <-doesnt work well,fix!

    %sigma{i} = 1.4826*mad((W{i}(:))',1);       %median absolute deviation
    %t{i} = Kthreshold*sigma{i};
    %Wth{i} = ((W{i})>=t{i}).*W{i};               %bright spots of dark background
    %Wth{i} = ((W{i})<=t{i}).*W{i};                %dark spots of bright background
    %Wth{i} = (abs(W{i})<=t{i}).*W{i};                %both %does not work well
    %sng_lim(Wth{i})
                
   Wth = abs(((W)<= 0).*W);   %dark spots of bright background

if image_TF
 
    ss = ceil(sqrt(ScaleLevels+1)) %subplot size
    
    figure;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.01 0.01]);
    for i = 1:ScaleLevels+1
        subplot(ss,ss,i);imagesc(A(:,:,i)); colormap(gray);
        set(gca,'xtick',[],'ytick',[]);
        text(sy*0.90,sx*0.05,['A' , num2str(i-1)],'FontSize',20,'Color',[0 0 0])
        axis equal tight
    end
    truesize

    figure;
    for i = 1:ScaleLevels
        subplot(ss,ss,i+1);imagesc(W(:,:,i));colormap(gray);
        set(gca,'xtick',[],'ytick',[]);
        text(sy*0.90,sx*0.05,['W' , num2str(i)],'FontSize',20,'Color',[0 0 0])
        axis equal tight
    end
    truesize
    
    figure;
    for i = 1:ScaleLevels
        subplot(ss,ss,i+1);imagesc(Wth(:,:,i));colormap(gray);
        text(sy*0.90,sx*0.05,['W' , num2str(i)],'FontSize',20,'Color',[1 1 1])
        axis equal tight
        set(gca,'xtick',[],'ytick',[]);
    end
    truesize
end

%% Multiproduct

%OriginalImage = A(:,:,1) + sum(Wth(:,:,1),3);

MultiProduct = prod(Wth(:,:,MPlevels),3);
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