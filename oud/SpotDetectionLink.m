function [Ispots,Parameters] = SpotDetectionLink(Ibrain,Ialligned,zfinput)


if ~exist('zfinput','var');
    zfinput = struct;
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','RgbToGray','GrayMethod',2);  %select color channel
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleBase',0.5,'');  %select color channel
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleLevels',8,'');  %select color channel
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','Kthreshold',3,'');  %select color channel
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','CCthreshold',1,'');  %select color channel


end
sng_zfinputAssign(zfinput,'SpotDetection')

%TODO: suffisticated color transform of afterwards color detection
%% RGB to gray algorithm
if isnumeric(GrayMethod)

    GI = Ialligned(:,:,GrayMethod);
    %temporary choose green channel
    A0 = double(GI);
end



for i = 1:ScaleLevels
    gsigma{i} = 2^(ScaleBase*(i));
    %gaussian filter from A0
    H{i} = fspecial('Gaussian',2*ceil(2*gsigma{i})+1,gsigma{i});
    A{i} = imfilter(A0,H{i},'symmetric');
    % wavelet coefficients
    if i==1; W{i} = A0-A{i}; else W{i} = A{i-1}-A{i}; end;

%% Spot detection

%TODO: apply on original images which are not dept of field combined as 
%some spots appear due to dof artifacts. Combine images afterwards

    Wnt{i}=W{i};
    % Hard Thresholding (scale dependent)
    sigma{i} = 1.4826*mad((W{i}(:))',1);       %median absolute deviation
    t{i} = Kthreshold*sigma{i};
    %W{i} = ((W{i})>=t{i}).*W{i};               %bright spots of dark background
    W{i} = ((W{i})<=t{i}).*W{i};                %dark spots of bright background
    %W{i} = (abs(W{i})<=t{i}).*W{i};                %both

    % correlation image P{i} is a multiplication of Wavelate planes 1 to i
    if i==1; P{i} = W{i}; else P{i} = W{i}.*P{i-1};end;
    % correlation image P3{i} is a multiplication of Wavelate planes i-2 to i   
    if i>=3
        P3{i}=(W{i-2}.*W{i-1}.*W{i});
    end;
end

for i = 1:ScaleLevels
    % discrimination from background
    P{i} = 255*(abs(P{i})>=CCthreshold);
    if i>=3;P3{i} = 255*(abs(P3{i})>=CCthreshold);end;
end

%{
    figure;
    subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.01], [0.01 0.01], [0.01 0.01]);
    subplot(3,3,1);l=imagesc(A0);colormap(gray);
    set(gca,'xtick',[],'ytick',[]);
    axis equal tight
    text(sct(2)*0.90,sct(1)*0.05,'A0','FontSize',20,'Color',[0 0 0])
    for i = 1:lvls
        subplot(3,3,i+1);imagesc(A{i}); colormap(gray);
        set(gca,'xtick',[],'ytick',[]);
        text(sct(2)*0.90,sct(1)*0.05,['A' , num2str(i)],'FontSize',20,'Color',[0 0 0])
        axis equal tight
    end
    truesize


    figure;
    for i = 1:lvls
        subplot(3,3,i+1);imagesc(Wnt{i});colormap(gray);
        set(gca,'xtick',[],'ytick',[]);
        text(sct(2)*0.90,sct(1)*0.05,['W' , num2str(i)],'FontSize',20,'Color',[0 0 0])
        axis equal tight
    end
    truesize

    
    figure;
    for i = 1:lvls
        subplot(3,3,i+1);imagesc(-W{i});colormap(gray);
        text(sct(2)*0.90,sct(1)*0.05,['W' , num2str(i)],'FontSize',20,'Color',[1 1 1])
        axis equal tight
        set(gca,'xtick',[],'ytick',[]);
    end
    truesize
%}

P4_8=W{3}.*W{4}.*W{5}.*W{6}.*W{7}.*W{8};

%{
    figure;imagesc(P4_8);axis off tight equal
%}

P4_8th = 255*(abs(P4_8)>=3000);

%{
    figure;imagesc(P4_8th);axis off tight equal
    colormap gray
    set(gca,'position',[0 0 1 1],'units','normalized')
    truesize

    bwb = bwboundaries(Ibrain);

    figure;imagesc(P4_8th);axis off tight equal
    hold on;plot(bwb{1}(:,2),bwb{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',4)
    colormap gray
    set(gca,'position',[0 0 1 1],'units','normalized')
    truesize    

%}

MaskedIm = P4_8th & Ibrain;

%{ 
    figure;imagesc(MaskedIm);axis off tight equal
    hold on;plot(bwb{1}(:,2),bwb{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',4)
    colormap gray
    set(gca,'position',[0 0 1 1],'units','normalized')
    truesize
%}    

%measure spots
CC = bwconncomp(MaskedIm);
L = labelmatrix(CC);
nspots1 = CC.NumObjects;
Regions1 = regionprops(MaskedIm,'Area');
%{
figure;hist([Regions.Area],30)
%}

%remove spots with size lower than minspotsize
minspotsize = 30;
Ifilt = ismember(L,find([Regions1.Area] >= minspotsize));
Regions2 = regionprops(Ifilt,'Area','PixelIdxList','PixelList','Centroid','Image','BoundingBox');
nspots2 = numel(Regions2);

%compute color mean of found spots
for k = 1:numel(Regions2)
    spotimage = getfield(Regions2,{k},'Image');
    sbb = round(getfield(Regions2,{k},'BoundingBox')); %round because wc?
    spotcimage = Ialligned(sbb(2):sbb(2)+sbb(4)-1,sbb(1):sbb(1)+sbb(3)-1,1:3);
    spotcimage(repmat(spotimage,1,1,3) == 0) = 0; %color image of spot
    %color channels
    rc = spotcimage(:,:,1);gc = spotcimage(:,:,2);bc = spotcimage(:,:,3);
    %spot color mean
    colormean(k,:) = [mean(rc(spotimage == 1)),mean(gc(spotimage == 1)),mean(bc(spotimage == 1))];

 
    %figure;imagesc(spotimage);axis equal tight off;colormap gray
    %figure;imagesc(spotcimage);axis equal tight off
    
end


%{
figure;hist(colormean(:,1))
figure;hist(colormean(:,2))
figure;hist(colormean(:,3))

figure;imagesc(Ifilt);axis off tight equal
hold on;plot(bwb{1}(:,2),bwb{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',4)
colormap gray
set(gca,'position',[0 0 1 1],'units','normalized')
truesize



%}        

Ispots = Ialligned;
Ispots(repmat(Ifilt,1,1,3) == 0) = 0;

%{
    figure;imagesc(uint8(Ispots));axis off tight equal
    hold on;plot(bwb{1}(:,2),bwb{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
    colormap gray
    set(gca,'position',[0 0 1 1],'units','normalized')
    truesize
%}        

%TODO:  filter spots based on color
%TODO:  add more parameters


Parameters.Scalebase = ScaleBase;
Parameters.Levels = ScaleLevels;
Parameters.KThreshold = k;
Parameters.CorrelationThreshold = CCthreshold;

Parameters.SpotProperties = Regions2;
Parameters.Colormean = colormean;
Parameters.NumberofSpots = nspots2;






