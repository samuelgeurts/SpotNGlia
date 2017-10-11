function [Ispots,spoutput] = SpotDetectionLink2(Ialligned,CompleteTemplate,cmbr,zfinput)

%{

todo zfinput to zebrafish__
cmbr = fliplr(BrainAnn)
zfinput = obj.zfinput
Ialligned = Ialigned{k5};
%}


if ~exist('zfinput','var');
    ColorToGrayVector = CompleteTemplate.SpotContrastVector;
   
    zfinput = struct;
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','RgbToGray','ColorToGrayVector',ColorToGrayVector,'');  %select color channel [0 1 0],     
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleBase',0.5,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleLevels',9,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','Kthreshold',0,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPlevels',4:9,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPthreshold',50,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MinSpotSize',30,''); % size selection 
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MaxSpotSize',500,''); % size selection 
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MinProbability',0.01,''); %color selection
end

sng_zfinputAssign(zfinput,'SpotDetection')

%TODO: suffisticated color transform of afterwards color detection
%TODO: apply on original images which are not dept of field combined as :some spots appear due to dof artifacts. Combine images afterwards
%TODO:  filter spots based on color
%TODO:  add more parameters
%TODO:  find out how hard thresholding works


%% RGB to gray algorithm

GI = sng_RGB2Gray(Ialligned,ColorToGrayVector,false);


%% Generate Wavelets

%[MultiProductTh,MultiProduct,A,W,Wth] = sng_SpotWavelet(GI,ScaleLevels,ScaleBase,MPlevels,MPthreshold,false);
[MultiProductTh,MultiProduct] = sng_SpotWavelet(GI,ScaleLevels,ScaleBase,MPlevels,MPthreshold,false);

%{
figure;imagesc(Ialligned)

figure;imagesc(MultiProductTh)

%}
%% SpotMeasure

%function sng_spotproperties (Ialligned, MaskedIm)

CC = bwconncomp(MultiProductTh);
L = labelmatrix(CC);
nspots1 = CC.NumObjects;
Regions1 = regionprops(L,'Area','PixelIdxList','PixelList','Centroid','Image','BoundingBox');






%% Color mean of founded spots
colormean = zeros(numel(Regions1),3);
for k = 1:numel(Regions1)
    spotimage = getfield(Regions1,{k},'Image');
    sbb = round(getfield(Regions1,{k},'BoundingBox')); %round because wc?
    spotcimage = Ialligned(sbb(2):sbb(2)+sbb(4)-1,sbb(1):sbb(1)+sbb(3)-1,1:3);
    spotcimage(repmat(spotimage,1,1,3) == 0) = 0; %color image of spot
    %color channels
    rc = spotcimage(:,:,1);
    gc = spotcimage(:,:,2);
    bc = spotcimage(:,:,3);
    %spot color mean
    colormean(k,:) = [mean(rc(spotimage == 1)),mean(gc(spotimage == 1)),mean(bc(spotimage == 1))];
    %{
    figure;imagesc(spotimage);axis equal tight off;colormap gray
    figure;imagesc(spotcimage);axis equal tight off
    %}
end
temp = num2cell(colormean,2);[Regions1.ColorMean] = temp{:};

%color probability of mean spot color
vec = round(reshape([Regions1.ColorMean],3,numel(Regions1))');
colorindex = sub2ind([255,255,255],vec(:,1),vec(:,2),vec(:,3));


vec = round(reshape([Regions1.ColorMean],3,numel(Regions1))')+1;
colorindex = sub2ind([256,256,256],vec(:,1),vec(:,2),vec(:,3));






%ismember(uint32(colorindex),CompleteTemplate.SVAP_index)

probability = zeros(numel(Regions1),1);
for k2 = 1:numel(Regions1)
    in = find(CompleteTemplate.SVAP_index == colorindex(k2));
    if ~isempty(in)
        probability(k2,1) = CompleteTemplate.SpotVectorArrayProbability(in);
    end
end
temp = num2cell(probability);[Regions1.ColorProbability] = temp{:};

%{
[x1,y1,z1] = ind2sub([255,255,255],colorindex)
[x2,y2,z2] = ind2sub([255,255,255],CompleteTemplate.SVAP_index);

%figure;scatter3(x1,y1,z1)
kk = boundary(x2,y2,z2); %takes long ~1min
figure;t = trisurf(kk,x2,y2,z2,'Facecolor','red','FaceAlpha',0.1,'EdgeColor','none')
hold on
scatter3(colormean(:,1),colormean(:,2),colormean(:,3))
%}

%% other discriminations

cm = [Regions1.ColorProbability] >= MinProbability;
temp = num2cell(cm);[Regions1.MinProbability] = temp{:};


%% determine if insite brain polygon
[rc] = reshape([Regions1.Centroid],2,nspots1)';
[in,~] = inpolygon(rc(:,1),rc(:,2),cmbr(:,2),cmbr(:,1));
temp = num2cell(in);[Regions1.Insite] = temp{:};

%{
figure;
imagesc(Ialligned)
plot(cmbr(:,2),cmbr(:,1))
hold on
plot(rc(in,1),rc(in,2),'r+')
plot(rc(on,1),rc(on,2),'g+')
plot(rc(~in,1),rc(~in,2),'bo')

figure;hist([Regions1.Area],30)
%}

%% spotsize

%spots larger than MinSpotSize
lt = [Regions1.Area] >= MinSpotSize;
temp = num2cell(lt);[Regions1.LargerThan] = temp{:};

%spots smaller than MaxSpotSize
st = [Regions1.Area] <= MaxSpotSize;
temp = num2cell(st);[Regions1.SmallerThan] = temp{:};
 
%%
Ifilt = ismember(L,find(...
    [Regions1.Insite] == 1 &...
    [Regions1.LargerThan] == 1 &...
    [Regions1.SmallerThan] == 1 &...
    [Regions1.MinProbability] == 1));

nspots = sum([Regions1.Insite] == 1 &...
    [Regions1.LargerThan] == 1 &...
    [Regions1.SmallerThan] == 1 &...
    [Regions1.MinProbability] == 1);





%{
figure;imagesc(Ifilt);axis off tight equal
    hold on;plot(cmbr(:,2),cmbr(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
colormap gray
set(gca,'position',[0 0 1 1],'units','normalized')
truesize
figure;scatter3(colormean(:,1),colormean(:,2),colormean(:,3))
%}        

Ispots = Ialligned;
Ispots(repmat(Ifilt,1,1,3) == 0) = 0;

%{
    figure;imagesc(uint8(Ispots));axis off tight equal
    hold on;plot(cmbr(:,2),cmbr(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
    
    
    colormap gray
    set(gca,'position',[0 0 1 1],'units','normalized')
    truesize

   Ispotshsv =  rgb2hsv(Ispots);
   figure;imagesc(Ispotshsv);axis off tight equal

   
   Ihsv1 = Ispotshsv(:,:,1);
   Ihsv2 = Ispotshsv(:,:,2);
   Ihsv3 = Ispotshsv(:,:,3);
    figure;hist(Ihsv1(Ihsv1~=0))
    figure;hist(Ihsv2(Ihsv2~=0))
    figure;hist(Ihsv3(Ihsv3~=0))

    
    
    Ehsi=sng_rgb2hsi2(colormean')'
    pca1 = pca(colormean','NumComponents',1)
     
   
    
    SpotCom = cell2mat({Regions2.Centroid}')
    
    figure;imagesc(uint8(Ispots));axis off tight equal
    hold on;
    scatter(SpotCom(pca1<=0.07,1),SpotCom(pca1<=0.07,2),500,'blue')
    scatter(SpotCom(pca1>=0.25,1),SpotCom(pca1>=0.25,2),500,'yellow')


%}        




spoutput = struct('stage',[],'substage',[],'name',[]','value',[]);

spoutput = sng_StructFill(spoutput,{'SpotDetection','Wavelet','Multiproduct',MultiProduct});
spoutput = sng_StructFill(spoutput,{'SpotDetection','Wavelet','MultiproductThreshold',MultiProductTh});


spoutput = sng_StructFill(spoutput,{'SpotDetection','SpotSelection','SpotParameters',Regions1});
spoutput = sng_StructFill(spoutput,{'SpotDetection','SpotSelection','SpotFilter',Ifilt});
spoutput = sng_StructFill(spoutput,{'SpotDetection','SpotSelection','NumberOfSpots',nspots});








