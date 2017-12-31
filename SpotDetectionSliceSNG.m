function [SpotsDetected, SpotParameters, SpotDetectionInfo, SpotN] = SpotDetectionSliceSNG(Ialigned, CompleteTemplate, cmbr, ZFParameters)


%Version SpotDetectionLink3
% correctiong for color and index 255 vs 256

%{
   Ialigned = AlignedFish;
   CompleteTemplate = obj.CompleteTemplate;
   cmbr = BrainSegmentationInfo(k1).BrainEdge;
   ZFParameters = obj.ZFParameters;
   ZFParameters = ZFParametersTemp{1}
%}
%{
   Ialigned = AlignedFish{k1};
   CompleteTemplate = obj.CompleteTemplate;
   cmbr = MidBrain;
   ZFParameters =  ZFParametersTemp{1}
%}
%{
   Ialigned = AlignedSlice{1};
   cmbr = BrainSegmentationInfo(k1).BrainEdge;
   ZFParameters = obj.ZFParameters;
   CompleteTemplate = TEMPLATE
%}
limit = 5;
sng_zfinputAssign(ZFParameters, 'SpotDetection');
ScaleLevels = max(MPlevels);

%CubeSizeSpotErode = 4;
%SE1 = strel('cube',CubeSizeSpotErode); %spot structure element for erode


%TODO: suffisticated color transform or afterwards color detection
%TODO: add more parameters
%TODO: find out how hard thresholding works
%TODO: apply machinelearning on spot images or add backgroundinfo for new feature

%put image is cell if it is not
if ~iscell(Ialigned)
    sz = size(Ialigned)
    Ialigned = mat2cell(Ialigned, sz(1), sz(2), sz(3));
end


%% Apply Multiproduct on slices an achieve spot paramaters
for k3 = 1:numel(Ialigned)
    %RGB to gray algorithm
    GI = sng_RGB2Gray(Ialigned{k3}, ColorToGrayVector, false);
    %Generate Wavelets
    %[MultiProductTh,MultiProduct,A,W,Wth] = sng_SpotWavelet(GI,ScaleLevels,ScaleBase,MPlevels,MPthreshold,false);
    [MultiProductTh{k3}, MultiProduct] = sng_SpotWaveletQuick(GI, ScaleLevels, ScaleBase, MPlevels, MPthreshold, false);
    %{
       figure;imagesc(Ialigned{1})
       figure;imagesc(MultiProductTh)
 
       for k = 1:size(W,3)
       figure;imagesc(W(:,:,k))
       end
    %}
    %SpotMeasures
    CC = bwconncomp(MultiProductTh{k3});
    L = labelmatrix(CC);
    nspots1 = CC.NumObjects;
    Regions0{k3} = regionprops(L, 'Area', 'PixelIdxList', 'PixelList', 'Centroid', 'Image', 'BoundingBox');
    
    %add max value of the Multiproduct peak and slice number
    for k1 = 1:numel(Regions0{k3})
        Regions0{k3}(k1).MPPeak = max(MultiProduct(Regions0{k3}(k1).PixelIdxList));
        Regions0{k3}(k1).Slice = k3;
    end
end

%% sort out double counted spots

if numel(Ialigned) > 1   
    Regions1 = Regions0{k3};
    for k3 = 1:numel(Ialigned)-1
        Regions1 = RemoveDoubleSpots(Regions1, Regions0{k3+1});
    end
else
    Regions1 = Regions0;
end

    function RegNew = RemoveDoubleSpots(Reg1, Reg2)
        %{
        Reg1 = Regions1
        Reg2 = Regions0{k3+1}
        %}        
        rc1 = reshape([Reg1.Centroid], 2, numel(Reg1))';
        rc2 = reshape([Reg2.Centroid], 2, numel(Reg2))';
        src1 = size(rc1, 1);
        src2 = size(rc2, 1);
        
        distance = zeros(src1, src2);
        for k70 = 1:src1
            distance(k70, :) = sqrt((rc2(:, 1) - rc1(k70, 1)).^2+(rc2(:, 2) - rc1(k70, 2)).^2);
        end
        %indexrc2 gives for every rc1 spot the closest rc2 point
        [mindist, indexrc2] = min(distance, [], 2);
        
        %sort the minimal distance in accending order
        [sm, I] = sort(mindist);
        
        RegNew = Reg1;
        
        I(sm <= limit) %coupled indices of rc1
        indexrc2(I(sm <= limit)) %coupled indices of rc2
        
        %change spots in RegNew from Reg1 to Reg2 that has a higher peakvalue in Reg2
        for k71 = 1:src1
            if sm(k71) <= limit
                p1 = Reg1(I(k71)).MPPeak;
                p2 = Reg2(indexrc2(I(k71))).MPPeak;
                if p2 >= p1
                    RegNew(I(k71)) = Reg2(indexrc2(I(k71)));
                else
                    %RegNew(I(k71)) = Reg1(I(k71));
                end
                %Reg1(I(k71)).Centroid
                %Reg2(indexrc2(I(k71))).Centroid
            end
        end
        
        %add spots that only occur in Reg2
        rc2unique = 1:src2;
        rc2unique(indexrc2(I(sm <= limit))) = [];              
        for k3 = numel(rc2unique):-1:1   
            RegNew(end+1) = Reg2(rc2unique(k3));
        end
    end


%% Color mean of founded spots and neighbour background

%BackgroundEroded = imerode(~MultiProductTh,SE1);
%SpotEroded = imerode(MultiProductTh,SE1);
%figure;imagesc(BackgroundEroded)
%figure;imagesc(SpotEroded)

%As the boundary of the spot is often not well defined we want to give
%more weight to the inner pixel for computing the mean color. A distance
%map is a good solution. Erosion could leave non or to less pixels for
%small spots, while a distance map does not have that problem.

for k3 = 1:numel(Ialigned)
DistanceSpotMap{k3} = bwdist(~MultiProductTh{k3});
DistanceBackgroundMap{k3} = bwdist(MultiProductTh{k3});
end

%{
      figure;imagesc(MultiProductTh)
      figure;imagesc(DistanceSpotMap)
      figure;imagesc(DistanceBackgroundMap)
%}

SpotCmean = zeros(numel(Regions1), 3);
BackGroundCmean = zeros(numel(Regions1), 3);

for k1 = 1:numel(Regions1)
    %spotimage = getfield(Regions1, {k1}, 'Image');
    sbb = round(getfield(Regions1, {k1}, 'BoundingBox'));
    slc = round(getfield(Regions1, {k1}, 'Slice'));
    
    % get a slightly larger box to be able to compute the meanbackgroun more precise
    sbb = sbb + [-2, -2, 4, 4];
    
    spotcimage = Ialigned{slc}(sbb(2):sbb(2) + sbb(4) - 1, sbb(1):sbb(1) + sbb(3) - 1, 1:3);
    DSpot = DistanceSpotMap{slc}(sbb(2):sbb(2)+sbb(4)-1, sbb(1):sbb(1)+sbb(3)-1);
    DBackground = DistanceBackgroundMap{slc}(sbb(2):sbb(2)+sbb(4)-1, sbb(1):sbb(1)+sbb(3)-1);
    %computes a weighted mean of the spot color based on the distance map
    SpotCmean(k1, 1:3) = sum(sum((repmat(DSpot, 1, 1, 3) / sum(DSpot(:))).*single(spotcimage)));
    BackGroundCmean(k1, 1:3) = sum(sum((repmat(DBackground, 1, 1, 3) / sum(DBackground(:))).*single(spotcimage)));
    
    %{
          %use this for machine learning
          c = round(getfield(Regions1, {k1}, 'Centroid'));
          spotcimage = Ialigned{k3}(c(2)-9:c(2)+9, c(1)-9:c(1)+9, 1:3);
    %}
end
contrastvector = BackGroundCmean - SpotCmean;
[azimuth, elevation, r] = cart2sph(contrastvector(:, 1), contrastvector(:, 2), contrastvector(:, 3));

temp = num2cell(SpotCmean, 2); [Regions1.ColorMean] = temp{:};
temp = num2cell(BackGroundCmean, 2); [Regions1.BackRoundMean] = temp{:};
temp = num2cell(azimuth, 2); [Regions1.azimuth] = temp{:};
temp = num2cell(elevation, 2); [Regions1.elevation] = temp{:};
temp = num2cell(r); [Regions1.r] = temp{:};

%{
      figure;imagesc(spotimage);axis equal tight off;colormap gray
      figure;imagesc(spotcimage);axis equal tight off
%}


%color probability of mean spot color !!! old not well working feature
%because it is not complete, there are more batches with different lightning
%update template generation and use machine learning or background information
vec = round(reshape([Regions1.ColorMean], 3, numel(Regions1))');
%correction index vs color
colorindex = sub2ind([256, 256, 256], vec(:, 1)+1, vec(:, 2)+1, vec(:, 3)+1);
probability = zeros(numel(Regions1), 1);
for k2 = 1:numel(Regions1)
    in = find(CompleteTemplate.SVAP_index == colorindex(k2));
    if ~isempty(in)
        probability(k2, 1) = CompleteTemplate.SpotVectorArrayProbability(in);
    end
end
temp = num2cell(probability); [Regions1.ColorProbability] = temp{:};

%{
       %shows
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
temp = num2cell(cm); [Regions1.MinProbability] = temp{:};


%% determine if insite brain polygon
[rc] = reshape([Regions1.Centroid], 2, numel(Regions1))';
[in, ~] = inpolygon(rc(:, 1), rc(:, 2), cmbr(:, 2), cmbr(:, 1));
temp = num2cell(in); [Regions1.Insite] = temp{:};

%{
   figure;
   imagesc(Ialigned{1})
   hold on
   plot(cmbr(:,2),cmbr(:,1))
   plot(rc(in,1),rc(in,2),'r+')
   %plot(rc(on,1),rc(on,2),'g+')
   plot(rc(~in,1),rc(~in,2),'bo')
 
   figure;hist([Regions1.Area],30)
%}

%% spotsize

%spots larger than MinSpotSize
lt = [Regions1.Area] >= MinSpotSize;
temp = num2cell(lt); [Regions1.LargerThan] = temp{:};

%spots smaller than MaxSpotSize
st = [Regions1.Area] <= MaxSpotSize;
temp = num2cell(st); [Regions1.SmallerThan] = temp{:};

%%
%{
   Ifilt = ismember(L,find(...
       [Regions1.Insite] == 1 &...
       [Regions1.LargerThan] == 1 &...
       [Regions1.SmallerThan] == 1 &...
       [Regions1.MinProbability] == 1));
%}

%nspots = sum([Regions1.Insite] == 1 &...
%    [Regions1.LargerThan] == 1 &...
%    [Regions1.SmallerThan] == 1 &...
%    [Regions1.MinProbability] == 1);


SpotsDetected = Regions1([Regions1.Insite] == 1 & ...
    [Regions1.LargerThan] == 1 & ...
    [Regions1.SmallerThan] == 1 & ...
    [Regions1.MinProbability] == 1);


%{
   figure;imagesc(Ifilt);axis off tight equal
       hold on;plot(cmbr(:,2),cmbr(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
   colormap gray
   set(gca,'position',[0 0 1 1],'units','normalized')
   truesize
   figure;scatter3(colormean(:,1),colormean(:,2),colormean(:,3))
    %}
       
    %Ispots = Ialigned;
    %Ispots(repmat(Ifilt,1,1,3) == 0) = 0;  
    
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
    

    SpotParameters = Regions1;    
    SpotDetectionInfo.Multiproduct = MultiProduct;
    SpotDetectionInfo.MultiproductThreshold = MultiProductTh;
    %SpotDetectionInfo.SpotFilter = Ifilt;   
    SpotN = numel(SpotsDetected);
    
end