classdef SNGSpotDetection < handle
    %extern functions needed
    %sng_smoothzoom -for movie
    %sng_zoom
    %sng_NCC
    %IAT toolbox
    
    properties %(Access = private)
        %Object
        SpotNGliaObject
        %INPUTPARAMETERS
        ColorToGrayVector
        ScaleBase
        ScaleLevels
        Kthreshold
        MPlevels
        MPthreshold
        MinSpotSize
        MaxSpotSize
        MinProbability
        SpotDistLimit
        ComputeOnSlice
        %from CompleteTemplate
        template = []
        
        %from SpotNGlia
        iFish
        imageNames
        imagePath
    end
    properties
        %OUTPUTPARAMETERS
        
        
    end
    properties(Transient = true, Hidden = true)
        %IMAGES
        aligendImage = []
    end
    
    methods
        function objS = SNGSpotDetection(SpotNGliaObject) %constructor
            
            if exist('SpotNGliaObject', 'var')
                objS.SpotNGliaObject = SpotNGliaObject;
            end
            
            %set subfields of the registration field of the SngInputParameter
            %object as properties the obr i.e. SNGAllignment-object
            inputFields = fields(SpotNGliaObject.SngInputParameters.SpotDetection);
            for iInputField = 1:numel(inputFields)
                objS.(inputFields{iInputField}) = SpotNGliaObject.SngInputParameters.SpotDetection.(inputFields{iInputField});
            end
            
        end
        %         function set.imageSlice(obp, val)
        %             %automatically updating the property nSlices when updating imageSlice
        %             obp.imageSlice = val;
        %             obp.nSlices = numel(val);
        %         end
        
        
        function [SpotsDetected, SpotParameters, SpotDetectionInfo, SpotN] = SpotDetectionSNG(objS, Ialigned, CompleteTemplate, cmbr)
            
            
            Ialigned = objS.Ialigned;
            cmbr = objS.cmbr;
            
            
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
            
            
            %sng_zfinputAssign(ZFParameters,'SpotDetection')
            %obj.SngInputParameters.assign('SpotDetection')
            
            
            ScaleLevels = max(MPlevels);
            
            
            %TODO: suffisticated color transform of afterwards color detection
            %TODO: apply on original images which are not dept of field combined as :some spots appear due to dof artifacts. Combine images afterwards
            %TODO:  add more parameters
            %TODO:  find out how hard thresholding works
            
            
            %% RGB to gray algorithm
            
            GI = sng_RGB2Gray(Ialigned, ColorToGrayVector, false);
            
            
            %% Generate Wavelets
            
            %[MultiProductTh,MultiProduct,A,W,Wth] = sng_SpotWavelet(GI,ScaleLevels,ScaleBase,MPlevels,MPthreshold,false);
            [MultiProductTh, MultiProduct] = sng_SpotWavelet(GI, ScaleLevels, ScaleBase, MPlevels, MPthreshold, false);
            
            %{
  figure;imagesc(Ialigned)
  figure;imagesc(MultiProductTh)
 
  for k = 1:size(A,3)
  figure;imagesc(W(:,:,k))
  end
 
 
            %}
            
            %% SpotMeasure
            
            %function sng_spotproperties (Ialigned, MaskedIm)
            
            CC = bwconncomp(MultiProductTh);
            L = labelmatrix(CC);
            nspots1 = CC.NumObjects;
            Regions1 = regionprops(L, 'Area', 'PixelIdxList', 'PixelList', 'Centroid', 'Image', 'BoundingBox');
            
            
            %% Color mean of founded spots
            colormean = zeros(numel(Regions1), 3);
            for k = 1:numel(Regions1)
                spotimage = getfield(Regions1, {k}, 'Image');
                sbb = round(getfield(Regions1, {k}, 'BoundingBox')); %round because wc?
                spotcimage = Ialigned(sbb(2):sbb(2)+sbb(4)-1, sbb(1):sbb(1)+sbb(3)-1, 1:3);
                spotcimage(repmat(spotimage, 1, 1, 3) == 0) = 0; %color image of spot
                %color channels
                
                %TODO apply erosion, and solve isue with bouandary
                %erode spot for better colormean
                %CubeSizeSpotErode  =4
                %    SE1 = strel('cube',CubeSizeSpotErode); %spot structure element for erode
                %
                %    spotimage2 = imerode(spotimage,SE1)
                %
                %
                
                rc = spotcimage(:, :, 1);
                gc = spotcimage(:, :, 2);
                bc = spotcimage(:, :, 3);
                %spot color mean
                colormean(k, :) = [mean(rc(spotimage == 1)), mean(gc(spotimage == 1)), mean(bc(spotimage == 1))];
                %{
      figure;imagesc(spotimage);axis equal tight off;colormap gray
      figure;imagesc(spotcimage);axis equal tight off
                %}
            end
            temp = num2cell(colormean, 2); [Regions1.ColorMean] = temp{:};
            
            %color probability of mean spot color
            vec = round(reshape([Regions1.ColorMean], 3, numel(Regions1))');
            %correction index vs color
            %colorindex = sub2ind([255,255,255],vec(:,1),vec(:,2),vec(:,3));
            colorindex = sub2ind([256, 256, 256], vec(:, 1)+1, vec(:, 2)+1, vec(:, 3)+1);
            
            
            %ismember(uint32(colorindex),CompleteTemplate.SVAP_index)
            
            probability = zeros(numel(Regions1), 1);
            for k2 = 1:numel(Regions1)
                in = find(CompleteTemplate.SVAP_index == colorindex(k2));
                if ~isempty(in)
                    probability(k2, 1) = CompleteTemplate.SpotVectorArrayProbability(in);
                end
            end
            temp = num2cell(probability); [Regions1.ColorProbability] = temp{:};
            
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
            temp = num2cell(cm); [Regions1.MinProbability] = temp{:};
            
            
            %% determine if insite brain polygon
            [rc] = reshape([Regions1.Centroid], 2, nspots1)';
            [in, ~] = inpolygon(rc(:, 1), rc(:, 2), cmbr(:, 2), cmbr(:, 1));
            temp = num2cell(in); [Regions1.Insite] = temp{:};
            
            %{
  figure;
  imagesc(Ialigned)
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
    end
end

