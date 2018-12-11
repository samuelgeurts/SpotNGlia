classdef SNGMidBrainDetection < handle
    
    
    properties
        
        %input variables
        Method
        Method2
        
        %from template
        midbrainCenter
        
        edgeFilter
        
        SpotNGliaObject
    end
    
    properties(Transient = true)
        Alingenfish
        
    end
    
    
    methods
        function objB = SNGMidBrainDetection(SpotNGliaObject) %constructor
            %obr = obr.zfproperties(SpotNGliaObject.ZFParameters);
            %obj.SngInputParameters.assign('BrainSegmentation')
            %objB.template = SpotNGliaObject.CompleteTemplate.Template; ???
            
            %extract brain segmentation parameters from spotnglia object
            inputFields = fields(SpotNGliaObject.SngInputParameters.BrainSegmentation);
            for iInputField = 1:numel(inputFields)
                objB.(inputFields{iInputField}) = SpotNGliaObject.SngInputParameters.BrainSegmentation.(inputFields{iInputField});
            end
            
            %             %loads complete template, maybe not necesary
            %             if isempty(SpotNGliaObject.CompleteTemplate)
            %                 SpotNGliaObject = SpotNGliaObject.LoadTemplate;
            %             end
            
            objB.edgeFilter = SpotNGliaObject.CompleteTemplate.EdgeFilter;
            objB.SpotNGliaObject = SpotNGliaObject;
            objB.centerMidbrain = CompleteTemplate.CenterMidBrain;
        end
        
        
        function [Ibrain, Parameters] = MidBrainDetectionSNG(AlignedFish, CompleteTemplate, ZFParameters)
            %Detect midbrain
            %works in combination with other "Link" functions
            %based on MidBraindetection 2 and GeneratePolarTransforms2
            
            %Version MidBrainDetectionLink2
            %   applied on template aligned fishes
            %   add skeletonband for path selection
            %   put 1000x1000 selection in other function (sng_boxaroundcenter)
            %Version MidBrainDetectionLink3
            %   applicable for ZFParameters and broutput
            %Version MidBrainDetectionSNG
            %   applicable for SpotNGlia
            %20181105 Version SNGMidBrainDetection
            %   made a class added sng_cropAroundCoord
  
            
            %input parameters
            Method = objB.Method;
            Method2 = objB.Method;
            midbrainCenter = objB.midbrainCenter;
            
            
            %% polar transform
            [Isquare] = objB.sng_cropAroundCoord(AlignedFish, midbrainCenter);
            
            %extend with 200 pixels on every side
            siz = size(AlignedFish);
            cxy2 = round(midbrainCenter+200);
            rangex = cxy2(1) - 499:cxy2(1) + 500;
            rangey = cxy2(2) - 499:cxy2(2) + 500;
            
            %{
      sng_show(Isquare)
            %}
            
            %polar coordinates
            Ipolar = sng_Im2Polar3(Isquare);
            Ipolar = permute(Ipolar, [2, 1, 3]);
            
            %{
      %view(gca,[180 90]);
      sng_show(Ipolar)
            %}
            
            %TODO: normalized or not normalized correlation?
            
            %% correlation with Edgetemplates
            ICorr(:, :, 1) = normxcorr2_general(EdgeFilter(:, :, 1), Ipolar(:, :, 1), numel(EdgeFilter(:, :, 1))/2);
            ICorr(:, :, 2) = normxcorr2_general(EdgeFilter(:, :, 2), Ipolar(:, :, 2), numel(EdgeFilter(:, :, 2))/2);
            ICorr(:, :, 3) = normxcorr2_general(EdgeFilter(:, :, 3), Ipolar(:, :, 3), numel(EdgeFilter(:, :, 3))/2);
            
            %{
      figure;imagesc((ICorr+1)/2);sng_imfix
            %}
            
            %crop the correlated image to the original size
            stf = floor(size(EdgeFilter(:, :, 1))/2); %start coordinates to add
            stc = ceil(size(EdgeFilter(:, :, 1))/2); %end coordinates to add
            si = size(Ipolar); %image size
            crR = stc(1):si(1) + stf(1); %row coordinates
            crC = stc(2):si(2) + stf(2); %column coordinates
            ICorr2 = ICorr(stc(1):si(1)+stf(1), stc(2):si(2)+stf(2), :);
            
            %remove values below zero or lift above zero
            %ICorr2 (ICorr2 < 0) = 0;
            ICorr2 = (ICorr2 + 1) / 2; %function lift
            
            
            %{
      figure;imagesc(ICorr2);sng_imfix
            %}
            
            %combine channels
            %ISum = sum(ICorr2,3);
            INorm = sqrt(sum(ICorr2.^2, 3));
            
            %{
      figure;imshow(ISum)
      figure;imagesc(INorm);colormap gray;sng_imfix
      figure;imagesc(INorm);sng_imfix
            %}
            
            %% Polar transform distancemap
            %computation moved to complete template as it has to be don once
            
            
            %{
      %DONE: use more suffisticated method to exclude area in shortest path
      %finding, multiply shortest path with based on brains blurred distancemap
      %to exclude some area from pathfinding
      %TODO: add band to generateskeleton function output instead of computing from distancemap
      IsquareMB = sng_boxaroundcenter(CompleteTemplate.MidBrainDistanceMap,cxy);
 
            %{
      figure;imagesc(boolean(IsquareMB));sng_imfix;colormap gray
            %}
 
      IDistmap = sng_Im2Polar3(IsquareMB);
      IDistmap = permute(IDistmap,[2,1,3]);
      IDistmap(IDistmap ~= 0) = 1; %set every nonzero value 2 one (band)
 
      %gaussion blur on IDistmap, probably doesnt work as the main plan is
      %blurred
      IDistmap2 = imgaussfilt(IDistmap,[40,5]);
 
      %IDistmap3 contains added gausian filter at te boudary i.e. the values of
      %IDistmap which are 1 stays 1, i.e. the plane is preserved
      blurvar = 6;
      IDistmap3 = IDistmap;
      sz = size(IDistmap3);
      center = 501;
      for k2 = 1:sz(2)
      a = find(IDistmap3(:,k2) == 1,1,'first');
      b = find(IDistmap3(:,k2) == 1,1,'last');
      %y = gaussmf([-500:1:500],[(b-a)/blurvar 0])';
      %for versions without fuzzy logic toolbox
      y = exp(-((-500:1:500)-0).^2/(2*((b-a)/blurvar).^2))';
 
      IDistmap3(1:a,k2) = y(center-a+1:center);
      IDistmap3(b:sz(1),k2) = y(center:center+sz(1)-b);
      end
      %a bit of gaussian blur added
      IDistmap4 = imgaussfilt(IDistmap3,[5,5]);
      %figure;plot(y) % a gaussian
            %}
            
            
            %% Filter image with brain probability
            
            %INorm2 = INorm .* IDistmap4;
            INorm2 = INorm .* CompleteTemplate.polarMidbrainBandWithGaussian;
            
            %{
      figure;imagesc(IDistmap);sng_imfix;colormap gray
      figure;imagesc(IDistmap2);sng_imfix;colormap gray
      figure;imagesc(IDistmap3);sng_imfix;colormap gray
      figure;imagesc(IDistmap4);sng_imfix;colormap gray
 
      figure;imagesc(INorm2);sng_imfix;
            %}
            
            %% Find shortest Path
            
            %invert for lowest cost path and normalize
            INorm3 = -INorm2 / max(INorm2(:)) + 1;
            
            %{
      figure;imagesc(INorm3);sng_imfix;
            %}
            
            
            s = size(INorm3);
            %compute edges
            [edges] = sng_SquareGridDiGraph(si(1:2));
            %compute directer graph
            DGraph = digraph([edges(:, 1)], [edges(:, 2)], INorm3(edges(:, 2)));
            
            
            %compute shortest path start and end at equal row
            rows = 100:10:400;
            nrows = numel(rows);
            path = zeros(nrows, 1000);
            
            
            %TODO: check if area calculation is right
            %figure;imagesc(INorm3);sng_imfix
            a = arrayfun(@(x) sub2ind(si, x, 1), rows);
            b = arrayfun(@(x) sub2ind(si, x, s(2)), rows);
            
            for l = 1:nrows
                l
                [path, lengthd(l)] = shortestpath(DGraph, a(l), b(l), 'Method', 'acyclic');
                [I, J] = ind2sub(s, path);
                
                Area(l) = sum(pi*(I.^2)/1000);
                Hor(l) = I(end) + I(end/2);
                Ver(l) = I(round((3 * s(2))/4)) + I(round((s(2))/4));
                
                %hold on
                %plot(J,I,'r','LineWidth',2)
                %drawnow
            end
            
            %{
          figure;plot(Area)
          figure;plot(Hor)
          figure;plot(Ver)
          figure;plot(d)
            %}
            
            %% find lowestcost path from all lowest cost past with equal end and start row
            %TODO:  apply area/horizontal/vertical distance brain constrained
            %       find max value of hor and ver
            [mx, mr] = min(lengthd)
            [path, ~] = shortestpath(DGraph, a(mr), b(mr), 'Method', 'acyclic');
            [I, J] = ind2sub(s, path);
            
            %{
      figure;imagesc(INorm3);sng_imfix;hold on;plot(J,I,'r','LineWidth',2)
      figure;imagesc(INorm2);sng_imfix;hold on;plot(J,I,'r','LineWidth',2)
            %}
            
            
            %% create mask in cartesian space
            path2 = zeros(s);
            path2(path) = 1;
            %make everything below the path one
            for kk = 1:size(path2, 2);
                path2(1:I(kk), kk) = 1;
            end
            
            %{
          figure;imagesc(path2);sng_imfix
            %}
            
            
            path3 = permute(path2, [2, 1]); %restore original permutation
            Mask = sng_Polar2Im(path3);
            
            %{
          figure;imagesc(path3)
          figure;imagesc(uint8(Mask));sng_imfix
 
            %}
            
            frame = false(siz(1:2)+400);
            frame(rangey, rangex) = Mask;
            Ibrain = frame(201:end-200, 201:end-200);
            
            %{
          figure;imagesc(uint8(AlignedFish))
 
          figure;imagesc(uint8(Ibrain))
          figure;plot(J,I)
 
          IA = double(AlignedFish);
          IA(:,:,3) = IA(:,:,3) + 110 * double(IBrain);
          IA(:,:,2) = IA(:,:,2) + 50 * double(IBrain);
          figure;imagesc(uint8(IA));sng_imfix
 
      bc = bwboundaries(Ibrain)
      figure;imagesc(uint8(AlignedFish));sng_imfix
      hold on;plot(bc{1}(:,2),bc{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
 
            %}
            
            
            bwbb = bwboundaries(Ibrain);
            
            %%
            
            Parameters.EdgeFilterWidth = filterwidth;
            %Parameter.BrainArea
            %Parameter.BrainHorizontalsection
            %Parameter.BrainVerticalsection
            Parameters.ShortestPath = [J', I'];
            Parameters.ShortestPathValue = mx;
            
            Parameters.BrainEdge = bwbb{1};
            
            Parameters.PolarTransform = INorm3;
            
            
            %% the fish movie (not finished, look for mov_polartransform3)
            %{
  siz = size(AlignedFish)
  cxy2 = round(cxy + 200)
  rangex = cxy2(1)-499:cxy2(1)+500;
  rangey = cxy2(2)-499:cxy2(2)+500;
 
 
  figure;imagesc(uint8(Isquare));axis equal tight off
 
 
  tf = 200 %totalframes
  scale = 1
  viewrot = 90
  cutrot = 90
 
  Img = Isquare;
  Img = imrotate(Img,cutrot);
 
  h= figure(1);
 
  R=Img(:,:,1);
  G=Img(:,:,2);
  B=Img(:,:,3);
 
  s = size(Img);
  [X,Y] = meshgrid(1:s(2),1:s(1));
 
  %make center zero
  X2 = X - cxy2(2);
  Y2 = Y - cxy2(1);
 
  %image to polar
  Radius = sqrt(X2.^2+Y2.^2);
  Theta = atan2(Y2,X2); %branchcut on the left -pi -> pi
  r = scale * 500;%readius of circle
 
  %make area outsite radius white
  R2=double(R);R2(Radius >= r) = 256;
  G2=double(G);G2(Radius >= r) = 256;
  B2=double(B);B2(Radius >= r) = 256;
 
  writerObj1 = VideoWriter('polar tranform2','MPEG-4');
  writerObj1.FrameRate = 30;    % to perform realtime movie
  open(writerObj1);
 
  %the transformation parameters. s is for the angle transformation.
  s = linspace(1,0.01,tf).^2; %final s should be small but nog 0
  % s2 isfor the radius transformation based on the values of s
  s2 = (linspace(0,10,tf).^2)./s
  %used for the x-limit
  sxmin = linspace(1,0,tf).^2;
 
 
  HX = Radius .* cos(Theta);
  HY = Radius .* sin(Theta);
  figure(1);
  ScatH = scatter(HX(:),HY(:),100,double(cat(2,R2(:),G2(:),B2(:)))/255,'.','LineWidth',2);
  view(gca,[viewrot 90]);
  axis tight off equal
 
  for k = 1:numel(s)
 
  HX = (Radius+s2(k)) .* cos(Theta*s(k));
  HY = (Radius+s2(k)) .* sin(Theta*s(k));
 
  set(ScatH,'Xdata',HX(:))
  set(ScatH,'Ydata',HY(:))
 
  %xlimit such that the initial center of the image is smootly shifted to the
  %right. if k = 1, sxmin(1)=1 and xmin = 2*HX(end/2,end/2) - max(HX(R2 ~= 255)));
  % if k = 100, sxmin(1)=0 and xmin = HX(end/2,end/2), quadrating makes it smoother
  xmin = HX(floor(end/2),floor(end/2)) + sxmin(k)*(HX(floor(end/2),floor(end/2))- max(HX(R2 ~= 256)));
  xlim([xmin,max(HX(R2 ~= 256))]);
  ylim([min(HY(R2 ~= 256)),max(HY(R2 ~= 256))]);
 
  axis square off
 
 
  %xlim([xmin,HX(floor(end/2),end)]);
  %min(HX(R2~=256))
 
 
  %pause (0,1)
 
  frame = getframe(gcf);
  writeVideo(writerObj1,frame);
  end
 
  close(writerObj1);
 
  %% from fish to circel in fish
  k=1
  HX = (Radius+s2(k)) .* cos(Theta*s(k));
  HY = (Radius+s2(k)) .* sin(Theta*s(k));
 
  %inital fish
  figure;scatter(HX(:),HY(:),100,double(cat(2,R(:),G(:),B(:)))/255,'.','LineWidth',2);
  axis equal tight off;
  view(gca,[viewrot 90]);
 
  %with circle
  figure;scatter(HX(:),HY(:),100,double(cat(2,R(:),G(:),B(:)))/255,'.','LineWidth',2);
  axis equal tight off
  pos = [-500,-500, 1000, 1000]*scale;
  rectangle('Position',pos,'Curvature',[1 1]);
  view(gca,[viewrot 90]);
 
  %remove out circle
  figure;scatter(HX(:),HY(:),100,double(cat(2,R2(:),G2(:),B2(:)))/255,'.','LineWidth',2);
  axis equal tight off
  pos = [-500,-500, 1000, 1000]*scale;
  rectangle('Position',pos,'Curvature',[1 1]);
  view(gca,[viewrot 90]);
 
  %scale movie
  figure;scatter(HX(:),HY(:),100,double(cat(2,R2(:),G2(:),B2(:)))/255,'.','LineWidth',2);
  axis equal tight off
  pos = [-500,-500, 1000, 1000]*scale;
  rectangle('Position',pos,'Curvature',[1 1]);
  view(gca,[viewrot 90]);
 
 
  writerObj2 = VideoWriter('polar tranform zoom','MPEG-4');
  writerObj2.FrameRate = 30;    % to perform realtime movie
  open(writerObj2);
 
  z = linspace(750*scale,500*scale,50);
  xlim([-z(1),z(1)])
  for j = 1:numel(z)
  xlim([-z(j),z(j)])
  ylim([-500*scale,500*scale])
 
 
  frame = getframe(gcf);
  writeVideo(writerObj2,frame);
  end
 
  close(writerObj2);
 
            %}
            
        end
    end
    
    methods (Static = true)
        function [Img2] = sng_cropAroundCoord(Img1, coords1, boxsize, exteriorColor)
            %extends images and cut out a 1000x1000 image around a given center
            %extend with 200 pixels on every side
            
            %20181105 Version sng_cropAroundCoord comes from the function
                %sng_boxaroundcenter2. As it is a very specific function, it is
                %changed to a static function under the SNGMidBrinDetection class
            
            %{
         Img1 = CompleteTemplate.MidBrainDistanceMap;
         coords1 = cxy
 
         Img1 = objt.BandMidBrain;
         coords1 = objt.CenterMidBrain;
 
            %}
            siz = size(Img1);
            if ndims(Img1) == 2
                siz(3) = 1;
            end
            
            if ~exist('exteriorColor')
                exteriorColor = 255;
            end
            
            ext = exteriorColor * ones(siz(1)+400, siz(2)+400, siz(3));
            ext(201:siz(1)+200, 201:siz(2)+200, :) = Img1;
            
            
            %{
         figure;imagesc(uint8(ext))
         hold on
         scatter(cxy2(1),cxy2(2))
            %}
            
            
            %square 1000x1000 images with skeleton based center in the middle
            cxy2 = round(coords1+200);
            rangex = cxy2(1) - 499:cxy2(1) + 500;
            rangey = cxy2(2) - 499:cxy2(2) + 500;
            Img2 = ext(rangey, rangex, :);
            
            %{
         figure;imagesc(uint8(Img2))
            %}
            
            
            %{
         Img3 = ones(1000,1000)*255;
         c1 = round(coords1)
         r1 = max([c1 - 499;1,1])
         r2 = min([c1 + 500;siz(1),siz(2)])
 
         c1-499   substract the first part
 
         Img3(r1(1):r2(1),r1(2):r2(2)) = Img1(r1(1):r2(1),r1(2):r2(2))
 
         figure;imagesc(uint8(Img3))
            %}
        end
        
        
    end
    
    
end
