classdef SNGMidBrainDetection < handle
    
    
    properties
        %Object
        SpotNGliaObject
        iFish
        
        %input variables
        Method
        Method2
        
        %TODO input variables which has to be added to SNGInp...
        rows = 100:10:400;
        
        
    end
    properties
        siz = [];
        rangex = [];
        rangey = [];
        
        edges = []
        DGraph = []
        path = []
        path2 = []
        path3 = []
        imageSize = []
        
        EdgeFilterWidth = []
        BrainArea = []
        BrainHorizontalsection = []
        BrainVerticalsection = []
        ShortestPath = []
        ShortestPathValue = []
        BrainEdge = []
        
        ShortestPathIndex     
        Area
        Hor
        Ver 
    end
    properties(Transient = true, Hidden = true)
        %data from template
        midbrainCenter
        edgeFilter
        polarMidbrainBandWithGaussianWithBlurr
        
        %
        AllShortestPaths
        
        %IMAGES
        alignedImage
        Ibrain
        
        %
        Isquare = uint8([])
        Ipolar
        ICorr2
        INorm
        INorm2
        Mask
        PolarTransform
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
                   
            objB.SpotNGliaObject = SpotNGliaObject;
        end
        
        %template getter functions
        function value = get.midbrainCenter(objB)
            if isempty(objB.midbrainCenter)
                objB.midbrainCenter = objB.SpotNGliaObject.CompleteTemplate.CenterMidBrain;
            end
            value = objB.midbrainCenter;
        end    
        function value = get.edgeFilter(objB)
            if isempty(objB.edgeFilter)
                objB.edgeFilter = objB.SpotNGliaObject.CompleteTemplate.EdgeFilter;
            end
            value = objB.edgeFilter;
        end   
        function value = get.polarMidbrainBandWithGaussianWithBlurr(objB)
            if isempty(objB.polarMidbrainBandWithGaussianWithBlurr)
                objB.polarMidbrainBandWithGaussianWithBlurr = objB.SpotNGliaObject.CompleteTemplate.polarMidbrainBandWithGaussianWithBlurr;
            end
            value = objB.polarMidbrainBandWithGaussianWithBlurr;
        end
        %image getter functions
        function value = get.alignedImage(objB)
            if isempty(objB.alignedImage)
               warning('First load the input image alignedImage with BrainSegmentationObject(fn).alignedImage = obj.RegObject(fn).Ialigned') 
               %a separate load function?
               %pros: it prevents long computing time as previous algorithms has not to be uploaded
               %cons: a separate function has to be called to load images
               %pros: maybe parallel computer is available
               objB.alignedImage = objB.SpotNGliaObject.RegObject(objB.iFish).Ialigned;
            end
            value = objB.alignedImage;
        end  
        function value = get.Isquare(objB)
            if isempty(objB.Isquare)
                objB.polarTransform
            end
            value = objB.Isquare;
        end
        function value = get.Ipolar(objB)
            if isempty(objB.Ipolar)
                objB.polarTransform
            end
            value = objB.Ipolar;
        end        
        function value = get.ICorr2(objB)
            if isempty(objB.ICorr2)
                objB.edgeCorrelation
            end
            value = objB.ICorr2;
        end    
        function value = get.INorm(objB)
            if isempty(objB.INorm)
                objB.edgeCorrelation
            end
            value = objB.INorm;
        end    
        function value = get.INorm2(objB)
            if isempty(objB.INorm2)
                objB.multiplyWithProbabilityMap
            end
            value = objB.INorm2;
        end    
        function value = get.PolarTransform(objB)
            if isempty(objB.PolarTransform)
                objB.findShortestPath
            end
            value = objB.PolarTransform;
        end
        function value = get.AllShortestPaths(objB)
            if isempty(objB.AllShortestPaths)
                objB.findShortestPath
            end
            value = objB.AllShortestPaths;
        end
        function value = get.Mask(objB)
            if isempty(objB.Mask)
                %objB.reversePolarTransformPath
                objB.Mask = sng_Polar2Im(objB.path3);
            end
            value = objB.Mask;
        end
        function value = get.Ibrain(objB)
            if isempty(objB.Mask)
                objB.reversePolarTransformPath;
            end
            value = objB.Ibrain;
        end
        
        function objB = All(objB)
            %output variable is needed for parfor in SpotNGlia
            polarTransform(objB)
            edgeCorrelation(objB)
            multiplyWithProbabilityMap(objB)
            findShortestPath(objB)
            reversePolarTransformPath(objB)
        end
        function polarTransform(objB)
            
            
            % polar transform
            %Isquare = objB.sng_cropAroundCoord(objB.alignedImage, objB.midbrainCenter);
            Isquare = sng_CropExtendAroundCoord(objB.alignedImage, objB.midbrainCenter, [1000 1000], 'mode');
            
            %extend with 200 pixels on every side
            objB.siz = size(objB.alignedImage);
            cxy2 = round(objB.midbrainCenter+200);
            objB.rangex = cxy2(1) - 499:cxy2(1) + 500;
            objB.rangey = cxy2(2) - 499:cxy2(2) + 500;
            
            %{
          sng_show(Isquare)
            %}
            
            %polar coordinates
            Ipolar = sng_Im2Polar3(Isquare);
            objB.Ipolar = permute(Ipolar, [2, 1, 3]);
            
            objB.Isquare = Isquare;
            
            %{
          %view(gca,[180 90]);
          sng_show(Ipolar)
            %}
            
        end
        function edgeCorrelation(objB)
            
            %input parameters
            %Method = objB.Method;
            %Method2 = objB.Method;
            
            EdgeFilter = objB.edgeFilter;
            
            %TODO: normalized or not normalized correlation?
            
            Ipolar = objB.Ipolar;
            
            % correlation with Edgetemplates
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
            objB.ICorr2 = ICorr2;
            objB.INorm = INorm;
            objB.imageSize = si;
        end
        function multiplyWithProbabilityMap(objB)
            % Filter image with brain probability
            INorm = objB.INorm;
            
            %INorm2 = INorm .* IDistmap4;
            INorm2 = INorm .* objB.polarMidbrainBandWithGaussianWithBlurr;
            objB.INorm2 = INorm2;
            %{
          figure;imagesc(IDistmap);sng_imfix;colormap gray
          figure;imagesc(IDistmap2);sng_imfix;colormap gray
          figure;imagesc(IDistmap3);sng_imfix;colormap gray
          figure;imagesc(IDistmap4);sng_imfix;colormap gray
 
          figure;imagesc(INorm2);sng_imfix;
            %}
        end
        function findShortestPath(objB)
            
            % Find shortest Path
            si = objB.imageSize;
            INorm2 = objB.INorm2;
            rows = objB.rows;
            
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
                AllShortestPaths{l} = [J', I'];
            end
            
            %{
              figure;plot(Area)
              figure;plot(Hor)
              figure;plot(Ver)
              figure;plot(d)
            %}
            
            % find lowestcost path from all lowest cost past with equal end and start row
            %TODO:  apply area/horizontal/vertical distance brain constrained
            %       find max value of hor and ver
            [mx, mr] = min(lengthd);
            [path, ~] = shortestpath(DGraph, a(mr), b(mr), 'Method', 'acyclic');
            [I, J] = ind2sub(s, path);
            
            %{
          figure;imagesc(INorm3);sng_imfix;hold on;plot(J,I,'r','LineWidth',2)
          figure;imagesc(INorm2);sng_imfix;hold on;plot(J,I,'r','LineWidth',2)
            %}
            
            objB.ShortestPath = [J', I'];
            objB.ShortestPathValue = mx;
            objB.ShortestPathIndex = mx;
            objB.PolarTransform = INorm3;
            onjB.edges = edges;
            onjB.DGraph = DGraph;
            objB.path = path;
            objB.Area = Area;
            objB.Hor = Hor;
            objB.Ver = Ver;
            objB.AllShortestPaths = AllShortestPaths;
            
        end
        function reversePolarTransformPath(objB)
            
            I = objB.ShortestPath(:,2)';
            siz = objB.siz;
            rangex = objB.rangex;
            rangey = objB.rangey;
            
            s = size(objB.PolarTransform);
            path = objB.path;
            
            % create mask in cartesian space
            path2 = zeros(s);
            path2(path) = 1;
            %make everything below the path one
            for kk = 1:size(path2, 2)
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
              figure;imagesc(uint8(alignedImage))
 
              figure;imagesc(uint8(Ibrain))
              figure;plot(J,I)
 
              IA = double(alignedImage);
              IA(:,:,3) = IA(:,:,3) + 110 * double(IBrain);
              IA(:,:,2) = IA(:,:,2) + 50 * double(IBrain);
              figure;imagesc(uint8(IA));sng_imfix
 
          bc = bwboundaries(Ibrain)
          figure;imagesc(uint8(alignedImage));sng_imfix
          hold on;plot(bc{1}(:,2),bc{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
 
            %}
            
            
            bwbb = bwboundaries(Ibrain);
            
            %
            
            %objB.EdgeFilterWidth = filterwidth;
            %objB.BrainArea = BrainArea;
            %objB.BrainHorizontalsection = BrainHorizontalsection;
            %objB.BrainVerticalsection = BrainVerticalsection;
            
            objB.BrainEdge = bwbb{1};
            objB.Mask = Mask;
            objB.Ibrain = Ibrain;
            objB.path2 = path2;
            objB.path3 = path3;
            
            
            % the fish movie (not finished, look for mov_polartransform3)
            %{
      siz = size(alignedImage)
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
    
    methods(Static = true)
        %{
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
            
             Img1 = Image(:,:,5);
             coords1 = obj.BrainSegmentationObject(iFish).midbrainCenter
 
            %}
            siz = size(Img1);
            if ndims(Img1) == 2
                siz(3) = 1;
            end
            
            if ~exist('exteriorColor')
                exteriorColor = 255;
            end
            ext = (exteriorColor * ones(siz(1)+400, siz(2)+400, siz(3)));
            %ext = uint8(exteriorColor * ones(siz(1)+400, siz(2)+400, siz(3)));
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
        %}
        function [edges] = sng_SquareGridDiGraph(s)
            %this function generates nodes and edges from a rectangular grid (matrix)
            %direction is from left to right and connectivety is -45,0,45 degrees
            %also the values from the matrix input in terms of nodes are calculated
            %   s:  size of 2 dimensional input matrix
            
            
            %!!!fix the starting point
            
            %{
   ex
   matrix = rand(100,100);
   [indcom val] = sng_SquareGridDiGraph(matrix)
 
            %}
            
            %size
            %s = size(matrix)
            %indices
            [CC, RR] = meshgrid(1:s(2), 1:s(1));
            
            %RR = int32(RR);
            %CC = int32(CC);
            
            
            %{
   figure;scatter(CC(:),RR(:));axis ij
            %}
            
            %% first and last row
            RRfr = RR(1, 1:end-1);
            RRlr = RR(end, 1:end-1);
            CCfr = CC(1, 1:end-1);
            CClr = CC(end, 1:end-1);
            %edges of first row
            subfr = [RRfr(:), CCfr(:)];
            efr1 = subfr + repmat([0, 1], [length(subfr), 1]); %one step right
            efr2 = subfr + repmat([1, 1], [length(subfr), 1]); %one step down and right
            %edges of last row
            sublr = [RRlr(:), CClr(:)];
            elr1 = sublr + repmat([0, 1], [length(sublr), 1]); %one step right
            elr2 = sublr + repmat([-1, 1], [length(sublr), 1]); %one step up and right
            %sub to ind first row
            indfr = sub2ind(s, subfr(:, 1), subfr(:, 2));
            indefr1 = sub2ind(s, efr1(:, 1), efr1(:, 2));
            indefr2 = sub2ind(s, efr2(:, 1), efr2(:, 2));
            %sub to in last row
            indlr = sub2ind(s, sublr(:, 1), sublr(:, 2));
            indelr1 = sub2ind(s, elr1(:, 1), elr1(:, 2));
            indelr2 = sub2ind(s, elr2(:, 1), elr2(:, 2));
            %edges of first en last row
            
            %% indices except first and last row and last colum
            RRmr = RR(2:end-1, 1:end-1);
            CCmr = CC(2:end-1, 1:end-1);
            %indices middle rows
            submr = [RRmr(:), CCmr(:)];
            %edges
            e1 = submr + repmat([-1, 1], [length(submr), 1]);
            e2 = submr + repmat([0, 1], [length(submr), 1]);
            e3 = submr + repmat([1, 1], [length(submr), 1]);
            %indices of edges
            ind = sub2ind(s, submr(:, 1), submr(:, 2));
            inde1 = sub2ind([s(1), s(2)], e1(:, 1), e1(:, 2));
            inde2 = sub2ind([s(1), s(2)], e2(:, 1), e2(:, 2));
            inde3 = sub2ind([s(1), s(2)], e3(:, 1), e3(:, 2));
            
            %complete edges
            edges = [[indfr, indefr1]; [indfr, indefr2]; [indlr, indelr1]; [indlr, indelr2]; [ind, inde1]; [ind, inde2]; [ind, inde3]];
            %all values
            
            %values = matrix(edges(:,2));
            
            
            %{
   %the adjacentcy matrix doesn have to be determined as digraph generates a
   spart matrix
   adj = zeros(24)
   for j = 1:size(indc,1)
       adj(indc(j,1),indc(j,2)) = val(j)
   end
   %this does the same as accumarray
            %}
            
            %{
   %in addition the adjacentcy matrix or the graph cal be computer from the
   node and edges
 
   adj = accumarray([indcom;indflc],val,[24 24]) %adjacentcy matrix
   adj = accumarray(indcom,val,[prod(s),prod(s)],[],0,true); %sparse
 
 
   G3 = digraph([indcom(:,1)],[indcom(:,2)],val);
 
            %}
            
            
        end
    end
end
