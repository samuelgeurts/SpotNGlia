classdef SNGPreprocessing < handle
    %extern functions needed
    %sng_smoothzoom -for movie
    %sng_zoom
    %sng_NCC
    %IAT toolbox
    
    properties %(Access = private)
        %Objtect
        SpotNGliaObject
        %INPUTPARAMETERS
        onoff = []; %sets bandpassfilter on or off
        sigmalp = [];
        sigmahp = [];
        scaleC = [];
        levelsC = [];
        iterationsC = [];
        scaleS = [];
        levelsS = [];
        iterationsS = [];
        %extended dept of field
        variancedisksize = []; %sets bandpassfilter on or off
        
        %from CompleteTemplate
        template = []
        
        
        %ShowMovieProperties
        magnification = [1, 20];
        zoomlocation = [459, 768];
        nZoomInSteps = 60;
        nZoomOutSteps = 40;
        sliceNo = 1;
    end
    properties
        %OUTPUTPARAMETERS
        nSlices
        %imageFullPath
        iFish
        imageNames
        imagePath
        
        %rgbCorrection
        colorWarp; % = ECCWarp;
        NCC_Img_before;
        NCC_Img_after;
        NCC_BP_before;
        NCC_BP_after;
        rgbCropValues;
        
        %stackCorrection
        sliceWarp; % = ECCWarp2;
        
        %extended dept of field
        IndexMatrix = []
        
    end
    properties(Transient = true, Hidden = true)
        %IMAGES
        imageSlice = []
        imageSliceCor = []
        correctedSlice = []
        %RgbCorrection
        bandFiltered = []
        filteredImage = []
        
        %extended dept of field
        variance_sq = []
        mergedImage = []
    end
    
    methods
        function obp = SNGPreprocessing(SpotNGliaObject) %constructor
            
            if exist('SpotNGliaObject', 'var')
                obp.SpotNGliaObject = SpotNGliaObject;
            end
            
            %set subfields of the registration field of the SngInputParameter
            %object as properties the obr i.e. SNGAllignment-object
            inputFields = fields(SpotNGliaObject.SngInputParameters.Preprocessing);
            for iInputField = 1:numel(inputFields)
                obp.(inputFields{iInputField}) = SpotNGliaObject.SngInputParameters.Preprocessing.(inputFields{iInputField});
            end
            
        end
        function set.imageSlice(obp, val)
            %automatically updating the property nSlices when updating imageSlice
            obp.imageSlice = val;
            obp.nSlices = numel(val);
        end
        
        %{
        function obp = loadImageSlices(obp) %Load slices
            obp.imageSlice = cell(1, obp.nSlices); %preallocate for every new slice
            for iSlice = 1:obp.nSlices
                obp.imageSlice{iSlice} = imread([obp.imagePath, filesep, obp.imageNames{iSlice}]);
                obp.imageSlice{iSlice} = im2uint8(obp.imageSlice{iSlice}(:, :, 1:3));
            end
        end
        %}
        
        function value = get.imageSlice(obp)
            if isempty(obp.imageSlice)
                obp.imageSlice = cell(1, obp.nSlices); %preallocate for every new slice
                for iSlice = 1:obp.nSlices
                    obp.imageSlice{iSlice} = imread([obp.imagePath, filesep, obp.imageNames{iSlice}]);
                    obp.imageSlice{iSlice} = im2uint8(obp.imageSlice{iSlice}(:, :, 1:3));
                end
                disp('load images')  
            end
            value = obp.imageSlice;
        end        
        function value = get.imageSliceCor(obp)
            
            if isempty(obp.imageSliceCor)
                rgbWarp(obp);
                disp('apply rgb alignment on all images of a single fish')
            end
            value = obp.imageSliceCor;
        end
        function value = get.correctedSlice(obp)
            if isempty(obp.correctedSlice)
                obp.correctedSlice = obp.sng_stackWarp(obp.imageSliceCor, obp.sliceWarp);
                disp('apply stack alignment')
            end
            value = obp.correctedSlice;
        end
        function value = get.mergedImage(obp)
            if isempty(obp.mergedImage)
                obp.mergedImage = obp.sng_sliceCombine(obp.correctedSlice, obp.IndexMatrix);
                disp('apply image merge (extended dept of field)')
            end
            value = obp.mergedImage;
        end
        
        function obp = rgbCorrection(obp)
            %computes rgb correction i.e. translation of color channels
            
            %preallocation
            obp.colorWarp = cell(1, obp.nSlices);
            
            
            %does RGB correction on all images of one fish
            for iSlice = 1:obp.nSlices
                
                if obp.onoff
                    if isempty(obp.filteredImage)
                        obp = rgbBandpassFilter(obp);
                    end
                    
                    %computes warp on filtered image (green en blue channel)
                    [obp.colorWarp{iSlice}, ~] = obp.sng_rgbAlignment(obp.filteredImage{iSlice}(:, :, 1:3), obp.scaleC, 'translation', obp.levelsC, obp.iterationsC);
                else
                    %computes warp on unfiltered image (green en blue channel)
                    [obp.colorWarp{iSlice}, ~] = obp.sng_rgbAlignment(obp.imageSlice{iSlice}(:, :, 1:3), obp.scaleC, 'translation', obp.levelsC, obp.iterationsC);
                end
                
                %{
                             figure;imagesc((obp.filteredImage{iSlice}(:,:,1:3)));axis equal off tight
                             figure;imagesc(uint8(obp.imageSlice{iSlice}(:,:,1:3)));axis equal off tight
                %}
            end
        end
        function obp = rgbBandpassFilter(obp)
            
            %preallocation
            obp.filteredImage = cell(1, obp.nSlices);
            
            for iSlice = 1:obp.nSlices
                %bandpass filter focussed on spotsize passing
                [obp.filteredImage{iSlice}, obp.bandFiltered{iSlice}] = obp.sng_bandpassFilter(obp.imageSlice{iSlice}(:, :, 1:3), obp.sigmalp, obp.sigmahp); %%TODO put this function in this class
                
            end
        end
        function obp = rgbWarp(obp)
            
            %preallocation
            obp.NCC_Img_before = cell(1, obp.nSlices);
            obp.NCC_Img_after = cell(1, obp.nSlices);
            obp.imageSliceCor = cell(1, obp.nSlices);
            obp.rgbCropValues = cell(1, obp.nSlices);
            
            for iSlice = 1:obp.nSlices
                
                %colorwarp
                [obp.imageSliceCor{iSlice}, obp.rgbCropValues{iSlice}] = obp.sng_rgbWarp(obp.imageSlice{iSlice}(:, :, 1:3), obp.colorWarp{iSlice});
                
                %compute correlation
                %{
                             %normal correlation
                             CC1 = corr2(filteredImage(:,:,1),filteredImage(:,:,2));
                             CC2 = corr2(filteredImage2(:,:,1),filteredImage2(:,:,2));
                             CC3 = corr2(imageSlice{iSlice}(:,:,1),imageSlice{iSlice}(:,:,3));
                             CC4 = corr2(imageSliceCor{iSlice}(:,:,1),imageSliceCor{iSlice}(:,:,3));
                %}
                %original image before and after normalised correlation (green-red and blue-red channel)
                obp.NCC_Img_before{iSlice}(1) = sng_NCC(obp.imageSlice{iSlice}(:, :, 1), obp.imageSlice{iSlice}(:, :, 2));
                obp.NCC_Img_before{iSlice}(2) = sng_NCC(obp.imageSlice{iSlice}(:, :, 1), obp.imageSlice{iSlice}(:, :, 3));
                obp.NCC_Img_after{iSlice}(1) = sng_NCC(obp.imageSliceCor{iSlice}(:, :, 1), obp.imageSliceCor{iSlice}(:, :, 2));
                obp.NCC_Img_after{iSlice}(2) = sng_NCC(obp.imageSliceCor{iSlice}(:, :, 1), obp.imageSliceCor{iSlice}(:, :, 3));
            end
            %imwrite(uint8(Img2), fullFileName{l}, 'WriteMode', 'append', 'Compression','none');
        end
        function obp = rgbWarpFilteredImage(obp)
            
            %warps the bandpass filterd image which is only use for analysis
            if isempty(obp.filteredImage)
                obp = rgbBandpassFilter(obp);
            end
            
            obp.NCC_BP_before = cell(1, obp.nSlices);
            obp.NCC_BP_after = cell(1, obp.nSlices);
            
            for iSlice = 1:obp.nSlices
                
                [CorrectedfilteredImage, ~] = obp.sng_rgbWarp(obp.filteredImage, obp.colorWarp{iSlice});
                
                %{
                               figure;imagesc((imageSliceCor{iSlice}(:,:,1:3)));axis equal off tight
                               figure;imagesc((CorrectedfilteredImage(:,:,1:3)));axis equal off tight
                %}
                
                %bandpass filtered before and after normalised correlation  (green-red and blue-red channel)
                obp.NCC_BP_before{iSlice}(1) = sng_NCC(obp.filteredImage{iSlice}(:, :, 1), obp.filteredImage{iSlice}(:, :, 2));
                obp.NCC_BP_before{iSlice}(2) = sng_NCC(obp.filteredImage{iSlice}(:, :, 1), obp.filteredImage{iSlice}(:, :, 3));
                obp.NCC_BP_after{iSlice}(1) = sng_NCC(CorrectedfilteredImage(:, :, 1), CorrectedfilteredImage(:, :, 2));
                obp.NCC_BP_after{iSlice}(2) = sng_NCC(CorrectedfilteredImage(:, :, 1), CorrectedfilteredImage(:, :, 3));
            end
        end
        %{
function obp = stackCorrection(obp)
             obp.sliceWarp = obp.sng_stackAlignment(obp.imageSliceCor, obp.scaleS, 'translation', obp.levelsS, obp.iterationsS);
             %[obp.sliceWarp, obp.correctedSlice] = obp.sng_stackAlignment(obp.imageSliceCor, obp.scaleS, 'translation', obp.levelsS, obp.iterationsS);
 
         end
function obp = stackWarp(obp)
 
               [obp.correctedSlice] = obp.sng_stackWarp(obp.imageSliceCor, obp.sliceWarp);
         end
function obp = ExtendedDeptofField(obp)
             [obp.IndexMatrix, obp.variance_sq] = obp.sng_StackDOF2(obp.correctedSlice, obp.variancedisksize);
             %[obp.IndexMatrix, obp.variance_sq, obp.mergedImage] = obp.sng_StackDOF2(obp.correctedSlice, obp.variancedisksize);
             obp.mergedImage = sng_SliceCombine(obp.correctedSlice,obp.IndexMatrix);
         end
        %}
        
        function obp = CompletePreprocessing(obp)
            %performs all preproccesing algorithms
            %color alignment, stack alignment, image mergin
            
            %load image slices
%            obp = loadImageSlices(obp)
            %compute rgb alignment for all images of a single fish
            obp = rgbCorrection(obp);
            %apply rgb alignment on all images of a single fish
            obp = rgbWarp(obp);
            %compute stack alignment
            obp.sliceWarp = obp.sng_stackAlignment(obp.imageSliceCor, obp.scaleS, 'translation', obp.levelsS, obp.iterationsS);
            %apply stack alignment
            obp.correctedSlice = obp.sng_stackWarp(obp.imageSliceCor, obp.sliceWarp);
            %compute merge parameters
            [obp.IndexMatrix, obp.variance_sq] = obp.sng_StackDOF2(obp.correctedSlice, obp.variancedisksize);
            %apply image merge (extended dept of field)
            obp.mergedImage = obp.sng_sliceCombine(obp.correctedSlice, obp.IndexMatrix);
        end
        
        %%TODO add showfunctions for Indexmatrix,variance_sq
        function showExtendedDeptOfField(obp)
            sz = size(obp.correctedSlice{1});
            for k = 1:numel(obp.correctedSlice)
                figure; imagesc(uint8(obp.correctedSlice{k}))
                axis off tight equal
                set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
                text(sz(2)*0.95, sz(1)*0.05, num2str(k), 'FontSize', 20)
            end
            
%             for k = 1:numel(obp.correctedSlice)
%                 figure; imagesc(uint8(obp.variance_sq(:, :, k)));
%                 axis off tight equal
%                 set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
%                 text(sz(2)*0.95, sz(1)*0.05, num2str(k), 'FontSize', 20, 'Color', [1, 1, 1])
%             end
            
            SelectionImage = obp.correctedSlice;
            for k = 1:numel(obp.correctedSlice)
                SelectionImage{k}(repmat(obp.nSlices, 1, 1, 3) ~= k) = 0;
                figure; imagesc(uint8(SelectionImage{k}));
                axis off tight equal
                set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
                if k == 1;
                    text(sz(2)*0.95, sz(1)*0.05, num2str(k), 'FontSize', 20, 'Color', [0, 0, 0])
                else
                    text(sz(2)*0.95, sz(1)*0.05, num2str(k), 'FontSize', 20, 'Color', [1, 1, 1])
                end
            end
            
            figure; imagesc(obp.mergedImage)
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            axis off tight equal
        end
        function showRGBCorrectionZoomIn(obp)
            
            zoomn = obp.magnification;
            zoompoint2 = obp.zoomlocation;
            step = obp.nZoomInSteps;
            
            %RGB Correction figures and movies
            
            figure; imagesc(uint8(obp.imageSlice{obp.sliceNo}(:, :, 1:3))); axis equal off tight
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            
            writerObj1 = VideoWriter('zoom in initial slice', 'MPEG-4');
            writerObj1.FrameRate = 25; % to perform realtime movie
            open(writerObj1);
            figure; imagesc(uint8(obp.imageSlice{obp.sliceNo}(:, :, 1:3))); axis equal off tight
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            
            h = gca;
            siz = size((obp.imageSlice{1}(:, :, 1:3)));
            zoompoint1 = siz(1:2) / 2;
            [zoomp, yzp, xzp] = sng_smoothzoom(zoomn, zoompoint1, zoompoint2, step, h);
            
            for k = 1:step
                sng_zoom(zoomp(k), [yzp(k), xzp(k)], siz(1:2), h)
                frame = getframe(gcf);
                writeVideo(writerObj1, frame);
            end
            close(writerObj1);
            
            
        end
        function showRGBCorrectionZoomOut(obp)
            
            zoompoint1 = obp.zoomlocation - obp.rgbCropValues{obp.sliceNo}(1:2);
            zoomn = obp.magnification;
            step = obp.nZoomOutSteps;
            
            siz = size((obp.imageSliceCor{1}(:, :, 1:3)));
            zoompoint2 = siz(1:2) / 2;
            
            figure; imagesc(uint8(obp.imageSliceCor{obp.sliceNo})); axis equal off tight
            sng_zoom(zoomn(2), zoompoint1, siz(1:2), gca)
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            
            writerObj1 = VideoWriter('zoom out corrected slice', 'MPEG-4');
            writerObj1.FrameRate = 25; % to perform realtime movie
            open(writerObj1);
            figure; imagesc(uint8(obp.imageSliceCor{obp.sliceNo})); axis equal off tight
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            h = gca;
            [zoomp, yzp, xzp] = sng_smoothzoom(fliplr(zoomn), zoompoint1, zoompoint2, step, h);
            for k = 1:step
                sng_zoom(zoomp(k), [yzp(k), xzp(k)], siz(1:2), h);
                frame = getframe(gcf);
                writeVideo(writerObj1, frame);
            end
            close(writerObj1);
            
        end
        function showStackUncorrected(obp)
            
            %% Stack Correction figures and movies
            
            
            %frame parameters
            a = numel(obp.imageSlice); %scroll
            up = 1:a;
            down = linspace(a, 1, a);
            freeze = ones(1, numel(obp.imageSlice)); %freeze
            k = [up, a * freeze, down, freeze];
            k = repmat(k, 1, 4);
            
            writerObj1 = VideoWriter('stack scroll uncorrected', 'MPEG-4');
            writerObj1.FrameRate = 6; % to perform realtime movie
            open(writerObj1);
            figure; imagesc(uint8(obp.imageSliceCor{1})); axis equal off tight
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            sz = size(obp.imageSliceCor{1});
            
            for j = k
                imagesc(uint8(obp.imageSliceCor{j})); axis equal off tight;
                text(sz(2)*0.95, sz(1)*0.05, num2str(j), 'FontSize', 20);
                frame = getframe(gcf);
                writeVideo(writerObj1, frame);
            end
            close(writerObj1);
        end
        function showStackCorrected(obp)
            
            sz = size(obp.correctedSlice{1});
            %frame parameters
            a = numel(obp.imageSlice); %scroll
            up = 1:a;
            down = linspace(a, 1, a);
            freeze = ones(1, numel(obp.imageSlice)); %freeze
            k = [up, a * freeze, down, freeze];
            k = repmat(k, 1, 4);
            
            writerObj1 = VideoWriter('stack scroll corrected', 'MPEG-4');
            writerObj1.FrameRate = 6; % to perform realtime movie
            open(writerObj1);
            figure; imagesc(uint8(obp.correctedSlice{1})); axis equal off tight;
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            
            for j = k
                imagesc(uint8(obp.correctedSlice{j})); axis equal off tight
                text(sz(2)*0.95, sz(1)*0.05, num2str(j), 'FontSize', 20);
                frame = getframe(gcf);
                writeVideo(writerObj1, frame);
            end
            close(writerObj1);
            %}
            
            
            %{
 
             figure;imshow(filteredImage)
             figure;imshow(ImagesSliceCor{1})
             figure;imshow(Images2)
             figure;imshow(Icombined)
 
            %}
        end
    end
    
    methods(Static)
        function [ECCWarp, Img2] = sng_rgbAlignment(Img1, scale, transform, levels, iterations)
            %RGB2 does te same als RGB but the output functions are turned for speed
            %increase as you only want to now ECCWarp
            
            %log
            %20181022 renamed from sng_RGB2 to sng_rgbAlignment
            
            
            %set parameters
            %{
                      scale = 1/4;
                      par.transform = 'translation';
                      par.levels = 2;
                      par.iterations = 10; %iterations per level
                      Img1 = tempImage;
            %}
            
            %if ~exist(scale, 'var'); scale = 1 / 4; end
            %if ~exist(levels, 'var'); levels = 2; end
            %if ~exist(transform, 'var'); transform = 'translation'; end
            %if ~exist(iterations, 'var'); iterations = 10; end
            
            par.transform = transform;
            par.levels = levels;
            par.iterations = iterations; %iterations per level
            
            Img1sc = imresize(Img1, scale);
            
            ECCWarp{1} = iat_ecc(Img1sc(:, :, 2), Img1sc(:, :, 1), par) * (1 / scale);
            ECCWarp{2} = iat_ecc(Img1sc(:, :, 3), Img1sc(:, :, 1), par) * (1 / scale);
            % ECCWarp = iat_LucasKanade(Img1sc(:,:,j), Img1sc(:,:,1), par)*(1/scale);
            
            if nargout == 2
                [M, N, Oh] = size(Img1);
                [Img1(:, :, 2), supportECC{1}] = iat_inverse_warping(Img1(:, :, 2), ECCWarp{1}, 'translation', 1:N, 1:M);
                [Img1(:, :, 3), supportECC{2}] = iat_inverse_warping(Img1(:, :, 3), ECCWarp{2}, 'translation', 1:N, 1:M);
                support = supportECC{1} .* supportECC{2};
                %crop image, find the largest rectangle consist of ones
                A = numel(find(support(round(M/2), 1:end)));
                B = numel(find(support(1:end, round(N/2))));
                Img2 = reshape(Img1(repmat(logical(support), 1, 1, 3)), B, A, 3);
                %figure;imagesc(Img2);
            end
            
        end
        function [Img2, crop] = sng_rgbWarp(Img1, ECCWarp)
            %warp image according to warp info and remove not overlapping boundary
            
            %log
            %sng_RGB_IATwarp2: change cropping method and add output crop values
            %20181022 change sng_RGB_IATwarp2 to sng_rgbWarp
            
            %{
                         %Example
                         Img1 = ImageSlice{1}(:,:,1:3);
                         a = ECCWarp{1}{1}
                         b = ECCWarp{1}{2}
                         clear ECCWarp
                         ECCWarp{1} = a
                         ECCWarp{2} = b
            %}
            
            %preallocation
            [M, N, Oh] = size(Img1);
            
            [Img1(:, :, 2), supportECC{1}] = iat_inverse_warping(Img1(:, :, 2), ECCWarp{1}, 'translation', 1:N, 1:M);
            [Img1(:, :, 3), supportECC{2}] = iat_inverse_warping(Img1(:, :, 3), ECCWarp{2}, 'translation', 1:N, 1:M);
            support = supportECC{1} .* supportECC{2};
            
            %{
                        figure;imagesc(Img1)
                        figure;imagesc(support)
            %}
            
            %old crop method
            %{
                        %crop image, find the largest rectangle consist of ones
                        A = numel(find(support(round(M/2),1:end)));
                        B = numel(find(support(1:end,round(N/2))));
                        Img2 = reshape(Img1(repmat(logical(support),1,1,3)),B,A,3);
                        %figure;imshow(uint8(Img2{j}));
            %}
            
            xr1 = find(support(round(M/2), 1:end), 1, 'first');
            xr2 = find(support(round(M/2), 1:end), 1, 'last');
            yr1 = find(support(1:end, round(N/2)), 1, 'first');
            yr2 = find(support(1:end, round(N/2)), 1, 'last');
            
            Img2 = Img1(yr1:yr2, xr1:xr2, :);
            crop = [yr1, yr2; xr1, xr2];
        end
        function [Img2, Imgband] = sng_bandpassFilter(Img, sigma1, sigma2)
            %bandpassfilter RGB image to focus on spots
            %crops image to remove boundary and scale bar
            %calculates center of mass according bandpass filtered image
            
            %log
            %20181022 change name sng_RGBselection to sng_bandpassFilter
            
            %{
                       Img = ImageSlice{k}(:,:,1:3);
            %}
            %Img = imread(Imgloc1{1}{1});
            %Img = Img(:,:,1:3);
            
            %if ~exist('sigma1','var');sigma1 = 1;end
            %if ~exist('sigma2','var');sigma2 = 4;end
            
            
            %gaussian blur and difference image
            Img = double(Img(:, :, 1:3));
            %Imglow1 = imgaussfilt(Img,sigma1);
            %Imglow2 = imgaussfilt(Img,sigma2);
            %for older matlab versions:
            h1 = fspecial('gaussian', 2*ceil(2*sigma1)+1, sigma1);
            Imglow1 = imfilter(Img, h1, 'replicate');
            h2 = fspecial('gaussian', 2*ceil(2*sigma2)+1, sigma2);
            Imglow2 = imfilter(Img, h2, 'replicate');
            Imgband = (Imglow1 - Imglow2);
            
            %geen balkje en geen randen
            Img1 = Imgband(20:round(size(Img, 1)*(11 / 12)), 20:size(Img, 2)-20, 1:3);
            Img2a = Img1;
            Img2b = Img1;
            
            Img2a(Img1 <= 0) = 0;
            Img2b(Img1 >= 0) = 0;
            
            Img2 = cat(1, Img2a, -Img2b);
            
            
            %figure;imagesc(uint8(Img2*10));
            
            
            % %select area around center of mass detail image
            % %does not work for some fishes
            % [comx,comy] = sng_CenterOfMassColor(Img1,1);
            % x1 = round(comx-300);if x1 < 1; x1 = 1;end;
            % x2 = round(comx+300);if x2 > size(Img1,2); x2 = size(Img1,2);end;
            % y1 = round(comy-300);if y1 < 1; y1 = 1;end;
            % y2 = round(comy+300);if y2 > size(Img1,1); y2 = size(Img1,1);end;
            % Img2 = Img1(y1:y2,x1:x2,:);
            
            %{
                        figure;imagesc(uint8(Img(:,:,1:3)));
 
 
                        figure;imagesc(0.15*abs(Img1(:,:,1:3)))
                        figure;imagesc(0.15*Img2a(:,:,1:3))
                        figure;imagesc(-0.15*Img2b(:,:,1:3))
 
                        l = Img1-min(Img1(:));
 
                        figure;imagesc(QImg);
                        figure;imagesc(Img2(:,:,1));colormap('gray')
                        figure;imagesc(Img2(:,:,2));colormap('gray')
                        figure;imagesc(Img2(:,:,3));colormap('gray')
                        sng_figureslide
            %}
            
            % [Img2a,ECCWarp] = sng_RGB(Img1b,1/2,'translation',2,20);
            % %[Img2b,tform] = sng_RGB_imreg(Img1b,1,'translation');
            % figure;imagesc(10*Img2a)
            % figure;imagesc(10*Img2b)
            %
            % QImg = Img(:,:,1:3);
            % for j=2:3
            % QImg(:,:,j) = imwarp(QImg(:,:,j),tform{j},'OutputView',imref2d(size(QImg)));
            % end
            %
            
        end
        function [ECCWarp, cellImg1] = sng_stackAlignment(cellImg1, scale, transform, levels, iterations, varargin)
            %allign all fishes with different dept of field from location
            %the first fish acts as a template
            %imput is at least 2 images in a cell structure
            %iatool.net/ for more info
            %by the translation method, information has been lost, consider to use only
            %the translation values
            %examples
            %sng_AllignFishBatch3(cellImg1,scale,transform,levels,iterations,Fullpath)
            
            %log
            %20181022 change name sng_AllignFishBatch3 to sng_stackAlignment
            
            % scale = 1/4;
            % transform = 'affine'; %'translation';
            % levels = 3;
            % iterations = 20; %iterations per level
            % cellImg1 = Meanfish;
            
            %algorithm parameters
            par.transform = transform;
            par.levels = levels;
            par.iterations = iterations; %iterations per level
            
            Img1sc = imresize(cellImg1{1}, scale); %resize to speed up
            Img1sc = Img1sc(1:round(size(Img1sc, 1)*(11 / 12)), :, 1:3); %to exclude scalebar
            
            [M, N, Oh] = size(cellImg1{1});
            support = true(M, N);
            ECCWarp = cell(1, numel(cellImg1));
            
            for j = 2:numel(cellImg1)
                Img2sc = imresize(cellImg1{j}, scale);
                %to remove the scale bar
                Img2sc = Img2sc(1:round(size(Img2sc, 1)*(11 / 12)), :, 1:3);
                ECCWarp{j} = iat_ecc(Img2sc, Img1sc, par) * (1 / scale);
            end
            
            if nargin >= 6 || nargout >= 2
                for j = 2:numel(cellImg1)
                    [cellImg1{j}, supportECC] = iat_inverse_warping(cellImg1{j}, ECCWarp{j}, par.transform, 1:N, 1:M);
                    cellImg1{j} = uint8(cellImg1{j});
                    support = support & logical(supportECC);
                end
                %crop image, find the largest rectangle consist of ones
                A = numel(find(support(round(M/2), 1:end)));
                B = numel(find(support(1:end, round(N/2))));
                
                for j = 1:numel(cellImg1)
                    cellImg1{j} = reshape(cellImg1{j}(repmat(support, 1, 1, Oh)), B, A, Oh);
                    if nargin >= 6
                        Fullpath = varargin{1};
                        imwrite(uint8(cellImg1{j}), Fullpath, ...
                            'WriteMode', 'append', 'Compression', 'none');
                    end
                    %figure;imshow(uint8(cellImg1{j}));
                end
            end
            
            
            % for k=1:numel(cellImg1)
            %     figure;imagesc(uint8(Meanfish{k}))
            %     figure;imagesc(uint8(cellImg1{k}))
            %     figure;imagesc(support);
            % end
            % sng_figureslide
            
            
        end
        function [cellImg1] = sng_stackWarp(cellImg1, ECCWarp)
            
            par.transform = 'translation';
            [M, N, Oh] = size(cellImg1{1});
            support = true(M, N);
            
            for j = 2:numel(cellImg1)
                [cellImg1{j}, supportECC] = iat_inverse_warping(cellImg1{j}, ECCWarp{j}, par.transform, 1:N, 1:M);
                cellImg1{j} = uint8(cellImg1{j});
                support = support & logical(supportECC);
            end
            %crop image, find the largest rectangle consist of ones
            A = numel(find(support(round(M/2), 1:end)));
            B = numel(find(support(1:end, round(N/2))));
            
            for j = 1:numel(cellImg1)
                cellImg1{j} = reshape(cellImg1{j}(repmat(support, 1, 1, Oh)), B, A, Oh);
                %figure;imshow(uint8(cellImg1{j}));
            end
        end
        function [IndexMatrix, variance_sq, Icombined] = sng_StackDOF2(ImgCell, win_size, filt)
            %Version sng_StackDOF2
            %   change the slicen ouput to a logic variable IndexMatrix which costs
            %   less memory and can be applied faster to acchieve the Icombined image
            
            %{
   Img = CorrectedSlice;
   win_size = variancedisksize
            %}
            
            if ~exist('win_size', 'var')
                win_size = 3;
            end
            if ~exist('filt', 'var')
                filt = 'max';
            end
            
            n = numel(ImgCell);
            
            [M, N, ~] = size(ImgCell{1});
            variance_sq = zeros(M, N, size(ImgCell, 1));
            
            
            % Making spherical structuring element for determing the local variance
            % with radisu win_size
            struc_el = fspecial('disk', win_size) > 0;
            
            % Calculating the local variance for each pixel and each slice
            for i = 1:n
                variance = stdfilt(ImgCell{i}, struc_el);
                variance_sq(:, :, i) = sum(variance, 3);
                
                %variance_sq(:,:,i) = sqrt(sum(variance.^2,3));
                
            end
            
            % Making an MxN image with the slice indexes of the maximum/minimum local variance values
            if strcmp(filt, 'min')
                [~, slicen] = min(variance_sq, [], 3);
            elseif strcmp(filt, 'max')
                [~, slicen] = max(variance_sq, [], 3);
            else
                error('wronginput')
            end
            
            %create a logic image which selects pixels form the different slices
            IndexMatrix = false(size(variance_sq));
            Btemp = (slicen(:) - 1) * M * N;
            Ctemp = 1:M * N;
            IndexMatrix(Btemp+Ctemp') = true;
            
            if nargout >= 3
                Icombined = sng_SliceCombine(ImgCell, IndexMatrix);
            end
            
            
            %{
   figure;imagesc(Img{1})
   figure;imagesc(Img{2})
   figure;imagesc(Img{3})
   figure;imagesc(Img{4})
 
 
   figure;imagesc(variance_sq(:,:,1))
   figure;imagesc(variance_sq(:,:,2))
   figure;imagesc(variance_sq(:,:,3))
   figure;imagesc(variance_sq(:,:,4))
 
   figure;imagesc(Icombined)
 
            %}
            %stacked double used
            %value win_size is char, should be int
            
            
        end
        function Icombined = sng_sliceCombine(ImgCell, IndexMatrix)
            %creates combined image based on Indexmatrix called in sng_StackDof
            
            Icombined = zeros(size(ImgCell{1}), 'uint8');
            for k = 1:numel(ImgCell)
                ImgCell{k}(~repmat(IndexMatrix(:, :, k), 1, 1, 3)) = 0;
                Icombined = Icombined + ImgCell{k};
            end
            
        end
    end
    
end