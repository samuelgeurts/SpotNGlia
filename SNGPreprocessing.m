classdef SNGPreprocessing
    
    properties %(Access = private)
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
        
        
        nSlices = [];
        
        % RgbCorrection
        
        ColorWarp; % = ECCWarp;
        NCC_Img_before;
        NCC_Img_after;
        NCC_BP_before;
        NCC_BP_after;
        rgbCropValues;
        % step 2
        SliceWarp; % = ECCWarp2;
        
        
    end
    
    
    properties(Transient = true)
        %IMAGES
        ImageSlice = []
        ImageSliceCor = []
        CorrectedSlice = []
        %RgbCorrection
        FiltIm
    end
    
    
    methods
        function obp = SNGPreprocessing(SpotNGliaObject) %constructor
            %obr = obr.zfproperties(SpotNGliaObject.ZFParameters);
            
            %set subfields of the registration field of the SngInputParameter
            %object as properties the obr i.e. SNGAllignment-object
            inputFields = fields(SpotNGliaObject.SngInputParameters.Preprocessing);
            for iInputField = 1:numel(inputFields)
                obp.(inputFields{iInputField}) = SpotNGliaObject.SngInputParameters.Preprocessing.(inputFields{iInputField});
            end
        end
        
        function obp = RgbCorrection(obp, ImageSlice)
            %applys rgb correction i.e. translate channels and slices using iat
            
            if exist('ImageSlice', 'var')
                obp.ImageSlice = ImageSlice;
            else
                ImageSlice = obp.ImageSlice;
            end
            
            
            obp.nSlices = numel(obp.ImageSlice);
            
            %preallocation
            NCC_Img_before = cell(1, obp.nSlices);
            NCC_Img_after = cell(1, obp.nSlices);
            NCC_BP_before = cell(1, obp.nSlices);
            NCC_BP_after = cell(1, obp.nSlices);
            ECCWarp = cell(1, obp.nSlices);
            ImageSliceCor = cell(1, obp.nSlices);
            crop = cell(1, obp.nSlices);
            
            obp.FiltIm = cell(1, obp.nSlices);
            
            %% does RGB correction on all images of one fish
            for iSlice = 1:obp.nSlices
                
                if obp.onoff
                    %bandpass filter focussed on spotsize passing
                    FilteredImage = sng_RGBselection(ImageSlice{iSlice}(:, :, 1:3), obp.sigmalp, obp.sigmahp); %%TODO put this function in this class
                else
                    %no bandpass filtering
                    FilteredImage = obp.ImageSlice{iSlice}(:, :, 1:3);
                end
                
                obp.FiltIm{iSlice} = FilteredImage;
                
                
                %{
                   figure;imagesc((FilteredImage(:,:,1:3)));axis equal off tight
                   figure;imagesc(uint8(ImageSlice{iSlice}(:,:,1:3)));axis equal off tight
                %}
                
                %calculates the needed warp
                [ECCWarp{iSlice}, ~] = sng_RGB2(FilteredImage(:, :, 1:3), obp.scaleC, 'translation', obp.levelsC, obp.iterationsC);
                %warp images
                [ImageSliceCor{iSlice}, crop{iSlice}] = sng_RGB_IATwarp2(ImageSlice{iSlice}(:, :, 1:3), ECCWarp{iSlice});
                %imwrite(uint8(Img2), fullFileName{l}, 'WriteMode', 'append', 'Compression','none');
                
                %normal correlation
                %{
                   CC1 = corr2(FilteredImage(:,:,1),FilteredImage(:,:,2));
                   CC2 = corr2(FilteredImage2(:,:,1),FilteredImage2(:,:,2));
                   CC3 = corr2(ImageSlice{iSlice}(:,:,1),ImageSlice{iSlice}(:,:,3));
                   CC4 = corr2(ImageSliceCor{iSlice}(:,:,1),ImageSliceCor{iSlice}(:,:,3));
                %}
                
                [CorrectedFilteredImage, ~] = sng_RGB_IATwarp2(FilteredImage, ECCWarp{iSlice});
                
                %{
                     figure;imagesc((ImageSliceCor{iSlice}(:,:,1:3)));axis equal off tight
                     figure;imagesc((CorrectedFilteredImage(:,:,1:3)));axis equal off tight
                %}
                
                %original image before and after normalised correlation (green-red and blue-red channel)
                NCC_Img_before{iSlice}(1) = sng_NCC(ImageSlice{iSlice}(:, :, 1), ImageSlice{iSlice}(:, :, 2));
                NCC_Img_before{iSlice}(2) = sng_NCC(ImageSlice{iSlice}(:, :, 1), ImageSlice{iSlice}(:, :, 3));
                NCC_Img_after{iSlice}(1) = sng_NCC(ImageSliceCor{iSlice}(:, :, 1), ImageSliceCor{iSlice}(:, :, 2));
                NCC_Img_after{iSlice}(2) = sng_NCC(ImageSliceCor{iSlice}(:, :, 1), ImageSliceCor{iSlice}(:, :, 3));
                %bandpass filtered before and after normalised correlation  (green-red and blue-red channel)
                NCC_BP_before{iSlice}(1) = sng_NCC(FilteredImage(:, :, 1), FilteredImage(:, :, 2));
                NCC_BP_before{iSlice}(2) = sng_NCC(FilteredImage(:, :, 1), FilteredImage(:, :, 3));
                NCC_BP_after{iSlice}(1) = sng_NCC(CorrectedFilteredImage(:, :, 1), CorrectedFilteredImage(:, :, 2));
                NCC_BP_after{iSlice}(2) = sng_NCC(CorrectedFilteredImage(:, :, 1), CorrectedFilteredImage(:, :, 3));
                
            end
            
            obp.rgbCropValues = crop;
            obp.ColorWarp = ECCWarp;
            obp.NCC_Img_before = NCC_Img_before;
            obp.NCC_Img_after = NCC_Img_after;
            obp.NCC_BP_before = NCC_BP_before;
            obp.NCC_BP_after = NCC_BP_after;
            obp.ImageSliceCor = ImageSliceCor;
            
            
        end
        
        function obp = StackCorrection(obp)
            [ECCWarp2, CorrectedSlice] = sng_AllignFishBatch3(obp.ImageSliceCor, obp.scaleS, 'translation', obp.levelsS, obp.iterationsS);
            
            obp.SliceWarp = ECCWarp2;
            obp.CorrectedSlice = CorrectedSlice;
        end
        
        
        function ShowRGBCorrectionZoomIn(obp)
            
            zoomn = obp.magnification;
            zoompoint2 = obp.zoomlocation;
            step = obp.nZoomInSteps;
            
            %% RGB Correction figures and movies
            
            figure; imagesc(uint8(obp.ImageSlice{obp.sliceNo}(:, :, 1:3))); axis equal off tight
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
                       
            writerObj1 = VideoWriter('zoom in initial slice', 'MPEG-4');
            writerObj1.FrameRate = 25; % to perform realtime movie
            open(writerObj1);
            figure; imagesc(uint8(obp.ImageSlice{obp.sliceNo}(:, :, 1:3))); axis equal off tight
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            
            h = gca;
            siz = size((obp.ImageSlice{1}(:, :, 1:3)));
            zoompoint1 = siz(1:2) / 2;
            [zoomp, yzp, xzp] = sng_smoothzoom(zoomn, zoompoint1, zoompoint2, step, h);
            
            for k = 1:step
                sng_zoom(zoomp(k), [yzp(k), xzp(k)], siz(1:2), h)
                frame = getframe(gcf);
                writeVideo(writerObj1, frame);
            end
            close(writerObj1);
            
            
        end
        
        function ShowRGBCorrectionZoomOut(obp)
            
            zoompoint1 = obp.zoomlocation - obp.rgbCropValues{obp.sliceNo}(1:2);
            zoomn = obp.magnification;
            step = obp.nZoomOutSteps;
                      
            siz = size((obp.ImageSliceCor{1}(:, :, 1:3)));
            zoompoint2 = siz(1:2) / 2;
            
            figure; imagesc(uint8(obp.ImageSliceCor{obp.sliceNo})); axis equal off tight
            sng_zoom(zoomn(1), zoompoint1, siz(1:2), gca)
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
                    
            writerObj1 = VideoWriter('zoom out corrected slice', 'MPEG-4');
            writerObj1.FrameRate = 25; % to perform realtime movie
            open(writerObj1);
            figure; imagesc(uint8(obp.ImageSliceCor{obp.sliceNo})); axis equal off tight
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            h = gca;
            [zoomp, yzp, xzp] = sng_smoothzoom(zoomn, zoompoint1, zoompoint2, step, h);
            for k = 1:step
                sng_zoom(zoomp(k), [yzp(k), xzp(k)], siz(1:2), h);
                frame = getframe(gcf);
                writeVideo(writerObj1, frame);
            end
            close(writerObj1);
             
        end
        
        function ShowStackUncorrected(obp)
            
            %% Stack Correction figures and movies
            

            
            
            %frame parameters
            a = numel(obp.ImageSlice); %scroll
            up = 1:a;
            down = linspace(a, 1, a);
            freeze = ones(1, numel(obp.ImageSlice)); %freeze
            k = [up, a * freeze, down, freeze];
            k = repmat(k, 1, 4);
            
            writerObj1 = VideoWriter('stack scroll uncorrected', 'MPEG-4');
            writerObj1.FrameRate = 6; % to perform realtime movie
            open(writerObj1);
            figure; imagesc(uint8(obp.ImageSliceCor{1})); axis equal off tight
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            sz = size(obp.ImageSliceCor{1});
            
            for j = k
                imagesc(uint8(obp.ImageSliceCor{j})); axis equal off tight;
                t = text(sz(2)*0.95, sz(1)*0.05, num2str(j), 'FontSize', 20);
                frame = getframe(gcf);
                writeVideo(writerObj1, frame);
            end
            close(writerObj1);
        end
        
        function ShowStackCorrected(obp)
            
            
            %frame parameters
            a = numel(obp.ImageSlice); %scroll
            up = 1:a;
            down = linspace(a, 1, a);
            freeze = ones(1, numel(obp.ImageSlice)); %freeze
            k = [up, a * freeze, down, freeze];
            k = repmat(k, 1, 4);
            
            writerObj1 = VideoWriter('stack scroll corrected', 'MPEG-4');
            writerObj1.FrameRate = 6; % to perform realtime movie
            open(writerObj1);
            figure; imagesc(uint8(obp.CorrectedSlice{1})); axis equal off tight;
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            
            for j = k
                imagesc(uint8(obp.CorrectedSlice{j})); axis equal off tight
                t = text(sz(2)*0.95, sz(1)*0.05, num2str(j), 'FontSize', 20);
                frame = getframe(gcf);
                writeVideo(writerObj1, frame);
            end
            close(writerObj1);
            %}
            
            
            %{
 
   figure;imshow(FilteredImage)
   figure;imshow(ImagesSliceCor{1})
   figure;imshow(Images2)
   figure;imshow(Icombined)
 
            %}
        end
        
    end
end