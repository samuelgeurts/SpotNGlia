classdef SNGAlignment
    % this function is a class which makes use of properties from the class SpotNGlia but
    % is also called form the class SpotNGlia. I dont know if this is a nice solution but 
    % matlab dont know subclasses.
    % SNGalignment comes form the older function AlignmentSNG which is an ordinary function.
    
%{
    To do adjustment on the this class, use obr = SNGAlignment(obj)
    fn = 42
    CombinedFish = imread([obj.SavePath, '/', 'CombinedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
    obr.Icombined = CombinedFish;
%}  
    
    properties %(Access = private)
        %INPUTPARAMETERS
        %from CompleteTemplate
        template = []
        
        %BackgroundRemoval
        Method = []
        Smooth = []
        ChannelMethod = []
        %RotationalAlignment
        AngleSteps1 = []
        Scale1 = []
        AngleRange2 = []
        AngleSteps2 = []
        FinalCenteredYCrop = []
        %ScaleFlipAlignment
        ScaleRange = []
        ScaleSteps = []
        Scale2 = []
        %FinalTemplateMatching
        RemoveTail = []
        Scale3 = [] %lowered initial scale to increase computation speed for correlation
        AffineMethod = [] %affine or translation
        Levels = [] %number of scaled levels correlation is performed
        Iterations = [] %number of iterations iat is performed
    end
    properties
        %OUTPUTPARAMETERS
        
        %background removal
        BackgroundMask
        BackgroundThreshold
        stretchlim_1procent
        stretchlim_0procent
        
        %Rotation correlation output
        xcom2
        ycom2
        theta1 %rotation angle1
        maxangle1 %raw angle1
        peaksep1 %quality 1-peak2/peak1
        ncorr1 %correlation1
        angles1
        theta3 %rotation angle2'
        maxangle2 %raw angle2
        ncorr2 %correlation2
        angles2
        %
        ref25
        ref30
        tform1 %rotation correlation
        tform2 %scale and flip correlation
        
        %Scale correlation output        
        scalebest; %max correlation value
        Tflip; 
        Tscale;
        Ttranslate;
        Corrbest; %Best scale correlation image
        maxk; %max correlation plot
        csflip;
        
        %   
        tform3 %template matching tailcut
        tform4 %template matching translation
        CorCoef %Template correlation coefficient
        tform_1234 %complete transformation matrix
        MeanFishColor
        MaxFishColor
    end
    properties(Transient = true)
        %IMAGES
        Icombined
        Ialigned
        
        %BackgroundRemoval
        I20
        %RotationCorrelation
        I21
        I22
        I22_1
        I22_2
        I22_3
        I22_4
        I23
        I24
        I25
        I30
        %ScaleAlignment
        I31
        I32
        I40
        %Final Correlation with IAT
        I55
    end
    properties
        RotationIndices = [13, 32, 51, 70, 89]; %between 1 and 100
    end
    
    methods
        function obr = SNGAlignment(SpotNGliaObject) %constructor
            %obr = obr.zfproperties(SpotNGliaObject.ZFParameters);
            
            %set subfields of the registration field of the SngInputParameter 
            %object as properties the obr i.e. SNGAllignment-object
            inputFields = fields(SpotNGliaObject.SngInputParameters.Registration);
            for iInputField = 1:numel(inputFields)
                obr.(inputFields{iInputField}) = SpotNGliaObject.SngInputParameters.Registration.(inputFields{iInputField});
            end
            
%             if isempty(SpotNGliaObject.CompleteTemplate)    
%                 SpotNGliaObject = SpotNGliaObject.LoadTemplate;
%             end
            
            obr.template = SpotNGliaObject.CompleteTemplate.Template;
        end
        function obr = BackgroundRemoval(obr, Icombined)
            
            if exist('Icombined', 'var')
                obr.Icombined = Icombined;
            end
            
            
            %{
                   [xcom1,ycom1] = sng_CenterOfMassColor(imcomplement(I10),1);
                   shift = (fliplr(size(I10(:,:,1))/2)-[xcom1,ycom1]);
                   I11= imtranslate(I10,shift,'OutputView','same');
                   figure;imagesc(I11)
            %}
            
            %horizontal allignment
            %remove background to compute center of mass in horizontalfish well
            
            [obr.I20, ...
                obr.BackgroundMask, ...
                obr.BackgroundThreshold] ...
                = sng_RemoveBackgroundColor3( ...
                obr.Icombined, ...
                obr.Method, ...
                obr.Smooth, ...
                obr.ChannelMethod);
            
            %stretchlim parameters (not used for registration so far)
            for k1 = 1:3
                temp = obr.I20(:, :, k1);
                obr.stretchlim_1procent(1:2, k1) = stretchlim(temp(~obr.BackgroundMask), 0.01);
                obr.stretchlim_0procent(1:2, k1) = stretchlim(temp(~obr.BackgroundMask), 0);
            end
            
        end
        function obr = Rotationalalignment(obr, I20)
            
            if exist('I20', 'var')
                obr.I20 = I20;
            end
            
            obr.angles1 = linspace(2*pi/obr.AngleSteps1, 2*pi, obr.AngleSteps1);
            
            
            % % first iteration global angle calculation
            obr.I21 = imcomplement(uint8(obr.I20));
            %calculate center of mass.
            [xcom1, ycom1] = sng_CenterOfMassColor(obr.I21, 1);
            %scale image for global calculation
            [obr.I22, obr.xcom2, obr.ycom2] = sng_scale(obr.I21, xcom1, ycom1, obr.Scale1);
            %apply normalized rotation correlation
            %[obr.theta1, obr.maxangle1, obr.peaksep1, obr.ncorr1] = sng_NormalizedRotationCorrelation4(obr.I22, obr.angles1, [xcom2, ycom2]);
            [obr.theta1, obr.maxangle1, obr.peaksep1, obr.ncorr1, obr.I22_1, obr.I22_2, obr.I22_3, obr.I22_4] = ...
                SNGAlignment.NormalizedRotationCorrelation(obr.I22, obr.angles1, [obr.xcom2, obr.ycom2], obr.RotationIndices);
            
            
            %{
      figure;imagesc(obr.I20)
      figure;imagesc(obr.I21)
      figure;imagesc(obr.I22)
                %}
                
                %rotate
                Ttrans = [1, 0, 0; 0, 1, 0; -xcom1, -ycom1, 1];
                Trotate = [cos(obr.theta1), - sin(obr.theta1), 0; sin(obr.theta1), cos(obr.theta1), 0; 0, 0, 1];
                tform_rc1 = affine2d(Ttrans*Trotate);
                
                % % a more precise rotation calculation
                if ~isempty(obr.AngleRange2)
                    
                    obr.angles2 = pi + 2 * pi * linspace(-obr.AngleRange2, obr.AngleRange2, obr.AngleSteps2);
                        
                    [obr.I23, ref23] = imwarp(uint8(obr.I21), tform_rc1, 'FillValues', 0);
                    %figure;imagesc(Img3)
                    %crops image s.t. iamges become more symmetric
                    [obr.I24, ~] = sng_fishcrop3(obr.I23, ref23, [size(obr.I23, 2), size(obr.I23, 1) / 2]); %kan meer exact
                    [xcom4, ycom4] = sng_CenterOfMassColor(obr.I24, 1);
                    
                    %more precise angle calculation
                    [theta2, obr.maxangle2, ~, obr.ncorr2] = SNGAlignment.NormalizedRotationCorrelation(obr.I24, obr.angles2, [xcom4, ycom4]);
                    obr.theta3 = obr.theta1 + theta2 + pi;
                    
                    %final rotation
                    Trotate = [cos(obr.theta3), - sin(obr.theta3), 0; sin(obr.theta3), cos(obr.theta3), 0; 0, 0, 1];
                    tform_rc1 = affine2d(Ttrans*Trotate);
                end
                
                
                %{
      figure;imagesc(obr.I23);sng_imfix
      figure;imagesc(obr.I24);sng_imfix
      figure;plot(obr.ncorr1)
      figure;plot(obr.ncorr2)
                %}
                
                % % transform such that the center lays on the y axis
                %transformation with (0.0) as the center coordinate
                [obr.I25, obr.ref25] = imwarp(obr.I20, tform_rc1, 'FillValues', 255);
                %reference object with cropped Y-axis such that the center is the real center, x-reference stays the same
                ref26 = imref2d([obr.FinalCenteredYCrop, size(obr.I25, 2)], obr.ref25.XWorldLimits, [0.5 - obr.FinalCenteredYCrop / 2, 0.5 + obr.FinalCenteredYCrop / 2]);
                %cropped image
                %I26 = imwarp(obr.I20,tform_rc1,'FillValues',255,'OutputView',ref26);
                
                tform1ToPixel = affine2d([1, 0, 0; 0, 1, 0; -ref26.XWorldLimits(1), - ref26.YWorldLimits(1), 1]);
                
                tform_rc2 = affine2d(tform_rc1.T*tform1ToPixel.T);
                
                obr.tform1 = tform_rc2;
                [obr.I30, obr.ref30] = imwarp(obr.I20, obr.tform1, 'FillValues', 255, 'OutputView', imref2d([obr.FinalCenteredYCrop, size(obr.I25, 2)]));
                
                %{
      figure;imagesc(ref25.XWorldLimits,ref25.YWorldLimits,uint8(obr.I25));sng_imfix
      figure;imagesc(ref26.XWorldLimits,ref26.YWorldLimits,uint8(I26));
      figure;imagesc(obr.ref30.XWorldLimits,obr.ref30.YWorldLimits,uint8(obr.I30));sng_imfix
                %}
        end
        function obr = ShowRotationMovie(obr)
            
            % show fasthorizontal correlation with movie
            
            [obr.theta1, obr.maxangle1, obr.peaksep1, obr.ncorr1, obr.I22_1, obr.I22_2, obr.I22_3] = ...
                SNGAlignment.NormalizedRotationCorrelation(obr.I22, obr.angles1, [obr.xcom2, obr.ycom2], 1:100);
            
            figure; imagesc(obr.I20); axis equal off tight
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            figure; imagesc(obr.I21); axis equal off tight
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            figure; imagesc(obr.I22); axis equal off tight
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            figure; imagesc(obr.I22_1); axis off tight equal;
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            figure; imagesc(obr.I22_2); axis off tight equal;
            set(gca, 'position', [0, 0, 1, 1], 'units', 'normalized')
            
            figure;
            subplot = @(m, n, p) subtightplot(m, n, p, [0.01, 0.05], [0.1, 0.01], [0.1, 0.01]);
            set(gcf, 'color', 'white')
            
            writerObj1 = VideoWriter('rotation correlation', 'MPEG-4');
            writerObj1.FrameRate = 10; % to perform realtime movie
            open(writerObj1);
            
            subplot(2, 2, 1); imagesc(obr.I22_1); axis off tight equal;
            subplot(2, 2, 2); imagesc(obr.I22_2); axis off tight equal;
            for k = 1:100
                subplot(2, 2, 2); imagesc(obr.I22_3{k}); axis off tight equal
                subplot(2, 1, 2); plot(obr.angles1(1:k), obr.ncorr1(1:k));
                %xlim([0,360]);
                ylim([0, 1])
                xlim([0, 2 * pi])
                xlabel('angle')
                frame = getframe(gcf);
                writeVideo(writerObj1, frame);
            end
            close(writerObj1);
            
        end
        function obr = Scalealignment(obr, I30)
            
            if exist('I30', 'var')
                obr.I30 = I30;
            end
                         
%                function [tform_2,ScaleCorrelationOutput,I40] = sng_ScaleFish4(I40,template,ScaleRange_rg,ScaleSteps_rg,Scale2_rg)
%                function [obr.tform,obr.ScaleCorrelationOutput,obr.I40] = sng_ScaleFish4(obr.30,obr.template,obr.ScaleRange,obr.ScaleSteps,obr.Scale2)           

                    scales = linspace(obr.ScaleRange(1),obr.ScaleRange(2),obr.ScaleSteps);

                    I30Gray = rgb2gray(uint8(obr.I30));
                    templateGray = rgb2gray(uint8(obr.template));

                    [templateScaled,obr.I31] = sng_scale4(obr.Scale2,templateGray,I30Gray);

                    %minimal overlap number of pixels for cross correlation for speeding up
                    mins = min(size(templateScaled),size(obr.I31));
                    C = mins(1) * mins(2)/6;

                    %preallocation
                    obr.maxk = zeros(2,numel(scales));
                    ck = zeros(2,numel(scales));
                    siz = zeros(2,numel(scales));

                    %when only the calculation scale as output
                    for j=1:numel(scales)  
                        %templatevar = sng_scale2(template,scales(j));
                        I31Scaled = sng_scale4(scales(j),obr.I31);
                        ICorCoef = normxcorr2_general(I31Scaled,templateScaled,C);
                        ICorCoefFlip = normxcorr2_general(rot90(I31Scaled,2),templateScaled,C);
                        [obr.maxk(1,j),ck(1,j)] = max(ICorCoef(:)); %maxk is the max correlationcoef per scale
                        [obr.maxk(2,j),ck(2,j)] = max(ICorCoefFlip(:));
                        siz(:,j) = size(ICorCoef);
                    end

                    %find de maxima of all correlation image maxima
                    %csflip is 1 for no flip and 2 for a flipped template
                    %csj could be 1 to numel(scales) dependent on the best match
                    [~, cs] = max(obr.maxk(:));                %cs is de index with the highest correlation
                    [obr.csflip,csj] = ind2sub(size(obr.maxk),cs); %if csflip=1 than no flip has the highest correlation, if csflip=2 than the flipped image has the highest correlation
                    obr.scalebest = scales(csj);


                    %{
                    figure;plot(1:200,obr.maxk(1,:))
                    hold on; plot(1:200,obr.maxk(2,:))

                    figure;imagesc(I31Scaled)
                    figure;imagesc(ICorCoef)
                    figure;imagesc(ICorCoefFlip)    
                        I31Scaled = sng_scale4(obr.scalebest,obr.I31);
                        if obr.csflip == 1; I101 = normxcorr2_general(I31Scaled,templateScaled,C);end;
                        if obr.csflip == 2; I101 = normxcorr2_general(rot90(I31Scaled,2),templateScaled,C);end;
                        figure;imagesc(I101)
                    %}

                    %recalculates the best template match and correlation image
                        %I40(I40(:,:,1)==255) = 100; <skipped this line
                        %I40 is van obr.I30 dont know the sence
                        [sy,sx,~] = size(obr.I30);

                        I31Scaled = sng_scale4(obr.scalebest,obr.I31);
                        if obr.csflip == 1           
                            obr.Tflip = [1 0 0;0 1 0; 0 0 1];
                            obr.Corrbest = normxcorr2_general(I31Scaled,templateScaled,C);
                        elseif obr.csflip == 2 
                            obr.Tflip = [-1 0 0;0 -1 0; sx,sy 1];        
                            obr.Corrbest = normxcorr2_general(rot90(I31Scaled,2),templateScaled,C);
                        end

                    %peak coordinates in correlationimage
                    [y,x] = ind2sub(siz(:,csj),ck(obr.csflip,csj));
                    offset = [y,x] - size(I31Scaled);
                    scaledoffset = offset / obr.Scale2;

                    obr.Ttranslate = [1,0,0;0,1,0;scaledoffset(2),scaledoffset(1),1]; %transform back to pixel coordinates
                    obr.Tscale = [scales(csj) 0 0; 0 scales(csj) 0;0 0 1];    
                    obr.tform2 = affine2d(obr.Tflip*obr.Tscale*obr.Ttranslate);

                    obr.I40 = imwarp(obr.I30,obr.tform2,'FillValues',255,'OutputView',imref2d(size(obr.template(:,:,1))));

                    %{
                        %to test performance indiviual transforms
                        I41 = imwarp(I40,affine2d(obr.Tscale));figure;imagesc(uint8(I41))
                        I42 = imwarp(I41,affine2d(obr.Tflip));figure;imagesc(uint8(I42))
                        I43 = imwarp(I42,affine2d(obr.Ttranslate),'OutputView',imref2d(size(obr.template(:,:,1))));figure;imagesc(uint8(I43))

                        Tscale = eye(3)
                        Ttranslate = eye(3)
                        Tflip = eye(3)
                        %}
                        %{
                        figure;imagesc(obr.Corrbest)
                        figure;imagesc(obr.I40)
                        figure;imagesc(uint8(obr.I40))
                        figure;imshowpair(obr.template,uint8(obr.I40))
                        figure; surf(obr.Corrbest)
                    %}


                    %end     
        end 
        function obr = SubPixelAlignment(obr, I40)
            
            if exist('I40', 'var')
                obr.I40 = I40;
            end
            
            %template matching
            %the tail is cut of as it will inprove the template matching, the template has an not well difined tail
            obr.tform3 = affine2d([1, 0, 0; 0, 1, 0; -obr.RemoveTail, 0, 1]);
            obr.I55 = imwarp(obr.I40, obr.tform3, 'FillValues', 0, 'Interp', 'cubic', 'OutputView', imref2d(size(obr.I40(:, :, 1))-[0, 600]));
            
            %final subpixel template matching using 3th party iat software
            
            %AffineMethod = 'affine'
            %AffineMethod = 'translation'
            
            Initialization = [-obr.RemoveTail * obr.Scale3; 0];
            if strcmp(obr.AffineMethod, 'affine')
                Initialization = [eye(2), Initialization];
            end
            
            [obr.tform4, ~, obr.CorCoef] = sng_AllignFish2Template2( ...
                obr.template, ...
                obr.I55, ...
                obr.Scale3, ...
                obr.AffineMethod, ...
                obr.Levels, ...
                obr.Iterations, ...
                Initialization);%missingoutput=I60
            
            %CorCoef= sng_NCCoi(uint8(I65),template);
            
            obr.tform_1234 = affine2d(obr.tform1.T*obr.tform2.T*obr.tform3.T*obr.tform4.T);
            
            obr.Ialigned = imwarp(obr.Icombined, obr.tform_1234, 'FillValues', 255, 'OutputView', imref2d(size(obr.template(:, :, 1))));
        end      
        function obr = fishcolorparameters(obr)
            
            %% extra mean fish color parameter
            %get image - take front(head) - compute histogram - smooth - find second max peak (no background)
            for k3 = 1:3
                Img = obr.Ialigned(:, 700:end, k3);
                obr.MeanFishColor(k3) = mean(Img(:) < obr.BackgroundThreshold(k3));
                h = hist(Img(:), 0:1:255);
                %h2 = smoothdata(h,'gaussian',6);
                h2 = imgaussfilt(h, 1); %gaussian filter function for matlab older than 2017a
                [~, obr.MaxFishColor(k3)] = max(h2(1:floor(obr.BackgroundThreshold(k3))));
            end
        end
        function obr = ClearVars(obr)
            %INPUT/OUTPUT
            %obr.Icombined = [];
            %obr.Ialigned = [];
            
            %BackgroundRemoval
            obr.I20 = [];
            
            %RotationCorrelation
            obr.I21 = [];
            obr.I22 = [];
            obr.I22_1 = [];
            obr.I22_2 = [];
            obr.I22_3 = [];
            obr.I23 = [];
            obr.I24 = [];
            obr.I25 = [];
            obr.I30 = [];
            %
            obr.I40 = [];
            obr.I55 = [];
        end
        function rgoutput = oldoutputfunction(obr)
            
            %{
       figure;imagesc(I10)
       figure;imagesc(uint8(I55));
       figure;imagesc(uint8(I60));
       figure;imagesc(uint8(template));
       figure;imagesc(uint8(CompleteTemplate.MidBrainDistanceMap));
 
       figure;imagesc(uint8(Ialligned));colormap(gray);
       hold on; plot(CompleteTemplate.EyeRegion{1}{1}(:,1),CompleteTemplate.EyeRegion{1}{1}(:,2))
       hold on; plot(CompleteTemplate.EyeRegion{1}{2}(:,1),CompleteTemplate.EyeRegion{1}{2}(:,2))
       plot(CompleteTemplate.MeanMidBrain(:,2),CompleteTemplate.MeanMidBrain(:,1))
       scatter(CompleteTemplate.CenterMidBrain(1),CompleteTemplate.CenterMidBrain(2))
 
       figure;imshowpair(template,uint8(Ialligned))
       figure;imshowpair(template,I51)
 
            %}
            
            
            rgoutput = struct('stage', [], 'substage', [], 'name', []', 'value', []);
            
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Rotation Correlation', 'rotation angle1', obr.theta1});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Rotation Correlation', 'raw angle1', obr.maxangle1});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Rotation Correlation', 'quality 1-peak2/peak1', obr.peaksep1});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Rotation Correlation', 'correlation1', obr.ncorr1});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Rotation Correlation', 'angles1', obr.angles1});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Rotation Correlation', 'rotation angle2', obr.theta3});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Rotation Correlation', 'raw angle2', obr.maxangle2});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Rotation Correlation', 'correlation2', obr.ncorr2});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Rotation Correlation', 'angles2', obr.angles2});
            
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Scale Correlation', 'MaxCorrelationValue', obr.scalebest});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Scale Correlation', 'Tflip', obr.Tflip});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Scale Correlation', 'Tscale', obr.Tscale});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Scale Correlation', 'Ttranslate', obr.Ttranslate});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Scale Correlation', 'BestScaleCorrelationImage', obr.Corrbest});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Scale Correlation', 'MaxCorrelationPlot', obr.maxk});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Scale Correlation', 'FlipIsTwo', obr.csflip});
         
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Rotation Correlation', 'tform1_hor', obr.tform1});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'scale and flip Corr', 'tform2_scale', obr.tform2});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'template matching', 'tform3_tailcut', obr.tform3});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'template matching', 'tform4_trans', obr.tform4});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'final', 'tform_complete', obr.tform_1234});
            
            %rgoutput = sng_StructFill(rgoutput,{'Registration','Rotation Correlation','ref30',ref30});
            %rgoutput = sng_StructFill(rgoutput,{'Registration','scale and flip Corr','ref40',ref40});
            %rgoutput = sng_StructFill(rgoutput,{'Registration','scale and flip Corr','ref50',ref50});
            
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'template matching', 'TemplateCorrelation', obr.CorCoef});
            
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Background Removal', 'Background', obr.BackgroundMask});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Background Removal', 'BackgroundThreshold', obr.BackgroundThreshold});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Background Removal', 'stretchlim_0procent', obr.stretchlim_0procent});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Background Removal', 'stretchlim_1procent', obr.stretchlim_1procent});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Background Removal', 'MeanFishColor', obr.MeanFishColor});
            rgoutput = sng_StructFill(rgoutput, {'Registration', 'Background Removal', 'MaxFishColor', obr.MaxFishColor});
            
            %% show background removal
            %{
       figure;imagesc(I10);axis off equal tight
       set(gca,'position',[0 0 1 1],'units','normalized')
 
       figure;imagesc(I20);axis off equal tight
       set(gca,'position',[0 0 1 1],'units','normalized')
 
       Hr1=histcounts(double(I10),[0:1:255]);
       figure; bar(Hr1,'Facecolor',[0.8,0.4,0.4],'Edgecolor','none','BarWidth',1);
       set(gca,'XLim',[0 256],'YLim',[0 100000],'FontSize',20);
       set(gcf,'Color',[1 1 1])
 
       HrBR1=histcounts(double(I20),[0:1:255]);
       figure; bar(Hr1,'Facecolor',[0.8,0.4,0.4],'Edgecolor','none','BarWidth',1);
       set(gca,'XLim',[0 256],'YLim',[0 100000],'FontSize',20);
       set(gcf,'Color',[1 1 1])
            %}
            
            
            %% Horizontal allignent figures are programmed in sng_HorizontalFish3
            %{
            %}
            
            
            %% Final Template Matching
            %{
       figure;imagesc(uint8(Ialligned))
       set(gca,'position',[0 0 1 1],'units','normalized')
       axis off tight equal
 
       figure; imagesc(boolean(CompleteTemplate.MidBrainDistanceMap))
 
       B = bwboundaries(CompleteTemplate.MidBrainDistanceMap)
       figure;imagesc(uint8(Ialligned))
       hold on;
       plot(B{1}(:,2),B{1}(:,1),'Color',[0 176/255 240/255],'LineWidth',2)
       h = plot(B{2}(:,2),B{2}(:,1),'Color',[0 176/255 240/255],'LineWidth',2)
       plot(CompleteTemplate.MeanMidBrain(:,2),CompleteTemplate.MeanMidBrain(:,1),'Color',[130/255 213/255 247/255],'LineWidth',2)
       set(gca,'position',[0 0 1 1],'units','normalized')
       axis off tight equal
            %}
            
        end
        function obr = All(obr, Icombined)
            
            if exist('Icombined', 'var')
                obr.Icombined = Icombined;
            end
            
            obr = BackgroundRemoval(obr);
            obr = Rotationalalignment(obr);
            obr = Scalealignment(obr);
            obr = SubPixelAlignment(obr);
            obr = fishcolorparameters(obr);
        end
        function ShowRotationAlignment(obr)
            figure; imagesc(obr.I20)
            figure; imagesc(obr.I21)
            figure; imagesc(obr.I22)
            
            figure; imagesc(obr.I22_1)
            figure; imagesc(obr.I22_2)
            for k = 1:numel(obr.I22_3)
                figure; imagesc(uint8(obr.I22_3{k}))
            end
            figure; plot(obr.angles1, obr.ncorr1); xlim([0, 2 * pi])
            
            figure; imagesc(obr.I23); sng_imfix
            figure; imagesc(obr.I24); sng_imfix
            
            figure;
            hold on
            plot(obr.angles1, obr.ncorr1); xlim([0, 2 * pi]);
            plot(obr.maxangle1-pi+obr.angles2, obr.ncorr2); xlim([0, 2 * pi]);
            hold off
            
            figure; imagesc(obr.ref25.XWorldLimits, obr.ref25.YWorldLimits, uint8(obr.I25)); sng_imfix
            %figure;imagesc(ref26.XWorldLimits,ref26.YWorldLimits,uint8(I26));
            figure; imagesc(obr.ref30.XWorldLimits, obr.ref30.YWorldLimits, uint8(obr.I30)); sng_imfix
        end
    end
    
    methods(Hidden) %old functions not used anymore
        function obr = Rotationalalignment_old(obr, I20)
            
            if exist('I20', 'var')
                obr.I20 = I20;
            end
            
            
            %tform1 = sng_HorizontalFish3(I20,'precise'); %translation + rotation
            [obr.tform1, ...
                RotationCorrelationOutput, ...
                obr.I30, ...
                obr.ref30] ...
                = sng_HorizontalFish5( ...
                obr.I20, ...
                obr.AngleSteps1, ...
                obr.Scale1, ...
                obr.AngleRange2, ...
                obr.AngleSteps2, ...
                obr.FinalCenteredYCrop);
            
            
            obr.theta1 = RotationCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'rotation angle1')).value;
            obr.maxangle1 = RotationCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'raw angle1')).value;
            obr.peaksep1 = RotationCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'quality 1-peak2/peak1')).value;
            obr.ncorr1 = RotationCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'correlation1')).value;
            obr.angles1 = RotationCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'angles1')).value;
            obr.theta3 = RotationCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'rotation angle2')).value;
            obr.maxangle2 = RotationCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'raw angle2')).value;
            obr.ncorr2 = RotationCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'correlation2')).value;
            obr.angles2 = RotationCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'angles2')).value;
            
            
        end
        function obr = Scalealignment_old(obr, I30)
            
            if exist('I30', 'var')
                obr.I30 = I30;
            end
            
            %{
       figure;imagesc(I10)
       figure;imagesc(I20)
       figure;imagesc(I30)
       figure;imagesc(ref30.XWorldLimits,ref30.YWorldLimits,uint8(I30));
            %}
            
            %scale and flip
            %scale and leftrigh fitting with template
            %input is needs to be a cropped image for good performance
            
            [obr.tform2, ...
                ScaleCorrelationOutput, ...
                obr.I40] ...
                = sng_ScaleFish4( ...
                double(obr.I30), ...
                double(obr.template), ...
                obr.ScaleRange, ...
                obr.ScaleSteps, ...
                obr.Scale2);
            
            obr.scalebest = ScaleCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'MaxCorrelationValue')).value;
            obr.Tflip = ScaleCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'Tflip')).value;
            obr.Tscale = ScaleCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'Tscale')).value;
            obr.Ttranslate = ScaleCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'Ttranslate')).value;
            obr.Corrbest = ScaleCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'BestScaleCorrelationImage')).value;
            obr.maxk = ScaleCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'MaxCorrelationPlot')).value;
            obr.csflip = ScaleCorrelationOutput(strcmp({obr.RotationCorrelationOutput.name}, 'FlipIsTwo')).value;
     


            
            
            %{
       figure;imagesc(uint8(I40));
           tform_12 = affine2d(tform1.T*tform2.T);
       [I40,ref50] = imwarp(I10,tform_12,'FillValues',255,'OutputView',imref2d(size(template(:,:,1))));
 
       figure;imagesc(uint8(I40));
       figure;imagesc(ref50.XWorldLimits,ref50.YWorldLimits,uint8(I40));
 
       figure;imagesc(ScaleCorrelationOutput(5).value)
 
 
                %}
        end       
    end
    
    methods(Static)
        function [theta, varargout] = NormalizedRotationCorrelation(I30, angles1, centerofmass, varargin)
            %
            %Version sng_NormalizedRotationCorrelation3 20170616
            %   change in/output for better use in rgoutput
            %Version sng_NormalizedRotationCorrelation4 20180611
            %   similar to version 3 but embedded in SNGAlignment class
            %
            %shifts image to center and calculated normalized correlation between the
            %image and a flipped and rotated version of itself
            %varargin{1} the rotated images specified with number
            %
            %
            %Example [theta,maxangledeg,peaksep,ncorr] = sng_NormalizedRotationCorrelation4(I30,angles,[xcom2,ycom2])
            %Example [theta,maxangledeg,peaksep,ncorr,I40,I50] = sng_NormalizedRotationCorrelation4(I30,angles1,[xcom2,ycom2])
            %Example [theta,maxangledeg,peaksep,ncorr,I40,I50,I61] = sng_NormalizedRotationCorrelation4(Img2,angles1,[xcom2,ycom2],1:10:100)
            %
            %
            
            %{
               anglesteps = 100
               centerofmass = [xcom2,ycom2]
               I30 = Img2;
               nargin = 3
            %}
            
            if nargin == 4
                I61 = cell(numel(varargin{1}), 1);
            end
            
            shift = (fliplr(size(I30(:, :, 1))/2) - centerofmass);
            
            %translates and crops 1 times the shift ->maybe more accurate
            I40 = imtranslate(I30, shift, 'OutputView', 'same');
            %no translation and crop 2 times the shift -> smaller images
            %I40 = imcrop(I30, [abs(shift),-abs(shift)] - [shift, shift] + [0 0 size(I30(:,:,1)));
            
            %make circle mask
            cm = sng_MaxCircleMask(I40);
            
            
            %{
   figure;imagesc(c_mask)
   figure;imagesc(I40)
   figure;imagesc(I40.*uint8(cm))
            %}
            
            I50 = fliplr(I40);
            
            %{
   figure;imagesc(I40norm);colorbar
   figure;imagesc(inproductI40)
   figure;imagesc(I50.*uint8(cm))
 
            %}
            
            ncorr = zeros(numel(angles1), 1);
            jj = 1;
            for j = 1:numel(angles1)
                %I60 = imrotate(I50,angles(j),'bilinear','crop');
                
                I60 = sng_imrotate(I50, angles1(j));
                ncorr(j) = sng_NCCoi(I40, I60, cm);
                
                if nargin == 4 && ismember(j, varargin{1})
                    I61{jj} = I60;
                    jj = jj + 1;
                end
            end
            
            [~, maxcorr] = (max(ncorr));
            theta = -(angles1(maxcorr) / 2 + (pi / 2));
            
            
            I60 = sng_imrotate(I50, angles1(maxcorr)); %the best angle
            
            
            %% calculates the relative difference of the first and second maxima
            % to give an idication how sure we are of the maximum peak
            ncorrmax = imregionalmax(ncorr) .* ncorr;
            corrpeaks = sort(nonzeros(ncorrmax), 'descend');
            
            if numel(corrpeaks) >= 2
                peaksep = 1 - abs(corrpeaks(2)/corrpeaks(1));
            else
                peaksep = 1 - abs(min(ncorr)/corrpeaks(1));
            end
            
            maxangle = angles1(maxcorr);
            
            if nargout >= 2; varargout{1} = maxangle; end %maximum angle in degrees
            if nargout >= 3; varargout{2} = peaksep; end %relative difference two highes peaks
            if nargout >= 4; varargout{3} = ncorr; end
            if nargout >= 5; varargout{4} = I40; end
            if nargout >= 6; varargout{5} = I50; end
            if nargout >= 7; varargout{6} = I61; end
            if nargout >= 8; varargout{7} = I60; end
            
            %{
   figure;plot(ncorr)
   for k = 1:numel(I61)
       figure;imagesc(uint8(I61{k}))
   end
 
 
            %}
            
            
        end
    end
end