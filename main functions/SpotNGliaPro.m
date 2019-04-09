classdef SpotNGliaPro < SpotNGlia
    
    %TODO move all show methods to SpotNGliaPro and add to SpotNGlia-private
    %TODO make a method in SpotNGliaPro that copies all props from an SpotNGlia instance to a SpotNGliaPro instance
    
    
    properties
        %obj = []; %SpotNGlia Object
        fsxy = [6, 5];
        exportit = false;
        %subject = []
        %SavePath = []
        
        
    end
    
    methods
        
        function obj = SpotNGliaPro(obj) %, subject, exportit, fsxy)
            
            
            obj = obj@SpotNGlia;
            
            %sobj.obj = obj;
            %sobj.SavePath = [obj.SavePath, filesep, 'Images'];
            
            %             if exist('subject','var')
            %                 sobj.subject = subject;
            %             else
            %                 sobj.subject = [];
            %             end
            %
            %             if exist('fsxy', 'var')
            %                 sobj.fsxy = fsxy;
            %             end
            %
            %             if exist('exportit', 'var')
            %                 sobj.exportit = exportit;
            %             end
            %
            %             %ShowFishDiscrimination1 - histogram shows correlation coefficient adjacent images
            %             %ShowFishDiscrimination2 - scatterplot of translation
            %
            %             if ismember(1,sobj.subject); sobj.ShowFishDiscrimination1;end
            %             if ismember(2,sobj.subject); sobj.ShowFishDiscrimination2;end
            %             if ismember(3,sobj.subject); sobj.ShowFishDiscrimination3;end
            %             if ismember(4,sobj.subject); sobj.ShowFishDiscrimination4;end
        end
        
        
        %{
function ShowColorChannelCorrection(sobj)
 
 
         Im2{1} = imread([sobj.obj.,'/set2-Image053.tif']);
 
 
             a = [sobj.obj.ImageInfo.corcoef];
             d = [sobj.obj.ImageInfo.CorNextStack];
 
             %histogram shows correlation coefficientadjacent images of the same fish
             [h2, g2] = sobj.setfigax1;
             histogram(a(d == 0), 0.90:0.005:1);
             hold on
             histogram(a(d == 1), 0.90:0.005:1);
             sobj.setfigax2(g2);
             set(g2, 'XLim', [0.90, 1]);
             set(g2, 'YLim', [0.5, 300])
             set(g2, 'YScale', 'log');
             set(g2, 'YTick', [1, 2, 5, 10 ,100 200]);
             set(g2, 'YTickLabels', [1, 2, 5, 10, {'10^2'},{''}]);
             set(g2, 'YMinorGrid','off');
 
             set(g2, 'XTick',[0.9:0.02:1]);
             %set(g2, 'XTickLabels', [0.90,{''},0.92,{''},0.94,{''},0.96,{''},0.98,{''},1]);
 
 
             g2.XLabel.String = 'correlation coefficient';
             g2.YLabel.String = 'counts';
             sobj.realsizeandsave(h2);
         end
        %}
        
        
    end
    
    methods %all Show methods
        function show(obj, subject, exportit, fsxy, var)
            
            if ~exist('subject', 'var')
                subject = [];
            end
            if exist('exportit', 'var')
                obj.exportit = exportit;
            end
            if exist('fsxy', 'var') && ~isempty(fsxy)
                obj.fsxy = fsxy;
            end
            
            %showFishDiscrimination1 - histogram shows correlation coefficient adjacent images
            %showFishDiscrimination2 - scatterplot of translation
            
            if ismember(01, subject); obj.showFishDiscrimination1; end
            if ismember(02, subject); obj.showFishDiscrimination2; end
            if ismember(03, subject); obj.showFishDiscrimination3; end
            if ismember(04, subject); obj.showFishDiscrimination4; end
            if ismember(11, subject); obj.showBackgroundRemoval(var); end
            
            if ismember(a1, subject); obj.showRotationPreproc(var); end
            if ismember(a2, subject); obj.showFishRotationplot(var); end
            if ismember(a3, subject); obj.showRotatedFishes(var); end
            if ismember(a4, subject); obj.showScalePlot(var); end
            
            if ismember(b1, subject); obj.showMidBrainContours; end
            if ismember(b2, subject); obj.showForeBrainContours; end
            if ismember(b3, subject); obj.showMidBrainBand; end
            if ismember(b4, subject); obj.showMidBrainBand2; end
            if ismember(b5, subject); obj.showMidBrainBandIntersection; end
            
            
        end
        function showImageTemplate(obj, fn)
            
            %LOAD/COMPUTE VALUES
            %Image = obj.CompleteTemplate.Template;
            SIZE = size(Image);
            
            %IMAGE ENHANCEMENT
            
            %SET FIGURE SIZE
            %obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 4; %set height
            %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
            obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj);
            
            %PLOT IMAGE
            imagesc(Image)
            
            g.Position = [0, 0, 1, 1]; axis off;
            
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['defaultImage']);
            
        end
        function showPlotTemplate(obj, fn)
            %output "obj" is to store Registration parameters in object
            %which can be used for showMethods on the same fish
            
            obj.fsxy = ([5, 4]);
            
            %LOAD/COMPUTE VALUES
            %PresetShowRegistration(obj,fn);
            [h, g] = setfigax1(obj); %create figure with axis with real-size obj.fsxy
            
            %PLOT FIGURES
            
            %histogram(a(d == 0), 0.90:0.005:1);
            %hold on
            %histogram(a(d == 1), 0.90:0.005:1);
            
            
            setfigax2(obj, g); %preset axis with standard fond and size
            
            %SET AXIS
            
            %set(g, 'XLim', [0.90, 1]);
            %set(g, 'XTick', [0.9:0.02:1]);
            %set(g, 'XTickLabels', [0.90,{''},0.92,{''},0.94,{''},0.96,{''},0.98,{''},1]);
            %set(g, 'XTickLabels', [0,{'1/2\pi'},{'\pi'},{'3/2\pi'},{'2\pi'}]);
            %g.XLabel.String = 'correlation coefficient';
            
            %set(g, 'YLim', [0.5, 300]);
            %set(g, 'YScale', 'log');
            %set(g, 'YTick', [1, 2, 5, 10, 100, 200]);
            %set(g, 'YTickLabels', [1, 2, 5, 10, {'10^2'}, {''}]);
            %set(g, 'YMinorGrid', 'off');
            %g.YLabel.String = 'counts';
            
            
            realsizeandsave(obj, h, 'SaveName');
        end
        function showFishDiscrimination1(obj)
            obj.fsxy = [6, 5];
            
            a = [obj.ImageInfo.corcoef];
            d = [obj.ImageInfo.CorNextStack];
            
            %histogram shows correlation coefficientadjacent images of the same fish
            [h2, g2] = setfigax1(obj);
            histogram(a(d == 0), 0.90:0.005:1);
            hold on
            histogram(a(d == 1), 0.90:0.005:1);
            setfigax2(obj, g2);
            set(g2, 'XLim', [0.90, 1]);
            set(g2, 'YLim', [0.5, 300]);
            set(g2, 'YScale', 'log');
            set(g2, 'YTick', [1, 2, 5, 10, 100, 200]);
            set(g2, 'YTickLabels', [1, 2, 5, 10, {'10^2'}, {''}]);
            set(g2, 'YMinorGrid', 'off');
            
            set(g2, 'XTick', 0.9:0.02:1);
            %set(g2, 'XTickLabels', [0.90,{''},0.92,{''},0.94,{''},0.96,{''},0.98,{''},1]);
            
            
            g2.XLabel.String = 'correlation coefficient';
            g2.YLabel.String = 'counts';
            realsizeandsave(obj, h2, 'FishDiscrimination1');
        end
        function showFishDiscrimination2(obj)
            obj.fsxy = [6, 5];
            
            a = [obj.ImageInfo.corcoef];
            d = [obj.ImageInfo.CorNextStack];
            
            %histogram shows correlation coefficient adjacent images containing different fish
            [h1, g1] = setfigax1(obj);
            histogram(a(d == 0), 0:0.05:1);
            hold on
            hh = histogram(a(d == 1), 0:0.05:1);
            set(hh, 'FaceColor', [0.8500, 0.3250, 0.0980]);
            set(g1, 'YScale', 'log')
            set(g1, 'YLim', [0.5, 300])
            set(g1, 'XLim', [0, 1])
            set(g1, 'YTick', [1, 2, 5, 10, 100, 200])
            set(g1, 'YTickLabels', [1, 2, 5, 10, {'10^2'}, {''}])
            set(g1, 'YMinorGrid', 'off')
            
            set(g1, 'XTick', 0:0.2:1)
            
            
            setfigax2(obj, g1)
            g1.XLabel.String = 'correlation coefficient';
            g1.YLabel.String = 'counts';
            realsizeandsave(obj, h1, 'FishDiscrimination2')
            
            %{
                           %histogram of modulus shift around the mean of similar images
                           [h3, g3] = setfigax1;
                               histogram(c(d == 0),0:0.16:1.6);
                               hold on
                               hh=histogram(c(d == 1),0:0.16:2.5);
                               set(hh,'FaceColor',[0.8500    0.3250    0.0980]);
                               set(gca,'YLim',[0.6,301]);
                               set(gca,'YScale','log');
                               set(gca,'YTick',[1,10,100]);
                               set(gca,'YTickLabels',[1,10,{'10^2'}])
                           setfigax2(g3)
                           g3.XLabel.String = 'translation modulus [pix]';
                           g3.YLabel.String = 'counts';
                           realsizeandsave(h3)
 
                    %histogram of angular shift zoomed out
 
 
                           [h3, g3] = setfigax1;
                               h=histogram(c(d == 1),[-pi:pi/20:pi])
                               set(h,'FaceColor',[0.8500    0.3250    0.0980])
                               hold on
                               h=histogram(c(d == 0),[-pi:pi/20:pi])
            %}
        end
        function showFishDiscrimination3(obj)
            %scatterplot of translation zoomed out
            obj.fsxy = [6, 5];
            
            d = [obj.ImageInfo.CorNextStack];
            e = reshape([obj.ImageInfo.warp], 2, size(obj.ImageInfo, 1))';
            
            [h4, g4] = setfigax1(obj);
            scatter(e((d == 0), 1), e((d == 0), 2))
            hold on
            scatter(e((d == 1), 1), e((d == 1), 2))
            xlim([-20, 20])
            ylim([-20, 20])
            set(g4, 'XTick', -20:10:20)
            set(g4, 'YTick', -20:10:20)
            
            setfigax2(obj, g4)
            %make axis equal but leftdown axis intersection at same position
            pos = get(gca, 'Position');
            set(gca, 'Position', [pos(1:2), min(pos(3:4)), min(pos(3:4))])
            
            g4.XLabel.String = 'translation x-axis [pix]';
            g4.YLabel.String = 'translation y-axis [pix]';
            realsizeandsave(obj, h4, 'FishDiscrimination3')
        end
        function showFishDiscrimination4(obj)
            %scatterplot of translation zoomed in
            obj.fsxy = [6, 5];
            b = [obj.ImageInfo.moduluswarp];
            c = [obj.ImageInfo.anglewarp];
            d = [obj.ImageInfo.CorNextStack];
            e = reshape([obj.ImageInfo.warp], 2, size(obj.ImageInfo, 1))';
            
            %compute the location of the outer angles
            c1 = c(d == 0);
            c2 = c1(c1 < 0);
            [c3min, n1] = min(c2);
            [c3max, n2] = max(c2);
            e1 = e((d == 0), :);
            e2 = e1((c1 < 0), :);
            %maximum modulus for circle
            b1 = b(d == 0);
            b2 = max(b1);
            
            
            %maxdeg = round(360*-0.15/(2*pi),1
            %mindeg = round(360*-0.60/(2*pi),1)
            
            
            [h4, g4] = setfigax1(obj);
            scatter(e((d == 0), 1), e((d == 0), 2))
            hold on
            scatter(e((d == 1), 1), e((d == 1), 2))
            
            xlim([-1.6, 1.6])
            ylim([-1.6, 1.6])
            setfigax2(obj, g4)
            
            %make axis equal but leftdown axis intersection at same position
            pos = get(gca, 'Position');
            set(gca, 'Position', [pos(1:2), min(pos(3:4)), min(pos(3:4))])
            
            
            %gets position of plot in figure, set position
            %set(gca,'units','centimeters')
            %ppos = plotboxpos;
            %set(gca,'Position',[get(gca,'Position') - [ppos(1:2),0,0] + [1.2 1.2,0,0]]);
            
            
            g4.XLabel.String = 'translation x-axis [pix]';
            g4.YLabel.String = 'translation y-axis [pix]';
            
            %add angle lines
            line(10*[-e2(n1, 1), e2(n1, 1)], 10*[-e2(n1, 2), e2(n1, 2)], 'Color', [0, 0, 0], 'LineStyle', '--');
            line(10*[-e2(n2, 1), e2(n2, 1)], 10*[-e2(n2, 2), e2(n2, 2)], 'Color', [0, 0, 0], 'LineStyle', '--');
            r = b2; %radius
            %add circle
            pos = [-r, -r, 2 * r, 2 * r];
            rectangle('Position', pos, 'Curvature', [1, 1], 'LineStyle', '--');
            %add annotations
            maxdeg = round(360*c3max/(2 * pi), 1);
            mindeg = round(360*c3min/(2 * pi), 1);
            text(0, 1.3, ['r = ', num2str(round(r, 1))], 'FontName', 'arial', 'FontSize', 8, 'HorizontalAlignment', 'Center');
            text(1.65, -0.32, [num2str(maxdeg), '\circ'], 'FontName', 'arial', 'FontSize', 8);
            text(1.65, -0.95, [num2str(mindeg), '\circ'], 'FontName', 'arial', 'FontSize', 8);
            
            realsizeandsave(obj, h4, 'FishDiscrimination4');
            
        end
        function showOriginal(obj, fn)
            
            imageSlice = cell(1, obj.StackInfo(fn).stacksize); %preallocate for every new slice
            for iSlice = 1:obj.StackInfo(fn).stacksize
                imageSlice{iSlice} = imread([obj.FishPath, filesep, obj.StackInfo(fn).imagenames{iSlice}]);
                imageSlice{iSlice} = im2uint8(imageSlice{iSlice}(:, :, 1:3));
                
                
                %LOAD/COMPUTE VALUES
                %Image = obj.CompleteTemplate.Template;
                SIZE = size(imageSlice{iSlice});
                
                %IMAGE ENHANCEMENT
                
                %SET FIGURE SIZE
                %obj.fsxy(1) = 5; %set width
                obj.fsxy(2) = 3; %set height
                %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
                obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
                
                %CREATE FIGURE WITH AXIS
                [h, g] = setfigax1(obj);
                
                %PLOT IMAGE
                imagesc(imageSlice{iSlice})
                
                g.Position = [0, 0, 1, 1]; axis off;
                
                %EXPORT IMAGE
                realsizeandsave(obj, h, ['OriginalImageFish', num2str(fn), '_', num2str(iSlice)]);
                
            end
            
        end
        function showRGBCorrection1(obj, fn)
            
            %LOAD/COMPUTE VALUES
            if isempty(obj.PreprocessingObject(fn).bandFiltered)
                obj.PreprocessingObject(fn) = obj.PreprocessingObject(fn).loadImageSlices;
                obj.PreprocessingObject(fn) = obj.PreprocessingObject(fn).rgbBandpassFilter;
            end
            Image = obj.PreprocessingObject(fn).bandFiltered{obj.PreprocessingObject.sliceNo};
            
            %IMAGE ENHANCEMENT
            Image2 = abs(0.15*Image);
            
            SIZE = size(Image);
            
            %SET FIGURE SIZE
            %obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 10; %set height
            %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
            obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj);
            
            %PLOT IMAGE
            imagesc(Image2)
            
            g.Position = [0, 0, 1, 1]; axis off;
            
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['BandFiltered_1']);
            
        end
        function showRGBCorrection2(obj, fn)
            
            %LOAD/COMPUTE VALUES
            if isempty(obj.PreprocessingObject(fn).bandFiltered)
                obj.PreprocessingObject(fn) = obj.PreprocessingObject(fn).loadImageSlices;
                obj.PreprocessingObject(fn) = obj.PreprocessingObject(fn).rgbBandpassFilter;
            end
            Image = obj.PreprocessingObject(fn).filteredImage{obj.PreprocessingObject.sliceNo};
            
            
            %IMAGE ENHANCEMENT
            Image2 = abs(0.15*Image);
            
            SIZE = size(Image);
            
            %SET FIGURE SIZE
            %obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 8; %set height
            %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
            obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj);
            
            %PLOT IMAGE
            imagesc(Image2)
            
            g.Position = [0, 0, 1, 1]; axis off;
            
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['filteredImage']);
        end
        function showRGBCorrectionZoom(obj, fn)
            
            %LOAD/COMPUTE VALUES
            if isempty(obj.PreprocessingObject(fn).bandFiltered)
                obj.PreprocessingObject(fn) = obj.PreprocessingObject(fn).loadImageSlices;
                obj.PreprocessingObject(fn) = obj.PreprocessingObject(fn).rgbBandpassFilter;
            end
            Image = obj.PreprocessingObject(fn).filteredImage{obj.PreprocessingObject.sliceNo};
            
            
            %IMAGE ENHANCEMENT
            Image2 = abs(0.15*Image);
            
            SIZE = size(Image);
            
            Image2 = Image2(1:SIZE(2)/1, 1:end);
            
            %zoom parameters
            zoomn = 20;
            zoompoint = [1404, 643];
            
            %SET FIGURE SIZE
            %obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 8; %set height
            %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
            obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj);
            
            %PLOT IMAGE
            imagesc(Image2)
            
            sng_zoom(zoomn, zoompoint, SIZE, g)
            
            g.Position = [0, 0, 1, 1]; axis off;
            
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['filteredImage']);
        end
        
        function showBackgroundRemoval(obj, fn)
            
            %sng_zfinputAssign(obj.ZFParameters, 'Registration')
            obj.SngInputParameters.assign('Registration')
            
            
            if fn <= 0
                %open chosen images from separate folder
                PathName = '/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/presentation nov-2016/remove background test/'
                
                %find dropbox folder
                Dropboxfolder = obj.SourcePath;
                i = 0;
                b = '';
                while ~strcmp(b, 'Dropbox') && (i <= 5) %removes folder uptil dropbox folder is found
                    [Dropboxfolder, b] = fileparts(Dropboxfolder);
                    i = i + 1;
                end
                
                PathName = [Dropboxfolder, filesep, 'Dropbox', filesep, 'Images for Report', filesep];
                FileName1 = 'image0015.tif';
                FileName2 = 'image0049.tif';
                FileName3 = 'Image56.tif';
                
                if fn == -1
                    Img = imread([PathName, FileName1]);
                elseif fn == -2
                    Img = imread([PathName, FileName2]);
                elseif fn == -3
                    Img = imread([PathName, FileName3]);
                end
                Img = Img(:, :, 1:3);
                Img = im2uint8(Img);
            else
                CombinedFish = imread([obj.SavePath, filesep, 'CombinedFish', filesep, obj.StackInfo(fn).stackname, '.tif']);
                %[AlignedFishTemp, RegistrationInfo] = AlignmentSNG(CombinedFish, obj.CompleteTemplate, obj.ZFParameters);
                Img = CombinedFish;
            end
            
            
            %[I20, BackgroundMask, BackgroundThreshold] = sng_RemoveBackgroundColor3(Img, Method, Smooth, ChannelMethod);
            
            for j = 1:3
                [Img2(:, :, j), thr(j), histo{j}, smoothhisto{j}, Tables{j}] = sng_RemoveBackground(Img(:, :, j), Method, Smooth);
            end
            %coordinates form lines and arrows
            tc = Tables{1}.Threshold; %for arrow
            hc = Tables{1}.Halfway;
            pc = Tables{1}.Peak; %for line
            lc = Tables{1}.LastValue;
            textc = (tc - hc) / 2 + hc; %for halfway arro
            
            %{
                        %histogram of original image
                        Hr1 = histcounts(double(Img(:,:,1)), [0:1:255]);
                        figure; bar(Hr1, 'Facecolor', [0.5, 0.5, 0.5], 'Edgecolor', 'none', 'BarWidth', 1);
                        set(gca, 'XLim', [0, 256], 'YLim', [0, 100000], 'FontSize', 20);
                        set(gcf, 'Color', [1, 1, 1])
            %}
            
            %{
                        %histogram of image with removed background
                        Hr1 = histcounts(double(Img2(:,:,1)), [0:1:255]);
                        figure; bar(Hr1, 'Facecolor', [0.5, 0.5, 0.5], 'Edgecolor', 'none', 'BarWidth', 1);
                        set(gca, 'XLim', [0, 256], 'YLim', [0, 100000], 'FontSize', 20);
                        set(gcf, 'Color', [1, 1, 1])
            %}
            
            %histogram of the 50pixel boundary
            [f1, a1] = setfigax1(obj);
            
            bar(fliplr(histo{1}), 'Facecolor', [0.5, 0.5, 0.5], 'Edgecolor', 'none', 'BarWidth', 1);
            set(gcf, 'Color', [1, 1, 1])
            
            hold on
            %line
            plot(255-[pc(1), lc(1)], [pc(2), lc(2)], 'Color', [0, 0, 0])
            
            %text
            th = text(255-textc(1), textc(2), 'd');
            set(th, 'VerticalAlignment', 'top');
            set(th, 'HorizontalAlignment', 'right');
            
            setfigax2(obj, a1);
            
            set(a1, 'YTick', [5000, 10000]);
            set(a1, 'YTickLabels', [{'5k'}, {'10k'}]);
            set(a1, 'XLim', [0, 255])
            set(a1, 'YLim', [0, pc(2) * 1.1])
            set(a1, 'XTick', [0, 255]);
            set(a1, 'XTickLabels', [{'black'}, {'white'}]);
            %g1.XLabel.String = 'correlation coefficient';
            a1.YLabel.String = 'counts';
            
            set(a1, 'Units', 'Normalized');
            pbp = plotboxpos(a1);
            transformx = pbp(3) / (diff(get(a1, 'XLim')));
            transformy = pbp(4) / (diff(get(a1, 'YLim')));
            tcn(1) = ((255 - tc(1)) * transformx) + pbp(1);
            hcn(1) = ((255 - hc(1)) * transformx) + pbp(1);
            tcn(2) = (tc(2) * transformy) + pbp(2);
            hcn(2) = (hc(2) * transformy) + pbp(2);
            annotation('doublearrow', [tcn(1), hcn(1)], [tcn(2), hcn(2)]);
            
            realsizeandsave(obj, f1, 'BackGroundRemoval');
        end
        function PresetShowRegistration(obj, fn)
            
            LoadParameters(obj, 'RegObject');
            
            if isempty(obj.RegObject(fn).Icombined)
                CombinedFish = imread([obj.SavePath, filesep, 'CombinedFish', filesep, obj.StackInfo(fn).stackname, '.tif']);
                obj.RegObject(fn).Icombined = CombinedFish;
            end
            if isempty(obj.RegObject(fn).I20)
                obj.RegObject(fn) = BackgroundRemoval(obj.RegObject(fn));
            end
            if isempty(obj.RegObject(fn).I30)
                obj.RegObject(fn) = Rotationalalignment(obj.RegObject(fn));
            end
            if isempty(obj.RegObject(fn).I40)
                obj.RegObject(fn) = Scalealignment(obj.RegObject(fn));
            end
            if isempty(obj.RegObject(fn).Ialigned)
                obj.RegObject(fn) = SubPixelAlignment(obj.RegObject(fn));
            end
        end
        function showRotationPreproc(obj, fn)
            
            obj.fsxy(1) = [6]; %give only width
            
            
            %LOAD/COMPUTE VALUES
            PresetShowRegistration(obj, fn);
            
            I{1} = obj.RegObject(fn).I20;
            I{2} = obj.RegObject(fn).I21;
            I{3} = obj.RegObject(fn).I22;
            I{4} = obj.RegObject(fn).I22_1;
            I{5} = obj.RegObject(fn).I22_2;
            I{6} = obj.RegObject(fn).I22_4;
            
            I{7} = obj.RegObject(fn).I23; %global rotation
            I{8} = obj.RegObject(fn).I24;
            I{9} = obj.RegObject(fn).I25; %precise rotation
            I{10} = obj.RegObject(fn).I30;
            
            
            for k = 1:numel(I)
                
                sz = size(I{k});
                obj.fsxy(2) = sz(1) / sz(2) * obj.fsxy(1); %height is dependent on width
                
                [h, g] = setfigax1(obj); %create figure with axis with real-size obj.fsxy
                
                %PLOT IMAGE
                imagesc(I{k})
                
                g.Position = [0, 0, 1, 1]; axis off;
                
                realsizeandsave(obj, h, ['RotationImages_fish', num2str(fn), '_', num2str(k)]);
                
            end
        end
        function showFishRotationplot(obj, fn)
            
            obj.fsxy = [10, 8];
            
            %LOAD/COMPUTE VALUES
            PresetShowRegistration(obj, fn);
            
            [h, g] = setfigax1(obj); %create figure with axis with real-size obj.fsxy
            
            %PLOT FIGURES
            plot(obj.RegObject(fn).angles1, obj.RegObject(fn).ncorr1);
            
            setfigax2(obj, g); %preset axis with standard fond and size
            
            %SET AXIS
            set(g, 'XLim', [0, 2 * pi]);
            set(g, 'YLim', [-0.4, 1]);
            %set(g, 'YScale', 'log');
            %set(g, 'YTick', [1, 2, 5, 10, 100, 200]);
            %set(g, 'YTickLabels', [1, 2, 5, 10, {'10^2'}, {''}]);
            %set(g, 'YMinorGrid', 'off');
            set(g, 'XTick', 2*pi*linspace(0, 1, 5));
            set(g, 'XTickLabels', [0, 90, 180, 270, 360]);
            g.XLabel.String = 'rotation angle';
            g.YLabel.String = 'correlation coefficient';
            %}
            
            realsizeandsave(obj, h, 'Rotationplot');
            
        end
        function showRotatedFishes(obj, fn)
            
            obj.fsxy(1) = [6];
            
            %LOAD/COMPUTE VALUES
            PresetShowRegistration(obj, fn);
            
            for k = 1:numel(obj.RegObject(fn).I22_3)
                
                sz = size(obj.RegObject(fn).I22_3{k});
                obj.fsxy(2) = sz(1) / sz(2) * obj.fsxy(1);
                
                [h, g] = setfigax1(obj); %create figure with axis with real-size obj.fsxy
                
                %PLOT FIGURES
                imagesc(uint8(obj.RegObject(fn).I22_3{k}));
                axis off;
                
                %SET AXIS
                g.Position = [0, 0, 1, 1]; axis off;
                
                angle = uint16(obj.RegObject(fn).angles1(obj.RegObject(fn).RotationIndices)*(180 / pi));
                realsizeandsave(obj, h, ['RotatedReflection_fish', num2str(fn), '_ang', num2str(angle(k))]);
            end
        end
        function showScalePlot(obj, fn)
            %output "obj" is to store Registration parameters in object
            %which can be used for showMethods on the same fish
            
            %obj.fsxy = ([5, 4]);
            
            %LOAD/COMPUTE VALUES
            PresetShowRegistration(obj, fn);
            [h, g] = setfigax1(obj); %create figure with axis with real-size obj.fsxy
            
            %PLOT FIGURES
            xrange = obj.RegObject(fn).ScaleRange;
            xaxs = linspace(xrange(1), xrange(2), 200);
            plot(xaxs, obj.RegObject(fn).maxk(1, :))
            hold on; plot(xaxs, obj.RegObject(fn).maxk(2, :))
            
            
            setfigax2(obj, g); %preset axis with standard fond and size
            
            %SET AXIS
            set(g, 'XLim', xrange);
            set(g, 'XTick', xrange(1):0.2:xrange(2));
            %set(g, 'XTickLabels', [0 90 180 270 360]);
            g.XLabel.String = 'Scale factor';
            
            ylimtemp = get(g, 'YLim');
            set(g, 'YLim', [ylimtemp(1), 1]);
            %set(g, 'YScale', 'log');
            %set(g, 'YTick', [1, 2, 5, 10, 100, 200]);
            %set(g, 'YTickLabels', [1, 2, 5, 10, {'10^2'}, {''}]);
            set(g, 'YMinorGrid', 'off');
            g.YLabel.String = 'correlation coefficient';
            %}
            
            realsizeandsave(obj, h, 'ScaleCorrelation');
        end
        
        %brain template images
        function showMidBrainContours(obj)
            
            %LOAD/COMPUTE VALUES
            templateImage = obj.CompleteTemplate.Template;
            nMidbrainXCoordList = obj.CompleteTemplate.midbrainXCoordList;
            nMidbrainYCoordList = obj.CompleteTemplate.midbrainYCoordList;
            SIZE = obj.CompleteTemplate.Size;
            
            %SET FIGURE SIZE
            %obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 4; %set height
            %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
            obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj);
            
            %PLOT IMAGE
            nFishes = numel(nMidbrainXCoordList);
            imagesc(templateImage)
            for iFish = 1:nFishes
                hold on
                plot(nMidbrainXCoordList{iFish}, nMidbrainYCoordList{iFish})
            end
            
            g.Position = [0, 0, 1, 1]; axis off;
            
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['MidBrainContours']);
        end
        function showForeBrainContours(obj)
            
            %LOAD/COMPUTE VALUES
            templateImage = obj.CompleteTemplate.Template;
            nForebrainXCoordList = obj.CompleteTemplate.forebrainXCoordList;
            nForebrainYCoordList = obj.CompleteTemplate.forebrainYCoordList;
            SIZE = obj.CompleteTemplate.Size;
            
            %SET FIGURE SIZE
            %obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 4; %set height
            %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
            obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj); %create figure with axis with real-size obj.fsxy
            
            %PLOT IMAGE
            nFishes = numel(nForebrainXCoordList);
            imagesc(templateImage)
            for iFish = 1:nFishes
                hold on
                plot(nForebrainXCoordList{iFish}, nForebrainYCoordList{iFish})
            end
            
            g.Position = [0, 0, 1, 1]; axis off;
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['ForeBrainContours']); %export with savename
        end
        function showMidBrainBand(obj)
            
            %LOAD/COMPUTE VALUES
            Image = obj.CompleteTemplate.Template;
            
            contourCenter = obj.CompleteTemplate.MeanMidBrain;
            contourInner = obj.CompleteTemplate.midbrainInnerContour;
            contourOuter = obj.CompleteTemplate.midbrainOuterContour;
            centerPoint = obj.CompleteTemplate.CenterMidBrain;
            SIZE = obj.CompleteTemplate.Size;
            
            
            %zoom out to show whole rectangle image
            addToImage = round(0.04*SIZE); %crusial value determining the zoom
            newImage = uint8(255*ones(SIZE+2*addToImage));
            newImage(1+addToImage(1):SIZE(1)+addToImage(1), ...
                1+addToImage(2):SIZE(2)+addToImage(2), 1:3) = Image;
            
            %add to coordinates to fit zoom
            contourCenter = contourCenter + repmat(addToImage(1:2), size(contourCenter, 1), 1);
            contourInner = contourInner + repmat(addToImage(1:2), size(contourInner, 1), 1);
            contourOuter = contourOuter + repmat(addToImage(1:2), size(contourOuter, 1), 1);
            centerPoint = centerPoint + [addToImage(2), addToImage(1)];
            SIZE = size(newImage);
            
            %figure;imagesc(Image)
            %figure;imagesc(uint8(newImage))
            
            
            %SET FIGURE SIZE
            %obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 4; %set height
            %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
            obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj); %create figure with axis with real-size obj.fsxy
            
            %PLOT IMAGE
            
            imagesc(newImage)
            hold on
            plot(contourCenter(:, 2), contourCenter(:, 1), 'Color', [0, 166, 214]/255, 'LineWidth', 2)
            plot(contourInner(:, 2), contourInner(:, 1), 'Color', [110, 187, 213]/255, 'LineWidth', 1)
            plot(contourOuter(:, 2), contourOuter(:, 1), 'Color', [110, 187, 213]/255, 'LineWidth', 1)
            
            plot(centerPoint(1), centerPoint(2), 'Color', [226, 26, 26]/255, 'LineWidth', 2, 'Marker', 'X')
            
            g.Position = [0, 0, 1, 1]; axis off;
            
            
            %the rectangle is not exacly 1000 what it should be but it is just
            %for showing this figure that it is made a bit smaller to be able to
            %see whole square
            %rectangle('Position',[centerPoint-480 960 960],'EdgeColor',[226 26 26]/255,'LineWidth' , 1);
            
            rectangle('Position', [centerPoint - 500, 1000, 1000], 'EdgeColor', [226, 26, 26]/255, 'LineWidth', 1);
            
            
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['MidBrainBand']); %export with savename
        end
        function showMidBrainBand2(obj)
            %%%%%%%%%follow with other brain images
            %obj.LoadParameters('BrainSegmentationInfo')
            
            %LOAD/COMPUTE VALUES
            Images{1} = obj.CompleteTemplate.squareMidbrainBand;
            Images{2} = obj.CompleteTemplate.polarMidbrainBand;
            Images{3} = obj.CompleteTemplate.polarMidbrainBandWithGaussian;
            Images{4} = obj.CompleteTemplate.polarMidbrainBandWithGaussianWithBlurr;
            
            for iImage = 1:numel(Images)
                SIZE = size(Images{iImage});
                
                %SET FIGURE SIZE
                %obj.fsxy(1) = 4; %set width
                obj.fsxy(2) = 4; %set height
                %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
                obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
                
                %CREATE FIGURE WITH AXIS
                [h, g] = setfigax1(obj); %create figure with axis with real-size obj.fsxy
                
                %PLOT IMAGE
                imagesc(Images{iImage}); colormap gray;
                
                if iImage >= 2
                    hold on
                    plot([500, 500], [0, 500], 'Color', [0, 166, 214]/255, 'LineWidth', 2)
                end
                
                g.Position = [0, 0, 1, 1]; axis off;
                %EXPORT IMAGE
                realsizeandsave(obj, h, ['MidBrainBand', num2str(iImage)]); %export with savename
            end
            
            
        end
        function showMidBrainBandIntersection(obj)
            
            %LOAD/COMPUTE VALUES
            rectangle = obj.CompleteTemplate.polarMidbrainBand(:, 500);
            rectangleSmooth = obj.CompleteTemplate.polarMidbrainBandWithGaussianWithBlurr(:, 500);
            
            
            %SET FIGURE SIZE
            obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 4; %set height
            %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
            %obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj); %create figure with axis with real-size obj.fsxy
            hold on
            plot(rectangleSmooth, 'Color', [110, 187, 213]/255, 'LineWidth', 2)
            plot(rectangle, 'Color', [0, 166, 214]/255, 'LineWidth', 2)
            
            %SET AXIS
            setfigax2(obj, g);
            set(g, 'XLim', [1, 500]);
            set(g, 'XTick', [1, 500]);
            %set(g, 'XTickLabels', [0.90,{''},0.92,{''},0.94,{''},0.96,{''},0.98,{''},1]);
            %set(g, 'XTickLabels', [0,{'1/2\pi'},{'\pi'},{'3/2\pi'},{'2\pi'}]);
            g.XLabel.String = 'row pixels';
            
            set(g, 'YLim', [0, 1.1]);
            %set(g, 'YScale', 'log');
            set(g, 'YTick', [0, 1]);
            %set(g, 'YTickLabels', [1, 2, 5, 10, {'10^2'}, {''}]);
            %set(g, 'YMinorGrid', 'off');
            g.YLabel.String = 'filter value';
            
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['MidBrainBandIntersection']); %export with savename
        end
        function showMidBrainEdgeContour(obj)
            
            %LOAD/COMPUTE VALUES
            Image = obj.CompleteTemplate.EdgeFilter;
            %SIZE = size(Image);
            
            %SET FIGURE SIZE
            obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 4; %set height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj);
            hold on
            plot(Image(:, :, 1)/255, 'Color', 'red', 'LineWidth', 2)
            plot(Image(:, :, 2)/255, 'Color', 'green', 'LineWidth', 2)
            plot(Image(:, :, 3)/255, 'Color', 'blue', 'LineWidth', 2)
            
            %SET AXIS
            setfigax2(obj, g);
            set(g, 'XLim', [1, 31]);
            set(g, 'XTick', [1, 16, 31]);
            set(g, 'XTickLabels', [-15, 0, 15]);
            %set(g, 'XTickLabels', [0,{'1/2\pi'},{'\pi'},{'3/2\pi'},{'2\pi'}]);
            g.XLabel.String = 'crosssection [pix]';
            
            set(g, 'YLim', [0, 1]);
            %set(g, 'YScale', 'log');
            set(g, 'YTick', [0, 1]);
            %set(g, 'YTickLabels', [1, 2, 5, 10, {'10^2'}, {''}]);
            %set(g, 'YMinorGrid', 'off');
            g.YLabel.String = 'intensity';
            
            %legend('Red channel', 'Green channel', 'Blue channel')
            
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['MidBrainEdgeContour']);
        end
        function showMidBrainEdgeTemplate(obj)
            
            %LOAD/COMPUTE VALUES
            Image = obj.CompleteTemplate.EdgeFilter;
            SIZE = size(Image);
            
            %SET FIGURE SIZE
            %obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 4; %set height
            %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
            obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj);
            
            %PLOT IMAGE
            imagesc(uint8(Image))
            
            g.Position = [0, 0, 1, 1]; axis off;
            
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['MidBrainEdgeTemplate']);
        end
        
        %spot template images
        function showSpotContrast(obj, fn)
            
            %output "obj" is to store Registration parameters in object
            %which can be used for showMethods on the same fish
            
            obj.fsxy = ([10, 8]);
            
            %LOAD/COMPUTE VALUES
            
            spotColorsUnique = single(unique(obj.CompleteTemplate.spotcolorsTT, 'rows')) / 255;
            backrColorsUnique = single(unique(obj.CompleteTemplate.backrcolorsTT, 'rows')) / 255;
            
            %PresetShowRegistration(obj,fn);
            [h, g] = setfigax1(obj); %create figure with axis with real-size obj.fsxy
            
            %PLOT FIGURES
            
            scatter3(spotColorsUnique(:, 1), spotColorsUnique(:, 2), ...
                spotColorsUnique(:, 3), 'Marker', '.', 'MarkerEdgeAlpha', 0.5, ...
                'MarkerEdgeColor', [230, 70, 22]/255);
            hold on
            scatter3(backrColorsUnique(:, 1), backrColorsUnique(:, 2), ...
                backrColorsUnique(:, 3), 'Marker', '.', 'MarkerEdgeAlpha', 0.5, ...
                'MarkerEdgeColor', [0, 166, 214]/255);
            view(g, [-63.5, 18.8]);
            
            setfigax2(obj, g); %preset axis with standard fond and size
            %set(g, 'FontName', 'arial', 'FontSize', 8, 'XGrid', 'on', 'YGrid', 'on');
            %set(g, 'Units', 'centimeters', 'Position', [1.3, 1, obj.fsxy(1) - 2.2, obj.fsxy(2) - 1.2]);
            
            %SET AXIS
            set(g, 'XLim', [0, 1]);
            g.XLabel.String = 'Red';
            set(g, 'YLim', [0, 1]);
            g.YLabel.String = 'Green';
            set(g, 'ZLim', [0, 1]);
            g.ZLabel.String = 'Blue';
            
            realsizeandsave(obj, h, ['SpotTemplate', filesep, 'SpotBackgroundScatterPlot']);
        end
        function showSpotSegmentation(obj)
            %this showfunction is not automatic
            %the fish and spot number in the function TemplateSpot in SNGTemplate has to manually changed
            %TODO update SNGTemplate, than change input to ,fn and sp (fish and spotnumber input)
            
            image{1} = obj.CompleteTemplate.Spot1;
            image{2} = obj.CompleteTemplate.SpotF;
            image{3} = obj.CompleteTemplate.BackF;
            image{4} = obj.CompleteTemplate.Mask2;
            
            for iImage = 1:4
                
                %LOAD/COMPUTE VALUES
                %Image = obj.CompleteTemplate.Template;
                SIZE = size(image{iImage});
                
                %IMAGE ENHANCEMENT
                
                %SET FIGURE SIZE
                %obj.fsxy(1) = 5; %set width
                obj.fsxy(2) = 4; %set height
                %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
                obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
                
                %CREATE FIGURE WITH AXIS
                [h, g] = setfigax1(obj);
                
                %PLOT IMAGE
                imagesc(image{iImage})
                
                g.Position = [0, 0, 1, 1]; axis off;
                
                %EXPORT IMAGE
                realsizeandsave(obj, h, ['SpotTemplate', filesep, 'SpotSegmentation', num2str(iImage)]);
            end
            
        end
        function showSpotHistogram(obj, fn)
            
            %output "obj" is to store Registration parameters in object
            %which can be used for showMethods on the same fish
            
            obj.fsxy = ([10, 8]);
            
            %LOAD/COMPUTE VALUES
            S = double(obj.CompleteTemplate.spotcolorsTT);
            B = double(obj.CompleteTemplate.backrcolorsTT);
            
            bin = 255;
            
            
            Shist = sng_chistcount(S', linspace(0, 255, bin+1)); %create bins
            Bhist = sng_chistcount(B', linspace(0, 255, bin+1)); %
            
            Snorm = Shist / size(S, 1);
            Bnorm = Bhist / size(B, 1);
            
            [minC1, C1, Q1] = sng_threshold2(Snorm, Bnorm); %find threshold
            
            sng_chistplot(Snorm, Bnorm, minC1); set(gcf, 'numbertitle', 'off', 'name', 'RGB')
            
            %% a new way to transform to gray instead of taking the green values
            bm = mean(B)
            sm = mean(S)
            
            
            b = (bm - sm)
            bnorm = b / norm(b)
            
            mx = dot([255, 255, 255], bnorm)
            mn = dot([0, 0, 0], bnorm)
            
            st = double(S) * bnorm'; %transformed spot
            bt = double(B) * bnorm'; %transformed background
            
            st = st * 255 / mx;
            bt = bt * 255 / mx;
            
            SThist = histcounts(st, linspace(0, 255, bin+1));
            BThist = histcounts(bt, linspace(0, 255, bin+1));
            
            STnorm = SThist / size(S, 1);
            BTnorm = BThist / size(B, 1);
            
            [minC2, C2, Q2] = sng_threshold2(STnorm, BTnorm); %find threshold
            
            sng_chistplot(STnorm, BTnorm, minC2); set(gcf, 'numbertitle', 'off', 'name', 'RGB')
            
            
            Overlap = min(STnorm, BTnorm);
            Combination = max(STnorm, BTnorm);
            Jaccard = 1 - sum(Overlap) / sum(Combination)
            
            
            spotColorsUnique = single(unique(obj.CompleteTemplate.spotcolorsTT, 'rows')) / 255;
            backrColorsUnique = single(unique(obj.CompleteTemplate.backrcolorsTT, 'rows')) / 255;
            
            %PresetShowRegistration(obj,fn);
            [h, g] = setfigax1(obj); %create figure with axis with real-size obj.fsxy
            
            %PLOT FIGURES
            
            scatter3(spotColorsUnique(:, 1), spotColorsUnique(:, 2), ...
                spotColorsUnique(:, 3), 'Marker', '.', 'MarkerEdgeAlpha', 0.5, ...
                'MarkerEdgeColor', [230, 70, 22]/255);
            hold on
            scatter3(backrColorsUnique(:, 1), backrColorsUnique(:, 2), ...
                backrColorsUnique(:, 3), 'Marker', '.', 'MarkerEdgeAlpha', 0.5, ...
                'MarkerEdgeColor', [0, 166, 214]/255);
            view(g, [-63.5, 18.8]);
            
            setfigax2(obj, g); %preset axis with standard fond and size
            %set(g, 'FontName', 'arial', 'FontSize', 8, 'XGrid', 'on', 'YGrid', 'on');
            %set(g, 'Units', 'centimeters', 'Position', [1.3, 1, obj.fsxy(1) - 2.2, obj.fsxy(2) - 1.2]);
            
            %SET AXIS
            set(g, 'XLim', [0, 1]);
            g.XLabel.String = 'Red';
            set(g, 'YLim', [0, 1]);
            g.YLabel.String = 'Green';
            set(g, 'ZLim', [0, 1]);
            g.ZLabel.String = 'Blue';
            
            realsizeandsave(obj, h, ['SpotTemplate', filesep, 'SpotBackgroundScatterPlot']);
        end
        
        %brain images
        function showBrainSegmentation(obj, fn)
            %LOAD/COMPUTE VALUES
            Image{1} = obj.BrainSegmentationObject(fn).Isquare;
            Image{2} = obj.BrainSegmentationObject(fn).Ipolar;
            Image{3} = obj.BrainSegmentationObject(fn).ICorr2;
            Image{4} = obj.BrainSegmentationObject(fn).INorm;
            Image{5} = obj.BrainSegmentationObject(fn).INorm2;
            Image{6} = obj.BrainSegmentationObject(fn).Mask;
            Image{7} = obj.BrainSegmentationObject(fn).PolarTransform;
            Image{8} = obj.BrainSegmentationObject(fn).Ibrain;
            for iImage = [1, 2, 4, 5, 7]
                
                SIZE = size(Image{iImage});
                
                %IMAGE ENHANCEMENT
                
                %SET FIGURE SIZE
                %obj.fsxy(1) = 5; %set width
                obj.fsxy(2) = 4; %set height
                %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
                obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
                
                %CREATE FIGURE WITH AXIS
                [h, g] = setfigax1(obj);
                
                %PLOT IMAGE
                imagesc(Image{iImage})
                
                g.Position = [0, 0, 1, 1]; axis off;
                
                %EXPORT IMAGE
                realsizeandsave(obj, h, ['Brain Detection', filesep, 'MidBrainSegmentation_fish', num2str(fn), '_', num2str(iImage)]);
            end
        end
        function showBrainPath(obj, fn)
            %LOAD/COMPUTE VALUES
            Image = obj.BrainSegmentationObject(fn).PolarTransform;
            Path = obj.BrainSegmentationObject(fn).ShortestPath;
            Paths = obj.BrainSegmentationObject(fn).AllShortestPaths;
            
            SIZE = size(Image);
            
            %IMAGE ENHANCEMENT
            
            %SET FIGURE SIZE
            %obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 4; %set height
            %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
            obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj);
            
            %PLOT IMAGE
            imagesc(imcomplement(Image));
            %imagesc((Image));
            
            hold on
            for ipath = 1:numel(Paths)
                plot(Paths{ipath}(:, 1), Paths{ipath}(:, 2), 'Color', [226 / 255, 26 / 255, 26 / 255], 'LineWidth', 1)
            end
            plot(Path(:, 1), Path(:, 2), 'Color', [0, 0, 0], 'LineWidth', 1)
            
            g.Position = [0, 0, 1, 1]; axis off;
            
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['Brain Detection', filesep, 'MidBrainPaths', '_fish', num2str(fn)]);
        end
        function showBrainEdge(obj, fn)
            %LOAD/COMPUTE VALUES
            Image = obj.RegObject(fn).Ialigned;
            %Image = obj.BrainSegmentationObject(fn).Isquare;
            Path = obj.BrainSegmentationObject(fn).BrainEdge;
            %Paths = obj.BrainSegmentationObject(fn).AllShortestPaths;
            
            %hold on;plot(bc{1}(:,2),bc{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
            
            
            SIZE = size(Image);
            
            %IMAGE ENHANCEMENT
            
            %SET FIGURE SIZE
            %obj.fsxy(1) = 5; %set width
            obj.fsxy(2) = 4; %set height
            %obj.fsxy(2) = SIZE(1) / SIZE(2) * obj.fsxy(1); %height is dependent on width
            obj.fsxy(1) = SIZE(2) / SIZE(1) * obj.fsxy(2); %width is dependent on height
            
            %CREATE FIGURE WITH AXIS
            [h, g] = setfigax1(obj);
            
            %PLOT IMAGE
            imagesc(Image);
            %imagesc((Image));
            
            hold on
            plot(Path(:, 2), Path(:, 1), 'Color', [0, 0, 0], 'LineWidth', 1)
            
            g.Position = [0, 0, 1, 1]; axis off;
            
            %EXPORT IMAGE
            realsizeandsave(obj, h, ['Brain Detection', filesep, 'MidBrainEdge', '_fish', num2str(fn)]);
        end
        
        %supporting functions used for the show functions
        function [f1, a1] = setfigax1(obj)
            %fsx = 5;fsy = 10;
            f1 = figure('PaperUnits', 'centimeters', 'Color', [1, 1, 1]);
            a1 = gca;
            
            obj.sng_figcm(obj.fsxy(1), obj.fsxy(2));
        end
        function setfigax2(obj, a1)
            set(a1, 'FontName', 'arial', 'FontSize', 8, 'XGrid', 'on', 'YGrid', 'on');
            set(a1, 'Units', 'centimeters', 'Position', [1.2, 1.2, obj.fsxy(1) - 1.7, obj.fsxy(2) - 1.7]);
            %set(axishandle, 'YLim', [-0.02, 1.02]);
        end
        function realsizeandsave(obj, f1, name)
            if ~exist('name', 'var')
                name = 'noname';
            end
            
            [destinationFolder, name] = fileparts([obj.SavePath, filesep, 'Images', filesep, name])
            
            if ~exist(destinationFolder, 'dir')
                mkdir(destinationFolder)
            end
            
            if obj.exportit
                export_fig(f1, [destinationFolder, filesep, name], '-png', '-r600', '-nocrop');
            end
            
            %sets the right scaling
            User1 = getenv('USER');
            User2 = getenv('username');
            OS = getenv('OS');
            if strcmp(User1, 'samuelgeurts') && isempty(User2) && isempty(OS)
                ScaledFigure.calibrateDisplay(113.6); %113.6 for macbook pro screen, 96 inch for erasmus screen
            elseif strcmp(OS, 'Windows_NT') && strcmp(User2, '260018')
                ScaledFigure.calibrateDisplay(96); %113.6 for macbook pro screen, 96 inch for erasmus screen
            else
            end
            
            ScaledFigure(f1, 'reuse');
            set(f1, 'Units', 'Centimeters');
            set(f1, 'Position', (get(f1, 'Position') + [obj.fsxy(1), 0, 0, 0]));
        end
    end
    
    
    %{
     methods(Hidden = true)
 
function [figurehandle, axishandle] = setfigax1(sobj)
             %fsx = 5;fsy = 10;
             figurehandle = figure('PaperUnits', 'centimeters', 'Color', [1, 1, 1]);
             axishandle = gca;
             sng_figcm(sobj.fsx, sobj.fsy);
         end
function setfigax2(sobj, axishandle)
             set(axishandle, 'FontName', 'arial', 'FontSize', 8, 'XGrid', 'on', 'YGrid', 'on');
             set(axishandle, 'Units', 'centimeters', 'Position', [1.2, 1.2, sobj.fsx - 1.7, sobj.fsy - 1.7]);
             %set(axishandle, 'YLim', [-0.02, 1.02]);
 
         end
function realsizeandsave(sobj, figurehandle)
             if sobj.exportit
                 export_fig(figurehandle, [sobj.SavePath, num2str(figurehandle.Number)], '-png', '-r600', '-nocrop');
             end
 
             %sets the right scaling
             User1 = getenv('USER');
             User2 = getenv('username');
             OS = getenv('OS');
             if strcmp(User1,'samuelgeurts') && isempty(User2) && isempty(OS)
                 ScaledFigure.calibrateDisplay(113.6); %113.6 for macbook pro screen, 96 inch for erasmus screen
             elseif  strcmp(OS, 'Windows_NT') && strcmp(User2, '260018')
                 ScaledFigure.calibrateDisplay(96); %113.6 for macbook pro screen, 96 inch for erasmus screen
             else
             end
 
 
             ScaledFigure(figurehandle, 'reuse');
             set(figurehandle, 'Units', 'Centimeters');
             set(figurehandle, 'Position', (get(figurehandle, 'Position') + [sobj.fsx, 0, 0, 0]));
         end
 
     end
    %}
end