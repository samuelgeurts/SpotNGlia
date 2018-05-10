classdef Show < SpotNGlia
    
    properties
        obj = [];
        fsx = 6;
        fsy = 5;
        exportit = false;
        subject = []
        %SavePath = []
    end
    
    methods
        
        function sobj = Show(obj, subject, exportit, fsxy)
            sobj.obj = obj;
            sobj.SavePath = [obj.SavePath, filesep, 'ImageSorting'];
            
            if exist('subject','var')
                sobj.subject = subject;
            else
                sobj.subject = [];
            end
            
            if exist('fsxy', 'var')
                sobj.fsx = fsxy(1);
                sobj.fsy = fsxy(2);
            end
            
            if exist('exportit', 'var')
                sobj.exportit = exportit;
            end
            
            %ShowFishDiscrimination1 - histogram shows correlation coefficient adjacent images
            %ShowFishDiscrimination2 - scatterplot of translation
            
            if ismember(1,sobj.subject); sobj.ShowFishDiscrimination1;end
            if ismember(2,sobj.subject); sobj.ShowFishDiscrimination2;end
            if ismember(3,sobj.subject); sobj.ShowFishDiscrimination3;end
            if ismember(4,sobj.subject); sobj.ShowFishDiscrimination4;end

        end
        function ShowFishDiscrimination1(sobj)
            
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
        function ShowFishDiscrimination2(sobj)
            a = [sobj.obj.ImageInfo.corcoef];
            d = [sobj.obj.ImageInfo.CorNextStack];
            
            %histogram shows correlation coefficient adjacent images containing different fish
            [h1, g1] = sobj.setfigax1;
            histogram(a(d == 0), 0:0.05:1);
            hold on
            hh = histogram(a(d == 1), 0:0.05:1);
            set(hh, 'FaceColor', [0.8500, 0.3250, 0.0980]);
            set(g1, 'YScale', 'log')
            set(g1, 'YLim', [0.5, 300])
            set(g1, 'XLim', [0,1])
            set(g1, 'YTick', [1, 2, 5, 10 ,100 200])
            set(g1, 'YTickLabels', [1, 2, 5, 10, {'10^2'},{''}])
            set(g1, 'YMinorGrid','off')

            set(g1, 'XTick', 0:0.2:1)
            
            
            sobj.setfigax2(g1)
            g1.XLabel.String = 'correlation coefficient';
            g1.YLabel.String = 'counts';
            sobj.realsizeandsave(h1)
            
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
        function ShowFishDiscrimination3(sobj)
            %scatterplot of translation zoomed out
            
            d = [sobj.obj.ImageInfo.CorNextStack];
            e = reshape([sobj.obj.ImageInfo.warp], 2, size(sobj.obj.ImageInfo, 1))';
            
            [h4, g4] = sobj.setfigax1;
            scatter(e((d == 0), 1), e((d == 0), 2))
            hold on
            scatter(e((d == 1), 1), e((d == 1), 2))
            xlim([-20, 20])
            ylim([-20, 20])
            set(g4, 'XTick', -20:10:20)
            set(g4, 'YTick', -20:10:20)
            
            sobj.setfigax2(g4)
            %make axis equal but leftdown axis intersection at same position
            pos = get(gca,'Position');
            set(gca,'Position',[pos(1:2),min(pos(3:4)),min(pos(3:4))])
        
            g4.XLabel.String = 'translation x-axis [pix]';
            g4.YLabel.String = 'translation y-axis [pix]';
            sobj.realsizeandsave(h4)
        end
        function ShowFishDiscrimination4(sobj)
                        %scatterplot of translation zoomed in

            b = [sobj.obj.ImageInfo.moduluswarp];
            c = [sobj.obj.ImageInfo.anglewarp];
            d = [sobj.obj.ImageInfo.CorNextStack];
            e = reshape([sobj.obj.ImageInfo.warp], 2, size(sobj.obj.ImageInfo, 1))';
            
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

                        

            [h4, g4] = sobj.setfigax1;
            scatter(e((d == 0), 1), e((d == 0), 2))
            hold on
            scatter(e((d == 1), 1), e((d == 1), 2))
            
            xlim([-1.6, 1.6])
            ylim([-1.6, 1.6])
            sobj.setfigax2(g4)
            
            %make axis equal but leftdown axis intersection at same position
            pos = get(gca,'Position');
            set(gca,'Position',[pos(1:2),min(pos(3:4)),min(pos(3:4))])
            
            
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
            maxdeg = round(360*c3max/(2*pi),1);
            mindeg = round(360*c3min/(2*pi),1);          
            text(0,1.3,['r = ',num2str(round(r,1))],'FontName', 'arial', 'FontSize', 8,'HorizontalAlignment','Center');
            text(1.65,-0.32,[num2str(maxdeg),'\circ'],'FontName', 'arial', 'FontSize', 8);
            text(1.65,-0.95,[num2str(mindeg),'\circ'],'FontName', 'arial', 'FontSize', 8);
            
            sobj.realsizeandsave(h4);

        end       
        function ShowBackgroundRemoval(sobj,fn)
            
            CombinedFish = imread([obj.SavePath, '/', 'CombinedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
            [AlignedFishTemp, RegistrationInfo{k1, 1}] = AlignmentSNG(CombinedFish, obj.CompleteTemplate, obj.ZFParameters);



            
            
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
    
end