function ShowFishDiscrimination(obj,exportit,fsxy)

if exist('fsxy','var')
    fsx = fsxy(1);
    fsy = fsxy(2);
else
    fsx = 6; 
    fsy = 5;
end



if ~exist('exportit','var')
    exportit = false;
end



load([obj.SavePath, '/', obj.InfoName, '.mat'], 'ImageInfo')

a = [ImageInfo.corcoef];
b = [ImageInfo.anglewarp];
c = [ImageInfo.moduluswarp];
d = [ImageInfo.CorNextStack];
e = reshape([ImageInfo.warp],2,size(ImageInfo,1))'


%histogram shows correlation coefficientadjacent images of the same fish
[h2, g2] = setfigax1;
histogram(a(d == 0),0.98:0.001:1);
hold on
histogram(a(d == 1),0.98:0.001:1);
setfigax2(g2);
set(g2, 'XLim', [0.98, 1]);
g2.XLabel.String = 'correlation coefficient';
g2.YLabel.String = 'counts';
realsizeandsave(h2);

%histogram shows correlation coefficient adjacent images containing different fish
[h1, g1] = setfigax1;
histogram(a(d == 0),0:0.05:1);
hold on
hh = histogram(a(d == 1),0:0.05:1);
set(hh,'FaceColor',[0.8500    0.3250    0.0980]);
set(gca,'YScale','log')
set(gca,'YLim',[0.6,301])
set(gca,'YTick',[1,2,5,10,100])
set(gca,'YTickLabels',[1,2,5,10,{'10^2'}])
setfigax2(g1)
g1.XLabel.String = 'correlation coefficient';
g1.YLabel.String = 'counts';
realsizeandsave(h1)

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
%}
%histogram of angular shift zoomed out

%{
[h3, g3] = setfigax1;
    h=histogram(b(d == 1),[-pi:pi/20:pi])
    set(h,'FaceColor',[0.8500    0.3250    0.0980])  
    h=histogram(b(d == 0),[-pi:pi/20:pi])
%}


%histogram of translation zoomed out
[h4, g4] = setfigax1;
scatter(e((d==0),1),e((d==0),2))
hold on
scatter(e((d==1),1),e((d==1),2)) 
xlim([-20,20])
ylim([-20,20])
setfigax2(g4)
g4.XLabel.String = 'translation x-axis [pix]';
g4.YLabel.String = 'translation y-axis [pix]';
realsizeandsave(h4)

%histogram of translation zoomed in
[h4, g4] = setfigax1;
scatter(e((d==0),1),e((d==0),2))
hold on
scatter(e((d==1),1),e((d==1),2)) 
xlim([-1.5,1.5])
ylim([-1.5,1.5])
setfigax2(g4)
g4.XLabel.String = 'translation x-axis [pix]';
g4.YLabel.String = 'translation y-axis [pix]';
realsizeandsave(h4)






function [figurehandle, axishandle] = setfigax1
    %fsx = 5;fsy = 10;
    figurehandle = figure('PaperUnits', 'centimeters', 'Color', [1, 1, 1]);
    axishandle = gca;
    sng_figcm(fsx, fsy);
end
function setfigax2(axishandle)
    set(axishandle, 'FontName', 'arial', 'FontSize', 8, 'XGrid', 'on', 'YGrid', 'on');
    set(axishandle, 'Units', 'centimeters', 'Position', [1.2, 1.2, fsx - 1.7, fsy - 1.7]);
    %set(axishandle, 'YLim', [-0.02, 1.02]);
end
function realsizeandsave(figurehandle)
    if exportit
        export_fig(figurehandle, [obj.SavePath,filesep, 'ImageSorting', num2str(figurehandle.Number)], '-png', '-r600', '-nocrop');
    end
    %ScaledFigure.calibrateDisplay(113.6); %113.6 for macbook pro screen, 96 inch for erasmus screen
    ScaledFigure(figurehandle, 'reuse');
    set(figurehandle, 'Units', 'Centimeters');
    set(figurehandle, 'Position', (get(figurehandle, 'Position') + [fsx, 0, 0, 0]));
end








end