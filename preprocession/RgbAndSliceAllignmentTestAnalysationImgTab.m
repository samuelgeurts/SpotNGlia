%open result from computations on artificial images and shows images

%Variables


fsx = 8.5; %widht in cm 
fsy = 8.5; %height in cm
Color1 = [0 166/255 214/255]; %scatter color




PathsComplete('or')


load([OriginalPath,'/RgbAlignmentTest.mat'],'RgbAlTestTable')
load([OriginalPath,'/SliceAlignmentTest.mat'],'SliceAlTestTable')

%color shift in terms of angle and modulus
angi = atan2(RgbAlTestTable.InitialWarp(:,1),RgbAlTestTable.InitialWarp(:,2))
modi = sqrt(RgbAlTestTable.InitialWarp(:,1).^2+RgbAlTestTable.InitialWarp(:,2).^2)
ango = atan2(RgbAlTestTable.ComputedWarp(:,1),RgbAlTestTable.ComputedWarp(:,2))
modo = sqrt(RgbAlTestTable.ComputedWarp(:,1).^2+RgbAlTestTable.ComputedWarp(:,2).^2)
angof = atan2(RgbAlTestTable.ComputedWarpFilt(:,1),RgbAlTestTable.ComputedWarpFilt(:,2))
modof = sqrt(RgbAlTestTable.ComputedWarpFilt(:,1).^2+RgbAlTestTable.ComputedWarpFilt(:,2).^2)

%%
%color allignment without filter
difference1 = abs(modi - modo)
mean(difference1)
std(difference1)
var(difference1)
figure;histogram(difference1)

difference2 = abs(modi - modof)
mean(difference2)
std(difference2)
var(difference2)
figure;histogram(difference2)
sum(difference2>=1)/numel(difference2)

%{
figure;scatter(modi,difference2)
figure;scatter(modi,difference1)
%}

%{
difference = abs(angi - ango)
mean(difference)
std(difference)

[histo bin] = hist(difference,linspace(0,0.2,100))
figure;histogram(difference,[linspace(0,0.1,25),6])
xlim([0,0.15])
set(gca,'XTickLabel',{0,0.05,0.1,6})
%}

% scatter plot of original, and computed calculations
%{
figure;scatter(RgbAlTestTable.InitialWarp(:,1),RgbAlTestTable.InitialWarp(:,2))
hold on;scatter(RgbAlTestTable.ComputedWarp(:,1),RgbAlTestTable.ComputedWarp(:,2))
hold on;scatter(RgbAlTestTable.ComputedWarpFilt(:,1),RgbAlTestTable.ComputedWarpFilt(:,2))
line([-10,4],[0,0])
line([0,0],[-4,10])
legend
%}

% for error scatterplot

i = RgbAlTestTable.InitialWarp;
o = RgbAlTestTable.ComputedWarp;
f = RgbAlTestTable.ComputedWarpFilt;

nofilt = RgbAlTestTable.ComputedWarp-RgbAlTestTable.InitialWarp;
filt = RgbAlTestTable.ComputedWarpFilt - RgbAlTestTable.InitialWarp;



%% slice allignment
angsi = atan2(SliceAlTestTable.IniWarp(:,1),SliceAlTestTable.IniWarp(:,2))
modsi = sqrt(SliceAlTestTable.IniWarp(:,1).^2+SliceAlTestTable.IniWarp(:,2).^2)
angso = atan2(SliceAlTestTable.CalWarp(:,1),SliceAlTestTable.CalWarp(:,2))
modso = sqrt(SliceAlTestTable.CalWarp(:,1).^2+SliceAlTestTable.CalWarp(:,2).^2)
mean(abs(modsi-modso))
std(abs(modsi-modso))
hist(abs(modsi-modso))

%figure;scatter(modsi,abs(modsi-modso))
SliceError = SliceAlTestTable.IniWarp - SliceAlTestTable.CalWarp;

%%


h1 = figure('PaperUnits','centimeters','Color',[1 1 1]);
sng_figcm(fsx,fsy)

s1 = scatter(filt(:,1),filt(:,2),'MarkerEdgeColor',Color1,'SizeData',10,'LineWidth',0.5,'Marker','o')

xlabel('x-displacement error [pix]','FontSize',8,'FontName','arial');
ylabel('y-displacement error [pix]','FontSize',8,'FontName','arial');

set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
set(gca,'XLim',[-7,7])
set(gca,'YLim',[-7,7])
set(gca,'Units','centimeters','Position',[1.2 1.2 fsx-1.7 fsy-1.7])
set(gca,'XTick',[-6,-4,-2,0,2,4,6])
set(gca,'YTick',[-6,-4,-2,0,2,4,6])

export_fig(h1 ,['/Users/samuelgeurts/Desktop/','DisplacementErrorFilt'], '-png', '-r600', '-nocrop');

%%

%figure('InvertHardcopy','off','PaperUnits','centimeters','Color',[1 1 1]);
h2 = figure('PaperUnits','centimeters','Color',[1 1 1]);
sng_figcm(fsx,fsy)
scatter(nofilt(:,1),nofilt(:,2),'MarkerEdgeColor',Color1,'SizeData',10,'LineWidth',0.5,'Marker','o')

xlabel('x-displacement error [pix]','FontSize',8,'FontName','arial');
ylabel('y-displacement error [pix]','FontSize',8,'FontName','arial');

set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
set(gca,'XLim',[-0.25,0.25])
set(gca,'YLim',[-0.25,0.25])
set(gca,'Units','centimeters','Position',[1.3 1.2 fsx-1.7 fsy-1.7])
set(gca,'XTick',[-0.2,-0.1,0,0.1,0.2])
set(gca,'YTick',[-0.2,-0.1,0,0.1,0.2])

export_fig(h2 ,['/Users/samuelgeurts/Desktop/','DisplacementErrorNoFilt'], '-png', '-r600', '-nocrop');

%%

h3 = figure('PaperUnits','centimeters','Color',[1 1 1]);
sng_figcm(fsx,fsy)
scatter(SliceError(:,1),SliceError(:,2),'MarkerEdgeColor',Color1,'SizeData',10,'LineWidth',0.5,'Marker','o')

xlabel('x-displacement error [pix]','FontSize',8,'FontName','arial');
ylabel('y-displacement error [pix]','FontSize',8,'FontName','arial');

set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
set(gca,'XLim',[-0.12,0.12])
set(gca,'YLim',[-0.12,0.12])
set(gca,'Units','centimeters','Position',[1.4 1.2 fsx-1.7 fsy-1.7])
set(gca,'XTick',[-0.1,-0.05,0,0.05,0.1])
set(gca,'YTick',[-0.1,-0.05,0,0.05,0.1])

export_fig(h3 ,['/Users/samuelgeurts/Desktop/','SliceDisplacementError'], '-png', '-r600', '-nocrop');
