function [wp, a] = DefineSpotThresholds(obj, fishnumbers)
%compute spot thresholds based on fals en positives

INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo', 'BrainSegmentationInfo', 'SpotParameters');
TEMPLATE = load([obj.SourcePath, '/', 'Template3dpf.mat'], 'ref_temp', 'SVAP_index', 'SpotVectorArrayProbability');

if ~exist('fishnumbers', 'var') || isempty(fishnumbers)
    fishnumbers = 1:numel(obj.StackInfo);
elseif max(fishnumbers) > numel(obj.StackInfo)
    error('at least one fish does not exist in StackInfo')
end

            
if sum(~cellfun('isempty',{obj.Annotations.Spots})) ~= numel(SpotParameters)
    warning('not all fishes are annotated by hand')
end
fishnumbers = 1:numel(SpotParameters);           
%removes fishnumbers from the list if not all annotations are known
fishnumbers = fishnumbers(ismember(fishnumbers,find(~cellfun('isempty',{obj.Annotations.Spots}))));

nfishes = numel(fishnumbers);





%[MidBrain, AlignedFish] = loadFishAndBrain(obj, fishnumbers, INFO,TEMPLATE, 'Annotation', 1);
%[MidBrain, AlignedFish2] = loadFishAndBrain(obj, fishnumbers, INFO,TEMPLATE, 'Annotation', 0);
[MidBrain, ~] = loadFishAndBrain(obj, fishnumbers, INFO, TEMPLATE, 'Correction', 0);
%[MidBrain, ~] = loadFishAndBrain(obj, fishnumbers, INFO,TEMPLATE, 'Annotation', 0);


%if ~exist('INFO','var') || ~isfield(INFO,'SpotParameters')
%    INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters');
%end
%INFO.SpotParameters

if ~isfield(INFO, 'BrainSegmentationInfo')
    INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters');
end
SpotParameters = SpotParametersBrainSelect(INFO.SpotParameters, MidBrain);

%preallocation
SpotparC = cell(1, nfishes);
SpotparI = cell(1, nfishes);
sizeC = zeros(1, nfishes);
sizeI = zeros(1, nfishes);


%split up correct and incorrect found spots
for k1 = 1:nfishes
    fn = fishnumbers(k1);
    if ~isempty(MidBrain{fn})
        Spotpar = SpotParameters{fn};
        SpotCom = reshape([SpotParameters{fn}.Centroid], 2, numel(SpotParameters{fn}))';
        atemp = ismember(SpotCom, obj.SpotInfo(fn).CorrectSpots);
        SpotparC{k1} = Spotpar(atemp(:, 1));
        SpotparI{k1} = Spotpar(~atemp(:, 1));
        sizeC(k1) = sum(atemp(:, 1));
        sizeI(k1) = numel(Spotpar) - sizeC(k1);
    end
end


%start and end indices
endC = cumsum(sizeC);
startC = [1, endC(1:end-1) + 1];
endI = cumsum(sizeI);
startI = [1, endI(1:end-1) + 1];

%preallocation
AreaC = zeros(sum(sizeC), 1);
ColorC = zeros(sum(sizeC), 1);
SpotMeanC = zeros(sum(sizeC), 3);
SpotMeanhsvC = zeros(sum(sizeC), 3);
AzimC = zeros(sum(sizeC), 1);
ElevC = zeros(sum(sizeC), 1);
RC = zeros(sum(sizeC), 1);
AreaI = zeros(sum(sizeI), 1);
ColorI = zeros(sum(sizeI), 1);
SpotMeanI = zeros(sum(sizeI), 3);
SpotMeanhsvI = zeros(sum(sizeI), 3);
AzimI = zeros(sum(sizeI), 1);
ElevI = zeros(sum(sizeI), 1);
RI = zeros(sum(sizeI), 1);

for k1 = 1:nfishes
    subC = SpotparC{k1};
    subI = SpotparI{k1};
    
    indc = startC(k1):endC(k1);
    indi = startI(k1):endI(k1);
    
    if ~isempty(indc)
        AreaC(indc, 1) = [subC.Area];
        ColorC(indc, 1) = [subC.ColorProbability];
        SpotMeanC(indc, 1:3) = reshape([subC.ColorMean], 3, numel(subC))';
        SpotMeanhsvC = rgb2hsv(SpotMeanC/255);
        AzimC(indc, 1) = [subC.azimuth];
        ElevC(indc, 1) = [subC.elevation];
        RC(indc, 1) = [subC.r];
 
        %...
    end
    
    if ~isempty(indi)
        AreaI(indi, 1) = [subI.Area];
        ColorI(indi, 1) = [subI.ColorProbability];
        SpotMeanI(indi, 1:3) = reshape([subI.ColorMean], 3, numel(subI))';
        SpotMeanhsvI = rgb2hsv(SpotMeanI/255);
        AzimI(indi, 1) = [subI.azimuth];
        ElevI(indi, 1) = [subI.elevation];
        RI(indi, 1) = [subI.r];
        %...
    end
    
end








%{
     figure;imagesc(AlignedFish2{1})
     scI = scatter(0, 0, 400, ...
     'LineWidth', 2, ...
     'MarkerEdgeColor', 1/255*[75, 75, 255]);
     scI.XData = SpotCom(~atemp(:, 1),1);
     scI.YData = SpotCom(~atemp(:, 1),2);
     %)
 
     x = linspace(0, 1, 1000); %histogram bins
     hazimc = histcounts(SpotMeanhsvC(:,1), x); %correct spot area
     hazimi = histcounts(SpotMeanhsvI(:,1), x); %incorrect spot area
     sng_DistDifHistPlot((x(1:end-1)+((x(2)-x(1))/2)),hazimc,hazimi);
 
     hold on;line([MinSpotSize MinSpotSize],get(gca,'Ylim'),'color',[0 0 0]);
     hold on;line([MaxSpotSize MaxSpotSize],get(gca,'Ylim'),'color',[0 0 0]);
 
 end
%}
        
DoubleHist(SpotMeanC(:,1), SpotMeanI(:,1),'SpotMeanI(1)');set(gca,'XLim',[0 255]);
DoubleHist(SpotMeanC(:,2), SpotMeanI(:,2),'SpotMeanI(2)');set(gca,'XLim',[0 255]);
DoubleHist(SpotMeanC(:,3), SpotMeanI(:,3),'SpotMeanI(3)');set(gca,'XLim',[0 255]);
DoubleHist(SpotMeanhsvC(:,1), SpotMeanhsvI(:,1),'SpotMeanhsv(1)');set(gca,'XLim',[0 1]);
DoubleHist(SpotMeanhsvC(:,2), SpotMeanhsvI(:,2),'SpotMeanhsv(2)');set(gca,'XLim',[0 1]);
%DoubleHist(SpotMeanhsvC(:,3), SpotMeanhsvI(:,3),'SpotMeanhsv(3)');set(gca,'XLim',[0 1]);
DoubleHist(AzimC, AzimI,'Azim');set(gca,'XLim',[-pi pi]);
DoubleHist(ElevC, ElevI,'Elev');set(gca,'XLim',[-pi/2 pi/2]);
DoubleHist(RC, RI,'R');set(gca,'XLim',[0 150]);
DoubleHist(ColorC, ColorI,'ColorProb');set(gca,'XLim',[0 2]);
DoubleHist(AreaC, AreaI,'Area');set(gca,'XLim',[0 500]);




weightC = endC(end)/(endC(end) + endI(end));
weightI = endI(end)/(endC(end) + endI(end));

pdC = fitdist(SpotMeanC(:,2),'Normal')
pdI = fitdist(SpotMeanI(:,2),'Normal')

x = 0:1:255
yC = pdf(pdC,x) * weightC
yI = pdf(pdI,x) * weightI

plot(x,yC,x,yI)






%using a machine learning technique would be an improvement
%DS = [[AreaC;AreaI],[ColorC;ColorI]]
DS = [[SpotMeanhsvC(:, 1); SpotMeanhsvI(:, 1)], ...
    [SpotMeanhsvC(:, 2); SpotMeanhsvI(:, 2)], ...
    [AreaC; AreaI], ...
    [AzimC; AzimI], ...
    [ElevC; ElevI], ...
    [RC; RI]];


lab = [ones(numel(AreaC), 1); 2 * ones(numel(AreaI), 1)];
a = prdataset(DS(:, 1:6), lab);
a = a * scalem(a, 'variance'); %scaling
wp = knnc(a);
err = a * wp * testc %error of testset c on trained classifier wb

end


%scatterd(a,2)
%plotc(wp,1)

