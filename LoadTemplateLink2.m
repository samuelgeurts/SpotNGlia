function [CompleteTemplate] = LoadTemplateLink2(age,TemplatePath,Basepath)


%version LoadTemplateLink2
%    add input PAth

%% Load template and template Parameters

%TODO:  auto determine fish age
%TODO:  make function which generates template and put all relevant data 
%       into one variable
%TODO:  generate template by alligning to complete fish image
%TODO:  find dependency between eyes variance and location brain
%TODO:  exclude weird fishes from template (highly axial rotated)
%       to acchieve a better locate quess in the polar transform based on 
%       the distancemap

%basepath = '/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Data Nynke/Template brain/3 dpf';

CompleteTemplate.Template = imread([TemplatePath,'/templates/fish template/edof.tif']);
%eye segmented on templateTe
CompleteTemplate.EyeRegion = sng_html2roi2([TemplatePath,'/templates/fish template/edof_eyes.html']);
%midbrain information
load([TemplatePath,'/templates/TemplateInfo.mat'])
CompleteTemplate.MidbrainInfo = TemplateInf.MidbrainRegionInfo; %averaged over all fishes
clear TemplateInf
%midbrain average
    %this one is generated for the snake, which contains the complete brain
    %load([basepath,'templates/fish template/skeletontemplatecoords.mat']) 
    %CompleteTemplate.AverageMidBrain = skeletontemplatecoords
load([TemplatePath,'/roi/html 2brain/skeleton.mat']) 
[cx,cy] = sng_CenterOfMass(skeleton);
[X,Y] = find(skeleton);
XY = sng_OrderContourCoordinates([X,Y]);
CompleteTemplate.MeanMidBrain = XY;
CompleteTemplate.CenterMidBrain = [cx,cy];
%clear skeleton

%adjust distancemap to template size
%later on do this in a new templategeneration function which contains
%the function GenerateSkeletonDistanceMap
sct = size(CompleteTemplate.Template);
cr = min([sct(1:2);size(distancemap)]);
dism = zeros(sct(1:2));
dism(1:cr(1),1:(cr(2))) = distancemap(1:cr(1),1:(cr(2)));

CompleteTemplate.MidBrainDistanceMap = dism;
%clear distancemap dism

%{
figure;imshow(CompleteTemplate.Template)
hold on
plot(XY(:,2),XY(:,1))
figure;imshow(CompleteTemplate.MidBrainDistanceMap)
%}
%midbrain edge template

load([TemplatePath,'/templates/EdgeTemplate.mat']);  %polar transforms and edge template
CompleteTemplate.EdgeTemplate1 = EdgeTemplate;
%clear EdgeTemplate


CompleteTemplate.Size = sct;

load([TemplatePath,'/templates/EdgeTemplate.mat']);  %polar transforms and edge template
%????


%% Template information based on newest 50fish stack

Basepath = '/Users/samuelgeurts/Dropbox/20170327_3dpf_test';
load([Basepath,'/Template Spot/SpotTemplate.mat'],'SpotTemplateVar');
CompleteTemplate.SpotContrastVector = SpotTemplateVar.SpotContrastVector;
CompleteTemplate.SpotVectorArrayProbability = SpotTemplateVar.SpotVectorArrayProbability;
CompleteTemplate.SVAP_index = SpotTemplateVar.SVAP_index;
%clear SpotTemplateVar






end