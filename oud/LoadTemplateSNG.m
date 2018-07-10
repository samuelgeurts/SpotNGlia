%it seems that this functions is replaced by MakeTemplateMat en no longer neede to be used


function [CompleteTemplate] = LoadTemplateSNG(TemplatePath3dpf)


%version LoadTemplateLink2
%    add input PAth
%version LoadTemplateLink3
%    Template information based on newest 50 fish 3dpf stack


%% Load template and template Parameters

%TODO:  auto determine fish age
%TODO:  make function which generates template and put all relevant data into one variable
%TODO:  generate template by alligning to complete fish image
%TODO:  find dependency between eyes variance and location brain
%TODO:  exclude weird fishes from template (highly axial rotated)
%       to acchieve a better locate quess in the polar transform based on 
%       the distancemap


%% Mean Fish Registration Template
%based on older template, can be updated with 50 fish batch but lead
%probably not to improvement. Source is not sure, TemplateGeneration?
CompleteTemplate.Template = imread([TemplatePath3dpf,'/template.tif']);

CompleteTemplate.Size = size(CompleteTemplate.Template);
CompleteTemplate.ref_temp = imref2d(CompleteTemplate.Size);


%Eye region, annotated by Anouk H
%annotated on the older template
temp = sng_html2roi2([TemplatePath3dpf,'/edof_eyes.html']);
CompleteTemplate.EyeRegion1 = temp{1}{2};
CompleteTemplate.EyeRegion2 = temp{1}{1};

    eye1 = CompleteTemplate.EyeRegion2;
    eye2 = CompleteTemplate.EyeRegion1;

CompleteTemplate.EyeMask1 = poly2mask(eye1(:,1),eye1(:,2),CompleteTemplate.Size(1),CompleteTemplate.Size(2));
CompleteTemplate.EyeMask2 = poly2mask(eye2(:,1),eye2(:,2),CompleteTemplate.Size(1),CompleteTemplate.Size(2));

    [x1c,y1c] = sng_CenterOfMass(CompleteTemplate.EyeMask1);
    [x2c,y2c] = sng_CenterOfMass(CompleteTemplate.EyeMask2);

CompleteTemplate.EyeCenter1 = [x1c,y1c];
CompleteTemplate.EyeCenter2 = [x2c,y2c];

%{
figure;imagesc(eyem1+eyem2)
%}

%% Brain 
%generated with the function TemplateBrainParameters

%Midbrain info (poly,area,centroid,box) of individual fishes
load([TemplatePath3dpf,'/MidbrainInfo.mat']);
CompleteTemplate.MidbrainInfo = MidbrainInfo;

%Forbrain info (poly,area,centroid,box) of individual fishes
load([TemplatePath3dpf,'/ForbrainInfo.mat']);
CompleteTemplate.ForbrainInfo = ForbrainInfo;

%MeanMidbrain 
temp= load([TemplatePath3dpf,'/MeanMidbrain.mat']);
CompleteTemplate.MeanMidBrain = temp.MeanMidBrain;
CompleteTemplate.CenterMidBrain = temp.Centroidm;
CompleteTemplate.MidBrainDistanceMap = temp.distancemapm;
CompleteTemplate.BandMidBrain = temp.bandm;
%Meanforbrain makes not much sense as it variates largely
%it must be dependt on the midbrean and eyes

%{
figure;image(CompleteTemplate.Template)
hold on
plot(CompleteTemplate.MeanMidBrain(:,2),CompleteTemplate.MeanMidBrain(:,1))
figure;imshow(CompleteTemplate.MidBrainDistanceMap)
%}

%% load old brain info
%I think it is not used already so could be deleted if so
%containes average informatie about brain location and std.
%{
load([TemplatePath3dpf,'/templates/TemplateInfo.mat'])
%midbrain information
load([TemplatePath,'/templates/TemplateInfo.mat'])
CompleteTemplate.MidbrainInfo = TemplateInf.MidbrainRegionInfo; %averaged over all fishes
clear TemplateInf
%} 


%% Brain Edge Template
% has te be removed later on
%a more normal function would be better
temp = load([TemplatePath3dpf,'/EdgeTemplate.mat']);  %polar transforms and edge template
CompleteTemplate.EdgeTemplate1 = temp.EdgeTemplate;


%% Spot Parameters
% generated with SpotTemplate2
load([TemplatePath3dpf,'/SpotTemplate.mat'],'SpotTemplateVar');
CompleteTemplate.SpotContrastVector = SpotTemplateVar.SpotContrastVector;
CompleteTemplate.SpotVectorArrayProbability = SpotTemplateVar.SpotVectorArrayProbability;
CompleteTemplate.SVAP_index = SpotTemplateVar.SVAP_index;
%clear SpotTemplateVar






end