%Combines informatie from Source folder with different informatie about the template ands stores it in a single mat file also
%in the source folder.



SourcePath = uigetdir;

%{
SourcePath = obj.SourcePath
%}
TemplatePath3dpf = [SourcePath, '/', 'Template 3 dpf'];

%% Mean Fish Registration Template
%based on older template, can be updated with 50 fish batch but lead
%probably not to improvement. Source is not sure, TemplateGeneration?

Template = imread([TemplatePath3dpf,'/template.tif']);
Size = size(Template);
ref_temp = imref2d(Size);


%Eye region, annotated by Anouk H
%annotated on the older template
temp = sng_html2roi2([TemplatePath3dpf,'/edof_eyes.html']);
EyeRegion1 = temp{1}{2};
EyeRegion2 = temp{1}{1};

    eye1 = EyeRegion2;
    eye2 = EyeRegion1;

EyeMask1 = poly2mask(eye1(:,1),eye1(:,2),Size(1),Size(2));
EyeMask2 = poly2mask(eye2(:,1),eye2(:,2),Size(1),Size(2));
    clear eye1 eye2

    [x1c,y1c] = sng_CenterOfMass(EyeMask1);
    [x2c,y2c] = sng_CenterOfMass(EyeMask2);

EyeCenter1 = [x1c,y1c];
EyeCenter2 = [x2c,y2c];


%{
figure;imagesc(eyem1+eyem2)
%}

%% Brain 
%generated with the function TemplateBrainParameters

%Midbrain info (poly,area,centroid,box) of individual fishes
load([TemplatePath3dpf,'/MidbrainInfo.mat']);
MidbrainInfo = MidbrainInfo;

%Forbrain info (poly,area,centroid,box) of individual fishes
load([TemplatePath3dpf,'/ForbrainInfo.mat']);
ForbrainInfo = ForbrainInfo;

%MeanMidbrain 
temp= load([TemplatePath3dpf,'/MeanMidbrain.mat']);
MeanMidBrain = temp.MeanMidBrain;
CenterMidBrain = temp.Centroidm;
MidBrainDistanceMap = temp.distancemapm;
BandMidBrain = temp.bandm;
%Meanforbrain makes not much sense as it variates largely
%it must be dependt on the midbrean and eyes




%% Brain Edge Template
% has te be removed later on
%a more normal function would be better
temp = load([TemplatePath3dpf,'/EdgeTemplate.mat']);  %polar transforms and edge template
EdgeTemplate1 = temp.EdgeTemplate;


%% Spot Parameters
% generated with SpotTemplate2
load([TemplatePath3dpf,'/SpotTemplate.mat'],'SpotTemplateVar');
SpotContrastVector = SpotTemplateVar.SpotContrastVector;
SpotVectorArrayProbability = SpotTemplateVar.SpotVectorArrayProbability;
SVAP_index = SpotTemplateVar.SVAP_index;
%clear SpotTemplateVar

%% Machine learning SpotParameters
load([TemplatePath3dpf,'/SpotTemplate.mat'],'SpotTemplateVar');
ClassifierInfo = load([SourcePath, filesep, 'Microglia classifiers', filesep ,'pksvc5050.mat'],'wp');
Classifier = ClassifierInfo.wp;



%%

filename = 'Template3dpf'

save([SourcePath, '/', filename, '.mat'], 'Template')
save([SourcePath, '/', filename, '.mat'], 'Size', '-append')
save([SourcePath, '/', filename, '.mat'], 'ref_temp', '-append')

save([SourcePath, '/', filename, '.mat'], 'EyeRegion1', '-append')
save([SourcePath, '/', filename, '.mat'], 'EyeRegion2', '-append')
save([SourcePath, '/', filename, '.mat'], 'EyeMask1', '-append')
save([SourcePath, '/', filename, '.mat'], 'EyeMask2', '-append')
save([SourcePath, '/', filename, '.mat'], 'EyeCenter1', '-append')
save([SourcePath, '/', filename, '.mat'], 'EyeCenter2', '-append')

save([SourcePath, '/', filename, '.mat'], 'MidbrainInfo', '-append')
save([SourcePath, '/', filename, '.mat'], 'ForbrainInfo', '-append')
save([SourcePath, '/', filename, '.mat'], 'MeanMidBrain', '-append')
save([SourcePath, '/', filename, '.mat'], 'CenterMidBrain', '-append')
save([SourcePath, '/', filename, '.mat'], 'MidBrainDistanceMap', '-append')
save([SourcePath, '/', filename, '.mat'], 'BandMidBrain', '-append')
save([SourcePath, '/', filename, '.mat'], 'EdgeTemplate1', '-append')

save([SourcePath, '/', filename, '.mat'], 'SpotContrastVector', '-append')
save([SourcePath, '/', filename, '.mat'], 'SpotVectorArrayProbability', '-append')
save([SourcePath, '/', filename, '.mat'], 'SVAP_index', '-append')

save([SourcePath, '/', filename, '.mat'], 'Classifier', '-append')

CompleteTemplate = load([SourcePath, '/', 'Template3dpf', '.mat'])



