classdef SNGTemplate
    %class for generating template object. Based on the function MakeTemplateMat which it replaces. Running the function load the needed template
    %variables and creates the template on a chosen location. In future the load function has to be replaced by the functions
    %really generating the template instead of loading.
    
    properties
        
        fileName = 'Template3dpf'
        sourcePath
        TemplatePath3dpf
        
        %Template
        Template
        Size
        ref_temp
        
        %eyes
        EyeRegion1
        EyeRegion2
        EyeMask1
        EyeMask2
        EyeCenter1
        EyeCenter2
        
        %Brain
        MidbrainInfo
        ForbrainInfo
        MeanMidBrain
        CenterMidBrain
        MidBrainDistanceMap
        BandMidBrain
        EdgeTemplate1
        
        %Spot
        SpotContrastVector
        SpotVectorArrayProbability
        SVAP_index
        Classifier
        
    end
    
    
    methods
        function objt = SNGTemplate
        end      
        function objt = load(objt)
            %Based on the function MakeTemplateMat. Running the function load the needed template
            %variables and creates the template on a chosen location. In future the load function has to be replaced by the functions
            %really generating the template instead of loading.
            
            objt.sourcePath = uigetdir;
            objt.TemplatePath3dpf = [objt.sourcePath, '/', 'Template 3 dpf'];
            
            %% Mean Fish Registration Template
            %based on older template, can be updated with 50 fish batch but lead
            %probably not to improvement. Source is not sure, TemplateGeneration?
            
            objt.Template = imread([objt.TemplatePath3dpf, '/template.tif']);
            objt.Size = size(objt.Template);
            objt.ref_temp = imref2d(objt.Size);
            
            
            %Eye region, annotated by Anouk H
            %annotated on the older template
            temp = sng_html2roi2([objt.TemplatePath3dpf, '/edof_eyes.html']);
            objt.EyeRegion1 = temp{1}{2};
            objt.EyeRegion2 = temp{1}{1};
            
            eye1 = objt.EyeRegion2;
            eye2 = objt.EyeRegion1;
            
            objt.EyeMask1 = poly2mask(eye1(:, 1), eye1(:, 2), objt.Size(1), objt.Size(2));
            objt.EyeMask2 = poly2mask(eye2(:, 1), eye2(:, 2), objt.Size(1), objt.Size(2));
            clear eye1 eye2
            
            [x1c, y1c] = sng_CenterOfMass(objt.EyeMask1);
            [x2c, y2c] = sng_CenterOfMass(objt.EyeMask2);
            
            objt.EyeCenter1 = [x1c, y1c];
            objt.EyeCenter2 = [x2c, y2c];
            
            
            %{
         figure;imagesc(eyem1+eyem2)
            %}
            
            %% Brain
            %generated with the function TemplateBrainParameters
            
            %Midbrain info (poly,area,centroid,box) of individual fishes
            temp = load([objt.TemplatePath3dpf, '/MidbrainInfo.mat']);
            objt.MidbrainInfo = temp.MidbrainInfo;
            
            %Forbrain info (poly,area,centroid,box) of individual fishes
            temp = load([objt.TemplatePath3dpf, '/ForbrainInfo.mat']);
            objt.ForbrainInfo = temp.ForbrainInfo;
            
            %MeanMidbrain
            temp = load([objt.TemplatePath3dpf, '/MeanMidbrain.mat']);
            objt.MeanMidBrain = temp.MeanMidBrain;
            objt.CenterMidBrain = temp.Centroidm;
            objt.MidBrainDistanceMap = temp.distancemapm;
            objt.BandMidBrain = temp.bandm;
            %Meanforbrain makes not much sense as it variates largely
            %it must be dependt on the midbrean and eyes
            
            %% Brain Edge Template
            % has te be removed later on
            %a more normal function would be better
            temp = load([objt.TemplatePath3dpf, '/EdgeTemplate.mat']); %polar transforms and edge template
            objt.EdgeTemplate1 = temp.EdgeTemplate;

            %% Spot Parameters
            % generated with SpotTemplate2
            temp = load([objt.TemplatePath3dpf, '/SpotTemplate.mat'], 'SpotTemplateVar');
            objt.SpotContrastVector = temp.SpotTemplateVar.SpotContrastVector;
            objt.SpotVectorArrayProbability = temp.SpotTemplateVar.SpotVectorArrayProbability;
            objt.SVAP_index = temp.SpotTemplateVar.SVAP_index;
            %clear SpotTemplateVar
            
            %% Machine learning SpotParameters
            ClassifierInfo = load([objt.sourcePath, filesep, 'Microglia classifiers', filesep, 'pksvc5050.mat'], 'wp');
            objt.Classifier = ClassifierInfo.wp;
        end      
        function save(objt)
            save(strcat(objt.sourcePath, filesep, objt.fileName, '.mat'), 'objt');
        end
        
        function BrainTemplate
        end
        
    end
end
    

    
    
    