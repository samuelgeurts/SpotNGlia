
%% select load and save paths
%clear all; 
%close all;

% SET specific image number to compute
stacknumbers = 0; %0 for processing all images
savefig_TF = 0;
savepar_TF = 0;
parforArg = 0; %0 for no parrallel computing, 2 for two core parallel computing
wbar_TF = 0;


if wbar_TF
    h = waitbar(0,'load files','Name','ZebraFish Complete Algorithm');
end
%assigns all paths needed
PathsComplete('bp','tp','pp','ed','rg','br','sp','rois','roib','tp3')

% Load all Template Information
%CompleteTemplate = LoadTemplateLink2(3,TemplatePath,Basepath);
CompleteTemplate = LoadTemplateLink3(TemplatePath3dpf);
ref_temp = imref2d(CompleteTemplate.Size);

%% step2: Extended Dept of Field
if wbar_TF;waitbar(1/5,h,'Extended Dept of Field');end


sng_MultiDir(who,'PreprocessionPath')
load([PreprocessionPath,'/stackinfo.mat'],'stackinfo')
load([PreprocessionPath,'/zfinput.mat'],'zfinput')
edoutput = cell(1,numel(stackinfo));

    zfinput = sng_zfinput(zfinput,0,'ExtendedDeptOfField','edof','variancedisksize',7,'moderate'); %sigma for bandpassfilter (lowpass)    

if stacknumbers == 0 ;stacknumbers = 1:numel(stackinfo);end
parfor (k2 = stacknumbers,parforArg)
    CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k2).stackname,'.tif']);       
    
    [Icombined,edoutput{k2}] = ExtendedDeptofFieldLink2(CorrectedSlice,zfinput);
    
    if savefig_TF
        imwrite(Icombined,[ExtendedDeptOfFieldPath,'/',stackinfo(k2).stackname,'-ExtendedDeptOfField.tif'], 'WriteMode', 'overwrite',  'Compression','none'); 
    end
end
[stackinfo.ExtendedDeptOfField] = edoutput{:};

if savepar_TF
    save([ExtendedDeptOfFieldPath,'/stackinfo.mat'],'stackinfo');
    save([ExtendedDeptOfFieldPath,'/zfinput.mat'],'zfinput');    
end

%% step 3: Registration
if wbar_TF;waitbar(2/5,h,'Registration');end

sng_MultiDir(who,'ExtendedDeptOfFieldPath','PreprocessionPath','TemplatePath')
load([ExtendedDeptOfFieldPath,'/stackinfo.mat'],'stackinfo')
load([ExtendedDeptOfFieldPath,'/zfinput.mat'],'zfinput')
rgoutput = cell(1,numel(stackinfo));

    zfinput = sng_zfinput(zfinput,0,'Registration','BackgroundRemoval','Method','TriangleSmooth','?');
    zfinput = sng_zfinput(zfinput,0,'Registration','BackgroundRemoval','Smooth',1,'?');  
    zfinput = sng_zfinput(zfinput,0,'Registration','BackgroundRemoval','ChannelMethod','cuboid','?');  
    zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','AngleSteps1',100,'?');  
    zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','Scale1',1/16,'?');  
    zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','AngleRange2',0.01,'?');  
    zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','AngleSteps2',50,'?');
    zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','FinalCenteredYCrop',1000,'?');
    zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','ScaleRange',[0.6,1.3],'?');      
    zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','ScaleSteps',200,'?');      
    zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','Scale2',1/16,'?');      
    zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','RemoveTail',600,'');   
    zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','Scale3',1/4,''); %lowered initial scale to increase computation speed for correlation
    zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','AffineMethod','translation',''); %affine or translation
    zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','Levels',3,''); %number of scaled levels correlation is performed
    zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','Iterations',40,''); %number of iterations iat is performed

if stacknumbers == 0 ;stacknumbers = 1:numel(stackinfo);end
parfor (k3 = stacknumbers,parforArg)
    k3
    CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k3).stackname,'.tif']);           
    Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k3).ExtendedDeptOfField.IndexMatrix);
    %Icombined = imread([ExtendedDeptOfFieldPath,'/',stackinfo(k1).stackname,'-ExtendedDeptOfField.tif']);
    
    [Ialligned,rgoutput{k3}] = AllignmentLink5(Icombined,CompleteTemplate,zfinput);
    if savefig_TF 
    imwrite(Ialligned,[RegistratedPath,'/',stackinfo(k3).stackname,'-Registration.tif'], 'WriteMode', 'overwrite',  'Compression','none'); 
    end
end

[stackinfo.Registration] = rgoutput{:};
if savepar_TF
    save([RegistratedPath,'/stackinfo.mat'],'stackinfo')
    save([RegistratedPath,'/zfinput.mat'],'zfinput');    
end
    
%% step 4: BrainSegmentation
if wbar_TF;waitbar(3/5,h,'Brain Segmentation');end

sng_MultiDir(who,'RegistratedPath','PreprocessionPath','TemplatePath')
load([RegistratedPath,'/stackinfo.mat'],'stackinfo')
load([RegistratedPath,'/zfinput.mat'],'zfinput')
broutput = cell(1,numel(stackinfo));

if stacknumbers == 0 ;stacknumbers = 1:numel(stackinfo);end
parfor (k4 = stacknumbers,parforArg)

    CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k4).stackname,'.tif']);           
    Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k4).ExtendedDeptOfField.IndexMatrix);
    %Icombined = imread([ExtendedDeptOfFieldPath,'/',stackinfo(k4).stackname,'-ExtendedDeptOfField.tif']);
    tform_1234 = stackinfo(k4).Registration(strcmp({stackinfo(k4).Registration.name},'tform_complete')).value
    Ialligned = imwarp(Icombined,tform_1234,'FillValues',255,'OutputView',ref_temp);

    [Ibrain,broutput{k4}] = MidBrainDetectionLink3(Ialligned,CompleteTemplate,zfinput);
   
    
    if savefig_TF
        imwrite(Ibrain,[BrainPath,'/',stackinfo(k4).stackname,'-BrainSegmentation.tif'], 'WriteMode', 'overwrite',  'Compression','none'); 
    end
end

[stackinfo.BrainSegmentation] = broutput{:};
if savepar_TF
    save([BrainPath,'/stackinfo.mat'],'stackinfo')
    save([BrainPath,'/zfinput.mat'],'zfinput');    
end
    
%% step 5: Spot Detection
if wbar_TF;waitbar(4/5,h,'Spot Detection');end

SpotPath = uigetdir(SpotPath)

sng_MultiDir(who,'BrainPath','PreprocessionPath','TemplatePath','SpotPath')
load([BrainPath,'/stackinfo.mat'],'stackinfo')
load([BrainPath,'/zfinput.mat'],'zfinput')

ColorToGrayVector = CompleteTemplate.SpotContrastVector;
ColorToGrayVector = [0;1;0]

    zfinput = sng_zfinput(zfinput,0,'SpotDetection','RgbToGray','ColorToGrayVector',ColorToGrayVector,'');  %select color channel [0 1 0], 
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleBase',0.5,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleLevels',7,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','Kthreshold',0,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPlevels',5:7,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPthreshold',256,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MinSpotSize',14,''); % size selection 
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MaxSpotSize',430,''); % size selection 
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MinProbability',0.1,''); %color selection

spoutput = cell(1,numel(stackinfo));

if stacknumbers == 0 ;stacknumbers = 1:numel(stackinfo);end
parfor (k5 = stacknumbers,parforArg)
    
    disp(k5);
    
    CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k5).stackname,'.tif']);           
    Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k5).ExtendedDeptOfField.IndexMatrix);
    %Icombined = imread([ExtendedDeptOfFieldPath,'/',stackinfo(k5).stackname,'-ExtendedDeptOfField.tif']);
    tform_1234 = stackinfo(k5).Registration(strcmp({stackinfo(k5).Registration.name},'tform_complete')).value;
    Ialligned = imwarp(Icombined,tform_1234,'FillValues',255,'OutputView',ref_temp);
    %Ialligned = imread([RegistratedPath,'/',stackinfo(k5).stackname,'-Registration.tif']); 
    cmbr = stackinfo(k5).BrainSegmentation.BrainEdge;
    %Ibrain = poly2mask(cmbr(:,2),cmbr(:,1),size(Ialligned,1),size(Ialligned,2));    %
    %Ibrain = imread([BrainPath,'/',stackinfo(k5).stackname,'-BrainSegmentation.tif']); 

    [Ispots,spoutput{k5}] = SpotDetectionLink2(Ialligned,CompleteTemplate,cmbr,zfinput);   
    
    if savefig_TF
        imwrite(Ispots,[SpotPath,'/',stackinfo(k5).stackname,'-SpotSelection.tif'], 'WriteMode', 'overwrite',  'Compression','none'); 
    end
end

[stackinfo.SpotSelection] = spoutput{:};
if savepar_TF
    save([SpotPath,'/stackinfo.mat'],'stackinfo'); %#ok<*UNRCH>
    save([SpotPath,'/zfinput.mat'],'zfinput');    
end
%{
    
    NewSpotPath = [Basepath,'/5_spotdetected_colortest']
    mkdir(NewSpotPath)
    save([NewSpotPath,'/stackinfo.mat'],'stackinfo');
    save([NewSpotPath,'/zfinput.mat'],'zfinput');   

%}
    

%%
if wbar_TF;waitbar(1,h,'complete');end

if wbar_TF;wb = findall(0,'type','figure','tag','TMWWaitbar');
delete(wb)
end

%load([SpotPath,'/stackinfo.mat'],'stackinfo')
%SpotSelection = sng_SubStruct({stackinfo.SpotSelection});
