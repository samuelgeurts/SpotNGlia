%Spot Opimization wavelet parameters

%To test the spotdetection performance, we use the ground truth about the
%brain region, which is segmentated by hand.
clear all
close all

limit = 10
stacknumbers = 20; %0 for processing all images

PathsComplete('bp','tp','pp','ed','rg','br','sp','rois','roib','tp3')
sng_MultiDir(who,'BrainPath','PreprocessionPath','TemplatePath','SpotPath')
SpotPath = '/Users/samuelgeurts/Dropbox/20170327_3dpf_test/5_spotdetected_color';
SpotPath = '/Users/samuelgeurts/Dropbox/20170327_3dpf_test/5_spotdetected 20170927';


load([SpotPath,'/zfinput.mat'],'zfinput')
load([SpotPath,'/stackinfo.mat'],'stackinfo');

% Load all Template Information
%CompleteTemplate = LoadTemplateLink2(3,TemplatePath);
CompleteTemplate = LoadTemplateLink3(TemplatePath3dpf);
ColorToGrayVector = CompleteTemplate.SpotContrastVector;


ref_temp = imref2d(CompleteTemplate.Size);
load([Basepath,'/Template Spot/SpotTemplate.mat'],'SpotTemplateVar');
CompleteTemplate.SpotContrastVector = SpotTemplateVar.SpotContrastVector;
CompleteTemplate.SpotVectorArrayProbability = SpotTemplateVar.SpotVectorArrayProbability;
CompleteTemplate.SVAP_index = SpotTemplateVar.SVAP_index;
clear SpotTemplateVar


if ~exist('zfinput','var');
    zfinput = struct;
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','RgbToGray','ColorToGrayMethod',ColorToGrayVector,'');  %select color channel
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleBase',0.5,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleLevels',8,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','Kthreshold',0,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPlevels',4:9,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPthreshold',3000,'');  
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MinSpotSize',30,''); % size selection 
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MaxSpotSize',500,''); % size selection 
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MinProbability',0.01,''); %color selection

end




%


if stacknumbers == 0 ;stacknumbers = 1:numel(stackinfo);end;
for k5 = stacknumbers
    
    CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k5).stackname,'.tif']);           
    Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k5).ExtendedDeptOfField.IndexMatrix);
    %Icombined = imread([ExtendedDeptOfFieldPath,'/',stackinfo(k5).stackname,'-ExtendedDeptOfField.tif']);
    tform_1234 = stackinfo(k5).Registration(strcmp({stackinfo(k5).Registration.name},'tform_complete')).value;
    Ialligned = imwarp(Icombined,tform_1234,'FillValues',255,'OutputView',ref_temp);
    %Ialligned = imread([RegistratedPath,'/',stackinfo(k5).stackname,'-Registration.tif']); 
    %cmbr = stackinfo(k5).BrainSegmentation.BrainEdge;    
    ambr = stackinfo(k5).BrainSegmentation.AnnotatedMidBrainRegistrated;
 

    % annotated spots
    RoiMicroglia = ReadImageJROI([RoiMicrogliaPath,'/',stackinfo(k5).stackname,'.roi']);
    spota = RoiMicroglia.mfCoordinates;
    spots = RoiMicroglia.vnSlices;
    [SpotAnnX,SpotAnnY] = transformPointsForward(tform_1234,spota(:,1),spota(:,2));
    %remove annotated spots if they are no insite the by hand segmentated
    %brain region. Small error could occur because diferent people
    %annotated spots and brain region
    [in,~] = inpolygon(SpotAnnY,SpotAnnX,ambr(:,2),ambr(:,1));
    SpotAnnX = SpotAnnX(in);
    SpotAnnY = SpotAnnY(in);
       
    %% Spot comparison for differen wavelet parameters
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleLevels',10,'');      
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MinSpotSize',1,''); % size selection 
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','SpotSelection','MinProbability',0,''); %color selection
    zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPlevels',4:10,'');  
    
    MPthresholdList = [0.5 1 2 4 8 16 32 64 128 256 512 1024 2048 4096];
    for l1 = 1:numel(MPthresholdList)
    
        %zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPlevels',3:8,'');  
        zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPthreshold',MPthresholdList(l1),'');      
   
        %use by hand segemented brain region ambr instead of computed cmbr
        [~,spoutput{k5}] = SpotDetectionLink2(Ialligned,CompleteTemplate,fliplr(ambr),zfinput); 
        
        %% computed spots
        Spotinfo = spoutput{k5}(3).value;
        %SpotinfoS = Spotinfo([Spotinfo.Insite]); %selection of spotinfo based on conditions    
        SpotinfoS = Spotinfo([Spotinfo.MinProbability] &... %selection of spotinfo based on conditions
        [Spotinfo.Insite] &...
        [Spotinfo.LargerThan]); % &...
        %[Spotinfo.SmallerThan])
        SpotCom = reshape([SpotinfoS.Centroid],2,numel(SpotinfoS))';
        
        
        %Spot compare                
        [Correct,FalsePos,FalseNeg,~] = sng_CoordinateMatching(SpotCom,[SpotAnnX SpotAnnY],limit);
        TruePosList(l1) = size(Correct,1);
        FalsePosList(l1) = size(FalsePos,1);
        FalseNegList(l1) = size(FalseNeg,1);
         
        
        figure;imagesc(uint8(Ialligned));sng_imfix
        hold on
        plot(ambr(:,1),ambr(:,2),'LineWidth',2)
        %plot(cmbr(:,2),cmbr(:,1),'LineWidth',2)
        if exist('Correct');scatter(Correct(:,1), Correct(:,2),400,'LineWidth',2,'MarkerEdgeColor',1/255*[75 255 75]);end
        if exist('FalsePos');scatter(FalsePos(:,1), FalsePos(:,2),400,'LineWidth',2,'MarkerEdgeColor',1/255*[255 75 75]);end
        if exist('FalseNeg');scatter(FalseNeg(:,1), FalseNeg(:,2),400,'LineWidth',2,'MarkerEdgeColor',1/255*[75 75 255]);end
        drawnow        
        
        
        
    end
    idx = 1:numel(MPthresholdList);
    figure;plot(idx,TruePosList,idx,FalsePosList,idx,FalseNegList)
    
    
    
end





