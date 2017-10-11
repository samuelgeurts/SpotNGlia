PathsComplete('pp','rois','roib','tp3','bp','tp','sp')

% add custum input spotpath (stackinfo)

SpotPath = '/Users/samuelgeurts/Dropbox/20170327_3dpf_test/5_spotdetected_color';

%\sng_MultiDir(who,'PreprocessionPath','TemplatePath','SpotPath','RoiMicrogliaPath')

load([SpotPath,'/zfinput.mat'],'zfinput')
load([SpotPath,'/stackinfo.mat'],'stackinfo');

CompleteTemplate = LoadTemplateLink3(TemplatePath3dpf);

%IMPORTANT PARAMETER, pixel limit is set at 
limit = 10;
Save_TF = 0;

%%    


for k1 = 5%:numel(stackinfo);
    disp(k1);
    
    clearvars -except PreprocessionPath TemplatePath SpotPath RoiMicrogliaPath limit CompleteTemplate ref_temp stackinfo k1 ...
        LinkDistance CorrectSpots nCorrect FalsePosSpots nFalsePos FalseNegSpots nFalseNeg Save_TF
     
    RoiBrain = ReadImageJROI([RoiBrainPath,'/',stackinfo(k1).stackname,'.zip']);    
    tform_1234 = stackinfo(k1).Registration(strcmp({stackinfo(k1).Registration.name},'tform_complete')).value;
    cmbr{k1} = stackinfo(k1).BrainSegmentation.BrainEdge;
       
    CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k1).stackname,'.tif']);           
    Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k1).ExtendedDeptOfField.IndexMatrix);
    Ialligned = imwarp(Icombined,tform_1234,'FillValues',255,'OutputView',CompleteTemplate.ref_temp);    

    ambr = stackinfo(k1).BrainSegmentation.AnnotatedMidBrainRegistrated;
    %cmbr = stackinfo(k1).BrainSegmentation.ComputedMidBrainRegistration;
    cmbr = stackinfo(k1).BrainSegmentation.BrainEdge;
    
        !!! add sng_roi2poly

    
    
    %annotated spots
    RoiMicroglia = ReadImageJROI([RoiMicrogliaPath,'/',stackinfo(k1).stackname,'.roi']);
    spota = RoiMicroglia.mfCoordinates;
    spots = RoiMicroglia.vnSlices;
    [SpotAnn(:,1),SpotAnn(:,2)] = transformPointsForward(tform_1234,spota(:,1),spota(:,2));
    
    %computed spotsimage
    
    infoselect = stackinfo(k1).SpotSelection(find([imageinfo.nextstack]~=2))

    
    Spotinfo =  stackinfo(k1).SpotSelection(strcmp({stackinfo(k1).SpotSelection.name},'SpotParameters')).value;
    %Spotinfo = stackinfo(k1).SpotSelection(3).value    
    
    SpotinfoS = Spotinfo([Spotinfo.MinProbability] &... %selection of spotinfo based on conditions
        [Spotinfo.Insite] &...
        [Spotinfo.LargerThan] &...
        [Spotinfo.SmallerThan])

    SpotCom = reshape([Spotinfo.Centroid],2,numel(Spotinfo))';

    %{
        figure;imagesc(uint8(Ialligned));sng_imfix
        hold on
        plot(ambr(:,1),ambr(:,2),'LineWidth',2)
        scatter(SpotAnn(:,1),SpotAnn(:,2))
    
        figure;imagesc(uint8(Ialligned));sng_imfix
        hold on    
        plot(cmbr(:,2),cmbr(:,1),'LineWidth',2)
        scatter(SpotCom(:,1),SpotCom(:,2))
    %}
    
    
    
    %% compute distance betwee every computed and annotated point
    [Correct,FalsePos,FalseNeg,link] = sng_CoordinateMatching(SpotCom,SpotAnn,limit);

    
%{
    figure;imagesc(uint8(Ialligned));sng_imfix
    hold on
    plot(ambr(:,1),ambr(:,2),'LineWidth',2)
    plot(cmbr(:,2),cmbr(:,1),'LineWidth',2)
    if exist('Correct');scatter(Correct(:,1), Correct(:,2),400,'LineWidth',2,'MarkerEdgeColor',1/255*[75 255 75]);end
    if exist('FalsePos');scatter(FalsePos(:,1), FalsePos(:,2),400,'LineWidth',2,'MarkerEdgeColor',1/255*[255 75 75]);end
    if exist('FalseNeg');scatter(FalseNeg(:,1), FalseNeg(:,2),400,'LineWidth',2,'MarkerEdgeColor',1/255*[75 75 255]);end
    drawnow
%}

%%
    %coordinates
    
    
    LinkDistance{k1} = link;

    if exist('Correct');
        CorrectSpots{k1} = Correct;
        nCorrect(k1) = size(Correct,1);
    else
        CorrectSpots{k1} = [];
        nCorrect(k1) = 0;
    end
    if exist('FalsePos');
        FalsePosSpots{k1} = FalsePos;
        nFalsePos(k1) = size(FalsePos,1);
    else
        FalsePosSpots{k1} = [];
        nFalsePos(k1) = 0;
    end
    if exist('FalseNeg');
        FalseNegSpots{k1} = FalseNeg;
        nFalseNeg(k1) = size(FalseNeg,1);
    else
        FalseNegSpots{k1} = [];
        nFalseNeg(k1) = 0;
    end
 
    
    
    
[stackinfo(k1).SpotAnalysis.LinkDistance] = link;
[stackinfo(k1).SpotAnalysis.CorrectSpots] = CorrectSpots{k1};
[stackinfo(k1).SpotAnalysis.nCorrect] = nCorrect(k1);
[stackinfo(k1).SpotAnalysis.FalsePosSpots] = FalsePosSpots{k1};
[stackinfo(k1).SpotAnalysis.nFalsePos] = nFalsePos(k1);
[stackinfo(k1).SpotAnalysis.FalseNegSpots] = FalseNegSpots{k1};
[stackinfo(k1).SpotAnalysis.nFalseNeg] = nFalseNeg(k1);



end

if Save_TF
    save([SpotPath,'/stackinfo.mat'],'stackinfo');
end

%SpotAnalysis = sng_SubStruct({stackinfo.SpotAnalysis})
%[[SpotAnalysis.nCorrect]',[SpotAnalysis.nFalsePos]',[SpotAnalysis.nFalseNeg]']


