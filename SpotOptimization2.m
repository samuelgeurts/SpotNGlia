%Spot Opimization wavelet parameters

%To test the spotdetection performance, we use the ground truth about the
%brain region, which is segmentated by hand.

obj = BrainSpotValidation;
obj = SpotVal(obj);


PathsComplete('tp3','pp')
load([obj.SpotPath,'/stackinfo.mat'],'stackinfo');
stacknumbers = 1:numel(stackinfo);
CompleteTemplate = LoadTemplateLink3(TemplatePath3dpf);

Ialigned = cell(obj.Numb,1);
spoutput = cell(obj.Numb,1);
%load all images
parfor k5 = stacknumbers
    disp(num2str(k5))
    CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k5).stackname,'.tif']);           
    Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k5).ExtendedDeptOfField.IndexMatrix);
    %Icombined = imread([ExtendedDeptOfFieldPath,'/',stackinfo(k5).stackname,'-ExtendedDeptOfField.tif']);
    tform_1234 = stackinfo(k5).Registration(strcmp({stackinfo(k5).Registration.name},'tform_complete')).value;
    Ialigned{k5} = imwarp(Icombined,tform_1234,'FillValues',255,'OutputView',CompleteTemplate.ref_temp);
end

MPthresholdList = [0.5 1 2 4 8 16 32 64 128 256 512 1024 2048 4096];
ColorToGrayVector = CompleteTemplate.SpotContrastVector;

    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','RgbToGray','ColorToGrayVector',ColorToGrayVector,'');  %select color channel [0 1 0], 
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','Wavelet','ScaleBase',0.5,'');      
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','Wavelet','ScaleLevels',9,'');  
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','Wavelet','Kthreshold',0,'');  
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','MultiProduct','MPlevels',4:9,'');  
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','MultiProduct','MPthreshold',50,'');  
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinSpotSize',30,''); % size selection 
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MaxSpotSize',500,''); % size selection 
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinProbability',0.01,''); %color selection

    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinSpotSize',1,''); % size selection 
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MaxSpotSize',2000,''); % size selection 
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinProbability',0,''); %color selection

for l1 = 7:numel(MPthresholdList)

    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','MultiProduct','MPthreshold',MPthresholdList(l1),'');      
    

    %compute Spots
    for k5 = stacknumbers
        disp([num2str(l1),' ',num2str(k5)])
        BrainAnn = obj.BrainInfo(k5).AnnotatedMidBrainRegistrated;
        [~,spoutput{k5}] = SpotDetectionLink3(Ialigned{k5},CompleteTemplate,fliplr(BrainAnn),obj.zfinput);         
        obj.SpotParameters{k5} = spoutput{k5}(3).value;              
    end
    
    obj = obj.SpotVal;

    %selects spot sizes and color probabilities, correct & incorrect
    AreaC = [];
    AreaI = [];
    ColorC = [];
    ColorI = [];
    
    for k5 = stacknumbers
        Spotpar = obj.SpotParameters{k5}([obj.SpotParameters{k5}.Insite]); %selects only insite brain region                
        SpotCom = reshape([Spotpar.Centroid],2,numel(Spotpar))';
        [a] = ismember(SpotCom,obj.SpotInfo(k5).CorrectSpots);

        AreaC = [AreaC,[Spotpar(a(:,1)).Area]];
        AreaI = [AreaI,[Spotpar(~a(:,1)).Area]];
        
        ColorC = [ColorC,[Spotpar(a(:,1)).ColorProbability]];
        ColorI = [ColorI,[Spotpar(~a(:,1)).ColorProbability]];   
    end
       
    %figure;imagesc(Spotpar(1).Image)

        x = linspace(0,500,100)
        Eh1=sng_chistcount(AreaC,x);      %create bins
        Bh1=sng_chistcount(AreaI,x);      %
        sng_DistDifHistPlot((x(1:end-1)+((x(2)-x(1))/2)),Eh1,Bh1);
         
        x = linspace(0,1.2,100)
        Eh1=sng_chistcount(ColorC,x);      %create bins
        Bh1=sng_chistcount(ColorI,x);      %
        sng_DistDifHistPlot((x(1:end-1)+((x(2)-x(1))/2)),Eh1,Bh1);
         

        
        DS = [[AreaC';AreaI'],[ColorC';ColorI']]
        lab = [ones(numel(AreaC),1);2*ones(numel(AreaI),1)]     
        a = prdataset(DS,lab)     
        a = a*scalem(a,'variance') %scaling
        
        figure;
        scatterd(a)
        wp = ldc(a)
        plotc(wp)
        
        err{l1} = a*wp*testc %error of testset c on trained classifier wb

        
        
obj.show
obj.saveit(['meas_',num2str(l1)])





%{
sng_zfinputAssign(obj.zfinput,'SpotDetection')
SpotOpt.MPlevels = MPlevels
SpotOpt.MPthreshold = MPthreshold
SpotOpt.MinSpotSize = MinSpotSize
SpotOpt.MaxSpotSize = MaxSpotSize
SpotOpt.MinProbability = MinProbability

SpotOpt.Precision = [obj.SpotInfo.Precision]
SpotOpt.Recall = [obj.SpotInfo.Recall]
SpotOpt.F1score = [obj.SpotInfo.F1score]
%}
end





