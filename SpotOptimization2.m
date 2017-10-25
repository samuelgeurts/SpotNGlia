

%Spot Opimization wavelet parameters

%To test the spotdetection performance, we use the ground truth about the
%brain region, which is segmentated by hand.

obj = BrainSpotValidation;
obj = SpotVal(obj);

fig_tf = false;
fig_tf2 = 1;



PathsComplete('tp3','pp')
load([obj.SpotPath,'/stackinfo.mat'],'stackinfo');
stacknumbers = 1:numel(stackinfo);
CompleteTemplate = LoadTemplateLink3(TemplatePath3dpf);

Ialigned = cell(obj.Numb,1);
spoutput = cell(obj.Numb,1);
%load all images
for k5 = stacknumbers
    disp(num2str(k5))
    CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k5).stackname,'.tif']);           
    Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k5).ExtendedDeptOfField.IndexMatrix);
    %Icombined = imread([ExtendedDeptOfFieldPath,'/',stackinfo(k5).stackname,'-ExtendedDeptOfField.tif']);
    tform_1234 = stackinfo(k5).Registration(strcmp({stackinfo(k5).Registration.name},'tform_complete')).value;
    Ialigned{k5} = imwarp(Icombined,tform_1234,'FillValues',255,'OutputView',CompleteTemplate.ref_temp);
end

MPthresholdList = [2 4 8 16 32 64 128 256];
MPlevelList = {3:7,3:8,3:9,4:7,4:8,4:9,5:7,5:8,5:9};

%MPthresholdList = [8 16];
%MPlevelList = {5:8};



ColorToGrayVector = CompleteTemplate.SpotContrastVector;
%ColorToGrayVector = [0;1;0];

    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','RgbToGray','ColorToGrayVector',ColorToGrayVector,'');  %select color channel [0 1 0], 
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','Wavelet','ScaleBase',0.5,'');      
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','Wavelet','ScaleLevels',9,'');  
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','Wavelet','Kthreshold',0,'');  
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','MultiProduct','MPlevels',4:9,'');  
    obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','MultiProduct','MPthreshold',50,'');  
    %obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinSpotSize',30,''); % size selection 
    %obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MaxSpotSize',500,''); % size selection 
    %obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinProbability',0.01,''); %color selection


%try out different values for Multiproduct
for l2 = 1:numel(MPlevelList)  
    for l1 = 1:numel(MPthresholdList)

        linen = ((l2 - 1) * numel(MPthresholdList) + l1); %variable 2 indicate at which line to save info
        
        obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinSpotSize',1,''); % size selection 
        obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MaxSpotSize',2000,''); % size selection 
        obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinProbability',0,''); %color selection
        obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','MultiProduct','MPthreshold',MPthresholdList(l1),'');      
        obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','MultiProduct','MPlevels',MPlevelList{l2},'');  

        %compute Spots
        for k5 = stacknumbers
            fprintf([num2str(l2),'/',num2str(l1),'/',num2str(k5),' '])
            BrainAnn = obj.BrainInfo(k5).AnnotatedMidBrainRegistrated;
            [~,spoutput{k5}] = SpotDetectionLink3(Ialigned{k5},CompleteTemplate,fliplr(BrainAnn),obj.zfinput);         
            obj.SpotParameters{k5} = spoutput{k5}(3).value;              
        end

        obj = obj.SpotVal;
        
        disp([num2str(l2),' ',num2str(l1),' ',num2str(mean([obj.SpotInfo.F1score]))]);
        
        load([obj.SavePath,'/','SpotOptList2','.mat'],'SpotOptList2');  
        SpotOptList2 = AddToSpotOptList(SpotOptList2,obj,(linen*2)-1); %stores all relevant info into the next field of SpotOptList
        save([obj.SavePath,'/','SpotOptList2','.mat'],'SpotOptList2');  

        if fig_tf    
            obj.show;
        end

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
        %{
        if fig_tf
            figure;imagesc(Spotpar(3).Image)    
        end
        %}

        %% Compute Threshold for area, larger than and smaller than!
        x = linspace(0,1000,1000);                        %histogram bins


        AreaHistCorrect = histcounts(AreaC,x);      %correct spot area
        AreaHistIncorrect = histcounts(AreaI,x);    %incorrect spot area
        %finds threshold, index is part of the right site of the histogram
        for j = 1:size(AreaHistCorrect,2)
            FR(j) = sum(AreaHistIncorrect(1:j-1)) + sum(AreaHistCorrect(j:end));
            FL(j) = sum(AreaHistCorrect(1:j-1)) + sum(AreaHistIncorrect(j:end));            
        end

        [mxR,indexR] = max(FR);
        MinSpotSize = x(indexR); %every area equal or larger than "LargerThan" is large enough                
        [mxL,indexL] = max(FL(indexR:end));
        MaxSpotSize = x(indexL+indexR); %every area equal or larger than "SmallerThan" is small enough

        if fig_tf2
            sng_DistDifHistPlot((x(1:end-1)+((x(2)-x(1))/2)),AreaHistCorrect,AreaHistIncorrect);      
            hold on;line([MinSpotSize MinSpotSize],get(gca,'Ylim'),'color',[0 0 0]);
            hold on;line([MaxSpotSize MaxSpotSize],get(gca,'Ylim'),'color',[0 0 0]);
        end

        %% Compute Threshold for color probability
        x = linspace(0,1.2,1000)
        ColorHistCorrect = histcounts(ColorC,x);      %create bins
        ColorHistIncorrect = histcounts(ColorI,x);      %

        for j = 1:size(ColorHistCorrect,2)
            FR(j) = sum(ColorHistIncorrect(1:j-1)) + sum(ColorHistCorrect(j:end));
        end
        [mxR,indexR] = max(FR);
        MinProbability = x(indexR); %every area equal or larger than "LargerThan" is large enough                

        if fig_tf2
            sng_DistDifHistPlot((x(1:end-1)+((x(2)-x(1))/2)),ColorHistCorrect,ColorHistIncorrect);
            hold on;line([MinProbability MinProbability],get(gca,'Ylim'),'color',[0 0 0]);
        end

     %using a machine learning technique would be an improvement
        %DS = [[AreaC';AreaI'],[ColorC';ColorI']]
        %lab = [ones(numel(AreaC),1);2*ones(numel(AreaI),1)]     
        %a = prdataset(DS,lab)     
        %a = a*scalem(a,'variance') %scaling        
        %figure;
        %scatterd(a)
        %wp = parzenc(a,0.25)
        %plotc(wp)
        %err{l1} = a*wp*testc %error of testset c on trained classifier wb

        %% assign new threshold logical array to obj.SpotParameters
        for k5 = stacknumbers
            temp = num2cell([obj.SpotParameters{k5}.Area] >= MinSpotSize);
            [obj.SpotParameters{k5}.LargerThan] = temp{:};

            temp = num2cell([obj.SpotParameters{k5}.Area] <= MaxSpotSize);
            [obj.SpotParameters{k5}.SmallerThan] = temp{:};

            temp = num2cell([obj.SpotParameters{k5}.ColorProbability] >= MinProbability);
            [obj.SpotParameters{k5}.MinProbability] = temp{:};
        end

        obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinSpotSize',MinSpotSize,'');      % size selection 
        obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MaxSpotSize',MaxSpotSize,'');      % size selection 
        obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinProbability',MinProbability,'');% color selection

        obj = obj.SpotVal;
        
        load([obj.SavePath,'/','SpotOptList2','.mat'],'SpotOptList2');  
        SpotOptList2 = AddToSpotOptList(SpotOptList2,obj,linen*2); %stores all relevant info into the next field of SpotOptList
        save([obj.SavePath,'/','SpotOptList2','.mat'],'SpotOptList2');    

        if fig_tf    
            obj.show
        end
        
        disp(num2str(mean([obj.SpotInfo.F1score])));
   
        
        %obj.saveit(['meas_',num2str(l1)])

        %{
        if ~exist('savename','var')
            savename = inputname(1);
        end
        save([obj.SavePath,'/',savename,'.mat'],'obj');  
        %}
    end
end

%save([obj.SavePath,'/','savename','.mat'],'obj');  

%{

%show SpotOptList results

SF1 = [SpotOptList.MeanF1score]
SF1s = [SpotOptList.StdF1score]

q = 2*(1:numel(SpotOptList)/2)

figure;
plot(SF1(q-1))
hold on;plot(SF1(q))
hold on;plot(SF1s(q-1))
hold on;plot(SF1s(q))
legend('show')

%}


%qq = SpotOptList([SpotOptList.MeanF1score] >= 0.7)












