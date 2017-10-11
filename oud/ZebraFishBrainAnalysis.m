
%sng_MultiDir(who,'PreprocessionPath','TemplatePath','BrainPath','RoiBrainPath')
PathsComplete('pp','roib','tp3','bp','br')


%add custom input brain
%
%function brainvalidation(stackinfo)
%add reference to stackinfo - template
%

load([BrainPath,'/stackinfo.mat'],'stackinfo');
load([BrainPath,'/zfinput.mat'],'zfinput');

CompleteTemplate = LoadTemplateLink3(TemplatePath3dpf);

saveTF = 0;

%preallocation
Jaccard = zeros(50,1);
cmbr = cell(50,1);
mbr2 = cell(50,1);


%%

for k1 = 1:numel(stackinfo)
    disp(k1);
    
    RoiBrain = ReadImageJROI([RoiBrainPath,'/',stackinfo(k1).stackname,'.zip']);
    tform_1234 = stackinfo(k1).Registration(strcmp({stackinfo(k1).Registration.name},'tform_complete')).value; 
    cmbr{k1} = stackinfo(k1).BrainSegmentation.BrainEdge;

    CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k1).stackname,'.tif']);           
    Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k1).ExtendedDeptOfField.IndexMatrix);
    Ialigned = imwarp(Icombined,tform_1234,'FillValues',255,'OutputView',CompleteTemplate.ref_temp);    

    
    mbr = sng_roicell2poly(RoiBrain,1);
    
    %transform coordinates to template 
    [mbr2{k1}(:,1),mbr2{k1}(:,2)] = transformPointsForward(tform_1234,mbr(:,1),mbr(:,2));
    mbr2{k1} = double(mbr2{k1});
    
    %{
    figure;imshow(uint8(Icombined))
    hold on;plot(mbr(:,1),mbr(:,2))
    figure;imshow(uint8(Ialigned))
    hold on;plot(mbr2{k1}(:,1),mbr2{k1}(:,2))
    %}    
    
    
%% compute Jaccard parameter
    
    a = poly2mask(cmbr{k1}(:,2),cmbr{k1}(:,1),1072,1424);
    b = poly2mask(mbr2{k1}(:,1),mbr2{k1}(:,2),1072,1424);

    %{
    figure;imagesc(a & b);
    figure;imagesc(a | b);
    figure;imagesc(1/2 * (a & b) + 1/4*(a - b) + 1/2*(a | b));sng_imfi
    figure;imagesc(a - b)
    %}
    %{
    q1 = ((a - b));
    q1(q1 < 0) = 0;
    q2 = ((b - a));
    q2(q2 < 0) = 0;    
    r = cat(3,(-60 * q2),(-40 * (q1+q2)),(-40 * q1));
    figure;imagesc(uint8(double(Ialligned) + r));sng_imfix
    hold on
    plot(mbr2(:,1),mbr2(:,2),'LineWidth',2)
    plot(cmbr(:,2),cmbr(:,1),'LineWidth',2)
%}
    
    num = (a & b);
    den = (a | b);

    Jaccard(k1) = sum(num(:))/sum(den(:));

end

Q = quantile(Jaccard,[.25 .5 .75]);
LF = Q(1) - 1.5 * (Q(3)-Q(1));
HF = Q(3) + 1.5 * (Q(3)-Q(1));

%fill the Validation structure
[BrainValidation.ImageInfo(1:numel(stackinfo)).ComputedMidBrain] = cmbr{:};
[BrainValidation.ImageInfo(1:numel(stackinfo)).AnnotatedMidBrainRegistrated] = mbr2{:};
temp = num2cell(Jaccard);[BrainValidation.ImageInfo(1:numel(stackinfo)).Jaccard] = temp{:};
BrainValidation.zfinput = zfinput
BrainValidation.FiveNumberSummaryJaccard = [LF,Q,HF]
BrainValidation.MeanJaccard = mean(Jaccard)


if saveTF
    save([Basepath,'/Validation Brain and Spot','/brainvalidation.mat'],'stackinfo')
end


figure;boxplot(Jaccard)
figure;plot(sort(Jaccard,'descend'))







%{
[stackinfo.SpotSelection] = spoutput{:};
brain = [stackinfo.BrainSegmentation]
Jaccard = [brain.Jaccard]'
%}

