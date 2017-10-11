function [Ialligned,Parameters] = AllignmentLink2(Icombined,CompleteTemplate)
% this function is a subfunction of "Main Function" which alligns a fish
% with a given template in a couple of steps.
% The function is strongly based on "HorizontalFish7"
%
% Version AllignmentLink2 
%       transforms fishes to template image
%       20170613 added crop before final translation 

template = CompleteTemplate.Template;

scales = linspace(0.6,1.3,200); %for scale and flip finding
I10 = Icombined;
        
%{
[xcom1,ycom1] = sng_CenterOfMassColor(imcomplement(I10),1);
shift = (fliplr(size(I10(:,:,1))/2)-[xcom1,ycom1]);
I11= imtranslate(I10,shift,'OutputView','same');
figure;imagesc(I11)
%}

%horizontal allignment
%remove background to compute center of mass in horizontalfish well
I20 = sng_RemoveBackgroundColor3(I10,'TriangleSmooth',1,'cuboid');
tform1 = sng_HorizontalFish3(I20,'precise'); %translation + rotation
[I30,ref30] = imwarp(I10,tform1,'FillValues',255);

%scale and leftrigh fitting with template
%cropping is needed to perform well for ScaleFish3 as it is assumed that fishes are centered.
%this crop is not interesting as the image I50 after scaling does not use the cropped image
[I40,ref40] = sng_fishcrop3(I30,ref30,[1000,1360]);
%crop1 = [ref40.XWorldLimits-ref30.XWorldLimits;ref40.YWorldLimits-ref30.YWorldLimits]

tform2 = sng_ScaleFish3(double(template),double(I40),scales,1/16);
    tform_12 = affine2d(tform1.T*tform2.T);
[I50,ref50] = imwarp(I10,tform_12,'FillValues',255);

[I51,ref51] = sng_fishcrop3(I50,ref50,[1000,1360]); %crop for better processing in iat algorithm
tform3 = affine2d([1,0,0;0,1,0;-ref51.XWorldLimits(1),-ref51.YWorldLimits(1),1]); %transform back to pixel coordinates

[tform4,I66,CorCoef] = sng_AllignFish2Template2(template,I51,1/4,'translation',3,40);

tform_1234 = affine2d(tform1.T*tform2.T*tform3.T*tform4.T)

Ialligned = imwarp(I10,tform_1234,'FillValues',255,'OutputView',imref2d(size(template(:,:,1))));


%{
figure;imagesc(I10)
figure;imagesc(ref30.XWorldLimits,ref30.YWorldLimits,uint8(I30));
figure;imagesc(ref40.XWorldLimits,ref40.YWorldLimits,uint8(I40));
figure;imagesc(uint8(I50));
figure;imagesc(ref50.XWorldLimits,ref50.YWorldLimits,uint8(I50));
figure;imagesc(uint8(I65));
figure;imagesc(uint8(template));
figure;imagesc(uint8(CompleteTemplate.MidBrainDistanceMap));

figure;
imagesc(uint8(Ialligned));colormap(gray);
hold on; plot(CompleteTemplate.EyeRegion{1}{1}(:,1),CompleteTemplate.EyeRegion{1}{1}(:,2))
hold on; plot(CompleteTemplate.EyeRegion{1}{2}(:,1),CompleteTemplate.EyeRegion{1}{2}(:,2))
plot(CompleteTemplate.MeanMidBrain(:,2),CompleteTemplate.MeanMidBrain(:,1))
scatter(CompleteTemplate.CenterMidBrain(1),CompleteTemplate.CenterMidBrain(2))

figure;imshowpair(template,uint8(Ialligned))

%}

Parameters.tform_hor = tform1;
Parameters.tform_scale = tform2;
Parameters.tform_transtopix = tform3;
Parameters.tform_trans = tform4;
Parameters.tform_complete = tform_1234;

Parameters.ref30 = ref30;
Parameters.ref40 = ref40;
Parameters.ref50 = ref50;

Parameters.TemplateCorrelation = CorCoef;







%% show background removal
%{
figure;imagesc(I10);axis off equal tight
set(gca,'position',[0 0 1 1],'units','normalized')

figure;imagesc(I20);axis off equal tight
set(gca,'position',[0 0 1 1],'units','normalized')

Hr1=histcounts(double(I10),[0:1:255]);
figure; bar(Hr1,'Facecolor',[0.8,0.4,0.4],'Edgecolor','none','BarWidth',1);
set(gca,'XLim',[0 256],'YLim',[0 100000],'FontSize',20);
set(gcf,'Color',[1 1 1])

HrBR1=histcounts(double(I20),[0:1:255]);
figure; bar(Hr1,'Facecolor',[0.8,0.4,0.4],'Edgecolor','none','BarWidth',1);
set(gca,'XLim',[0 256],'YLim',[0 100000],'FontSize',20);
set(gcf,'Color',[1 1 1])
%}


%% Horizontal allignent figures are programmed in sng_HorizontalFish3
%{
%}


%% Final Template Matching
%{
figure;imagesc(uint8(Ialligned))
set(gca,'position',[0 0 1 1],'units','normalized')
axis off tight equal

figure; imagesc(boolean(CompleteTemplate.MidBrainDistanceMap))

B = bwboundaries(CompleteTemplate.MidBrainDistanceMap)
figure;imagesc(uint8(Ialligned))
hold on; 
plot(B{1}(:,2),B{1}(:,1),'Color',[0 176/255 240/255],'LineWidth',2)
h = plot(B{2}(:,2),B{2}(:,1),'Color',[0 176/255 240/255],'LineWidth',2)
plot(CompleteTemplate.MeanMidBrain(:,2),CompleteTemplate.MeanMidBrain(:,1),'Color',[130/255 213/255 247/255],'LineWidth',2)
set(gca,'position',[0 0 1 1],'units','normalized')
axis off tight equal
%}




end

