function [Ialligned,rgoutput] = AllignmentLink4(Icombined,CompleteTemplate,zfinput)
% this function is a subfunction of "Main Function" which alligns a fish
% with a given template in a couple of steps.
% The function is strongly based on "HorizontalFish7"
%
% Version AllignmentLink2 
%       transforms fishes to template image
%       20170613 added crop before final translation
% Version AllignmentLink4
% make suitable for zfinput and rfoutput

if ~exist('zfinput','var');
    zfinput = struct;
    zfinput = sng_zfinput(zfinput,0,'Registration','BackgroundRemoval','Method','TriangleSmooth','?');
    zfinput = sng_zfinput(zfinput,0,'Registration','BackgroundRemoval','Smooth',1,'?');  
    zfinput = sng_zfinput(zfinput,0,'Registration','BackgroundRemoval','ChannelMethod','cuboid','?');  
    zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','AngleSteps1',100,'?');  
    zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','Scale1',1/16,'?');  
    zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','AngleRange2',0.01,'?');  
    zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','AngleSteps2',50,'?');
    zfinput = sng_zfinput(zfinput,0,'Registration','RotationalAlignment','FinalCenteredYCrop',1000,'?');
    
    zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','CropSize',[1000,1360],'?');      
    zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','ScaleRange',[0.6,1.3],'?');      
    zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','ScaleSteps',200,'?');      
    zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','Scale2',1/16,'?');      
end

rg = zfinput(strcmp({zfinput.stage},'Registration'));
Method_rg = rg(strcmp({rg.name},'Method')).value;
Smooth_rg = rg(strcmp({rg.name},'Smooth')).value;
ChannelMethod_rg = rg(strcmp({rg.name},'ChannelMethod')).value;
AngleSteps1_rg = rg(strcmp({rg.name},'AngleSteps1')).value;
Scale1_rg = rg(strcmp({rg.name},'Scale1')).value;
AngleRange2_rg = rg(strcmp({rg.name},'AngleRange2')).value;
AngleSteps2_rg = rg(strcmp({rg.name},'AngleSteps2')).value;
FinalCenteredYCrop_rg = rg(strcmp({rg.name},'FinalCenteredYCrop')).value;


CropSize_rg = rg(strcmp({rg.name},'CropSize')).value;
ScaleRange_rg = rg(strcmp({rg.name},'ScaleRange')).value;
ScaleSteps_rg = rg(strcmp({rg.name},'ScaleSteps')).value;
Scale2_rg = rg(strcmp({rg.name},'Scale2')).value;

template = CompleteTemplate.Template;

I10 = Icombined;
        
%{
[xcom1,ycom1] = sng_CenterOfMassColor(imcomplement(I10),1);
shift = (fliplr(size(I10(:,:,1))/2)-[xcom1,ycom1]);
I11= imtranslate(I10,shift,'OutputView','same');
figure;imagesc(I11)
%}

%% horizontal allignment
%remove background to compute center of mass in horizontalfish well
I20 = sng_RemoveBackgroundColor3(I10,Method_rg,Smooth_rg,ChannelMethod_rg);

%tform1 = sng_HorizontalFish3(I20,'precise'); %translation + rotation
[tform1,RotationCorrelationOutput,I30] = sng_HorizontalFish5(I20,AngleSteps1_rg,Scale1_rg,AngleRange2_rg,AngleSteps2_rg,FinalCenteredYCrop_rg);



[I30,ref30] = imwarp(I10,tform1,'FillValues',255);
    ref31 = imref2d([1000,size(I30,2)],ref30.XWorldLimits,[-499.5,500.5])
[I32,ref32] = imwarp(I10,tform1,'FillValues',255,'OutputView',ref31);


tform1A = affine2d([1,0,0;0,1,0;0,ref30.YWorldLimits(1)-ref32.YWorldLimits(1),1]);
tform1B = affine2d([1,0,0;0,1,0;tform1.T(3,1),tform1.T(3,2),1])

tformq = affine2d(tform1.T*tform1A.T*tform1B.T)


[I33,ref33] = imwarp(I10,affine2d(tform1.T*tform1A.T*tform1B.T),'FillValues',255);




figure;imagesc(I30)
figure;imagesc(ref30.XWorldLimits,ref30.YWorldLimits,uint8(I30));

figure;imagesc(I32)
figure;imagesc(ref32.XWorldLimits,ref32.YWorldLimits,uint8(I32));

figure;imagesc(I33)
figure;imagesc(ref33.XWorldLimits,ref33.YWorldLimits,uint8(I33));







%% scale and flip 
%scale and leftrigh fitting with template
%cropping is needed to perform well for ScaleFish3 as it is assumed that fishes are centered.
%this crop is not interesting as the image I50 after scaling does not use the cropped image
[I40,ref40] = sng_fishcrop3(I30,ref30,CropSize_rg);
%crop1 = [ref40.XWorldLimits-ref30.XWorldLimits;ref40.YWorldLimits-ref30.YWorldLimits]

tform1b = affine2d([1,0,0;0,1,0;-ref30.XWorldLimits(1),-ref30.YWorldLimits(1),1])
%[I35,ref35] = imwarp(I10,affine2d(tform1.T*tform1b.T),'FillValues',255);
%figure;imagesc(I40)

[tform2,ScaleCorrelationOutput,I45] = sng_ScaleFish4(double(I40),double(template),ScaleRange_rg,ScaleSteps_rg,Scale2_rg);

    tform_12 = affine2d(tform1.T*tform1b.T*tform2.T);
I50 = imwarp(I10,tform_12,'FillValues',255,'OutputView',imref2d(size(template(:,:,1))));
   
    
[I50,ref50] = imwarp(I10,tform_12,'FillValues',255);

%% template matching
%[I51,ref51] = sng_fishcrop3(I50,ref50,[1000,1360]); %crop for better processing in iat algorithm
%tform3 = affine2d([1,0,0;0,1,0;-ref51.XWorldLimits(1),-ref51.YWorldLimits(1),1]); %transform back to pixel coordinates
tform3 = affine2d([1,0,0;0,1,0;-ref50.XWorldLimits(1),-ref50.YWorldLimits(1),1]); %transform back to pixel coordinates



[tform4,I65,CorCoef] = sng_AllignFish2Template2(template,I51,1/4,'translation',3,40);

%CorCoef= sng_NCCoi(uint8(I65),template);


tform_1234 = affine2d(tform1.T*tform2.T*tform3.T*tform4.T);

Ialligned = imwarp(I10,tform_1234,'FillValues',255,'OutputView',imref2d(size(template(:,:,1))));


%{
figure;imagesc(I10)
figure;imagesc(ref30.XWorldLimits,ref30.YWorldLimits,uint8(I30));
figure;imagesc(ref40.XWorldLimits,ref40.YWorldLimits,uint8(I40));
figure;imshowpair(template,uint8(I45))
figure;imshowpair(template,uint8(I50))

figure;imagesc(uint8(I50));
figure;imagesc(ref50.XWorldLimits,ref50.YWorldLimits,uint8(I50));
figure;imagesc(uint8(I51));
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
figure;imshowpair(template,I51)

%}

rgoutput = struct('stage',[],'substage',[],'name',[]','value',[]);

for k11=1:numel(RotationCorrelationOutput)
    rgoutput(k11).stage = 'Registration';
    rgoutput(k11).substage = 'Rotation Correlation';
    rgoutput(k11).name = RotationCorrelationOutput(k11).name;
    rgoutput(k11).value = RotationCorrelationOutput(k11).value;
end

rgoutput = sng_StructFill(rgoutput,{'Registration','Rotation Correlation','tform1_hor',tform1});
rgoutput = sng_StructFill(rgoutput,{'Registration','scale and flip Corr','tform2_scale',tform2});
rgoutput = sng_StructFill(rgoutput,{'Registration','template matching','tform3_trans2pix',tform3});
rgoutput = sng_StructFill(rgoutput,{'Registration','template matching','tform4_trans',tform4});
rgoutput = sng_StructFill(rgoutput,{'Registration','final','tform_complete',tform_1234});

rgoutput = sng_StructFill(rgoutput,{'Registration','Rotation Correlation','ref30',ref30});
rgoutput = sng_StructFill(rgoutput,{'Registration','scale and flip Corr','ref40',ref40});
rgoutput = sng_StructFill(rgoutput,{'Registration','scale and flip Corr','ref50',ref50});

rgoutput = sng_StructFill(rgoutput,{'Registration','template matching','TemplateCorrelation',CorCoef});







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

%}