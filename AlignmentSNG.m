function [Ialligned,rgoutput] = AllignmentLink5(Icombined,CompleteTemplate,zfinput)
% this function is a subfunction of "Main Function" which alligns a fish
% with a given template in a couple of steps.
% The function is strongly based on "HorizontalFish7"
%
% Version AllignmentLink2 
%       transforms fishes to template image
%       20170613 added crop before final translation
% Version AllignmentLink4
%       make suitable for zfinput and rfoutput
% Version AllignmentLink5
%       put cropping into horizontal fish, so cropsize can be removed

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
    
    %zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','CropSize',[1000,1360],'?');      
    zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','ScaleRange',[0.6,1.3],'?');      
    zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','ScaleSteps',200,'?');      
    zfinput = sng_zfinput(zfinput,0,'Registration','ScaleFlipAlignment','Scale2',1/16,'?');

    zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','RemoveTail',600,'');   
    zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','Scale3',1/4,''); %lowered initial scale to increase computation speed for correlation
    zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','AffineMethod','translation',''); %affine or translation
    zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','Levels',3,''); %number of scaled levels correlation is performed
    zfinput = sng_zfinput(zfinput,0,'Registration','FinalTemplateMatching','Iterations',40,''); %number of iterations iat is performed
end

sng_zfinputAssign(zfinput,'Registration')

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
[I20,BackgroundMask,BackgroundThreshold] = sng_RemoveBackgroundColor3(I10,Method,Smooth,ChannelMethod);

%stretchlim parameters (not used for registration so far)
rtemp=I20(:,:,1);gtemp=I20(:,:,2);btemp=I20(:,:,3);
stretchlim_1procent(1:2,1) = stretchlim(rtemp(~BackgroundMask),0.01);    
stretchlim_1procent(1:2,2) = stretchlim(gtemp(~BackgroundMask),0.01);
stretchlim_1procent(1:2,3) = stretchlim(btemp(~BackgroundMask),0.01);
stretchlim_0procent(1:2,1) = stretchlim(rtemp(~BackgroundMask),0);    
stretchlim_0procent(1:2,2) = stretchlim(gtemp(~BackgroundMask),0);
stretchlim_0procent(1:2,3) = stretchlim(btemp(~BackgroundMask),0);


%tform1 = sng_HorizontalFish3(I20,'precise'); %translation + rotation
[tform1,RotationCorrelationOutput,I30,ref30] = sng_HorizontalFish5(I20,AngleSteps1,Scale1,AngleRange2,AngleSteps2,FinalCenteredYCrop);


%{
figure;imagesc(I10)
figure;imagesc(I20)
figure;imagesc(I30)
figure;imagesc(ref30.XWorldLimits,ref30.YWorldLimits,uint8(I30));
%}

%% scale and flip 
%scale and leftrigh fitting with template
%input is needs to be a cropped image for good performance


[tform2,ScaleCorrelationOutput,I50] = sng_ScaleFish4(double(I30),double(template),ScaleRange,ScaleSteps,Scale2);

%{
figure;imagesc(uint8(I50));
    tform_12 = affine2d(tform1.T*tform2.T);
[I50,ref50] = imwarp(I10,tform_12,'FillValues',255,'OutputView',imref2d(size(template(:,:,1))));
   
figure;imagesc(uint8(I50));
figure;imagesc(ref50.XWorldLimits,ref50.YWorldLimits,uint8(I50));

figure;imagesc(ScaleCorrelationOutput(5).value)


%}

%% template matching
%the tail is cut of as it will inprove the template matching, the template has an not well difined tail
tform3 = affine2d([1 0 0; 0 1 0; -RemoveTail 0 1]);
I55 = imwarp(I50,tform3,'FillValues',0,'Interp','cubic','OutputView',imref2d(size(I50(:,:,1))-[0 600]));

%final subpixel template matching using 3th party iat software

%AffineMethod = 'affine'
%AffineMethod = 'translation'

Initialization = [-RemoveTail * Scale3; 0 ];
if strcmp(AffineMethod,'affine')
    Initialization = [eye(2),Initialization];
end


[tform4,~,CorCoef] = sng_AllignFish2Template2(template,I55,Scale3,AffineMethod,Levels,Iterations,Initialization);%missingoutput=I60


%CorCoef= sng_NCCoi(uint8(I65),template);


tform_1234 = affine2d(tform1.T*tform2.T*tform3.T*tform4.T);

Ialligned = imwarp(I10,tform_1234,'FillValues',255,'OutputView',imref2d(size(template(:,:,1))));
%% extra mean fish color parameter
%get image - take front(head) - compute histogram - smooth - find second max peak (no background)
for k3 = 1:3
    Img = Ialligned(:,700:end,k3);                                        
    MeanFishColor(k3) = mean(Img(:) < BackgroundThreshold(k3));
    h = hist(Img(:),0:1:255);
    %h2 = smoothdata(h,'gaussian',6);
    h2 = imgaussfilt(h,1); %gaussian filter function for matlab older than 2017a
    [~,MaxFishColor(k3)] = max(h2(1:floor(BackgroundThreshold(k3))));    
end





%{
figure;imagesc(I10)
figure;imagesc(uint8(I55));
figure;imagesc(uint8(I60));
figure;imagesc(uint8(template));
figure;imagesc(uint8(CompleteTemplate.MidBrainDistanceMap));

figure;imagesc(uint8(Ialligned));colormap(gray);
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

for k12 = 1:numel(ScaleCorrelationOutput)
    rgoutput(k12+k11).stage = 'Registration';
    rgoutput(k12+k11).substage = 'Scale Correlation';
    rgoutput(k12+k11).name = ScaleCorrelationOutput(k12).name;
    rgoutput(k12+k11).value = ScaleCorrelationOutput(k12).value;
end

rgoutput = sng_StructFill(rgoutput,{'Registration','Rotation Correlation','tform1_hor',tform1});
rgoutput = sng_StructFill(rgoutput,{'Registration','scale and flip Corr','tform2_scale',tform2});
rgoutput = sng_StructFill(rgoutput,{'Registration','template matching','tform3_tailcut',tform3});
rgoutput = sng_StructFill(rgoutput,{'Registration','template matching','tform4_trans',tform4});
rgoutput = sng_StructFill(rgoutput,{'Registration','final','tform_complete',tform_1234});

%rgoutput = sng_StructFill(rgoutput,{'Registration','Rotation Correlation','ref30',ref30});
%rgoutput = sng_StructFill(rgoutput,{'Registration','scale and flip Corr','ref40',ref40});
%rgoutput = sng_StructFill(rgoutput,{'Registration','scale and flip Corr','ref50',ref50});

rgoutput = sng_StructFill(rgoutput,{'Registration','template matching','TemplateCorrelation',CorCoef});

rgoutput = sng_StructFill(rgoutput,{'Registration','Background Removal','Background',BackgroundMask});
rgoutput = sng_StructFill(rgoutput,{'Registration','Background Removal','BackgroundThreshold',BackgroundThreshold});
rgoutput = sng_StructFill(rgoutput,{'Registration','Background Removal','stretchlim_0procent',stretchlim_0procent});
rgoutput = sng_StructFill(rgoutput,{'Registration','Background Removal','stretchlim_1procent',stretchlim_1procent});
rgoutput = sng_StructFill(rgoutput,{'Registration','Background Removal','MeanFishColor',MeanFishColor});
rgoutput = sng_StructFill(rgoutput,{'Registration','Background Removal','MaxFishColor',MaxFishColor});

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