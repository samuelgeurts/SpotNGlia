function obj = TemplateSpot(obj)
%this function computes the template variables of the test batch with annotations.


CubeSizeSpotErode = 4; %spot structure element for erode
CubeSizeSpotDilate = 4; %background structure element for dilate

winsiz = 20; %size of spot box is 2*winsiz + 1 = 51x51
f = int16(-winsiz:winsiz);    

subplot = @(m,n,p) subtightplot (m, n, p, [0 0], [0 0], [0 0]);
image_TF = true;
Save_TF = true;
HistEq_TF = true;

obj = LoadTemplate(obj);

obj.StackInfo
obj.CompleteTemplate.ref_temp
obj.ZFParameters



load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo')


ColorToGrayVector = CompleteTemplate.SpotContrastVector;
%ColorToGrayVector = [0 1 0]
%ColorToGrayVector = [-0.378629920630393;0.865619921584407;-0.327630179561694]        

%{
zfinput = struct;
zfinput = sng_zfinput(zfinput,0,'SpotDetection','RgbToGray','ColorToGrayVector',ColorToGrayVector,'');  %select color channel [0 1 0], 
zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleLevels',9,'');  
zfinput = sng_zfinput(zfinput,0,'SpotDetection','Wavelet','ScaleBase',0.5,'');      
zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPlevels',4:6,'');  
zfinput = sng_zfinput(zfinput,0,'SpotDetection','MultiProduct','MPthreshold',1,'');
%}

sng_zfinputAssign(obj.ZFParameters,'SpotDetection')


spotcolorsTT = [];
backrcolorsTT = [];


for k1 = 1:numel(obj.StackInfo)
    
    %transformation matrix to aligned fish
    tform_1234 = RegistrationInfo{k1}(strcmp({RegistrationInfo{k1}.name},'tform_complete')).value;
    %open correctedfish
    CorrectedSlice = sng_openimstack2([obj.SavePath, '/', 'CorrectedFish', '/', obj.StackInfo(k1).stackname, '.tif']);
   %annotated spots
    RoiMicroglia = ReadImageJROI([obj.AnnotatedSpotPath,'/',obj.StackInfo(k1).stackname,'.roi']);
    spota = RoiMicroglia.mfCoordinates;
    spots = RoiMicroglia.vnSlices;
    [SpotAnn(:,1),SpotAnn(:,2)] = transformPointsForward(tform_1234,spota(:,1),spota(:,2));

    
    
    %for k2 = 1:4
    % mean(CorrectedSlice{k2}(CorrectedSlice{k2}(:) < bkt(1)))
    %end
    
        %{
        %histogram equalization
        %Wavelet spotdetection doesn work very well when used histogram equalization
        if HistEq_TF
                bkt = RegistrationInfo{k1}(strcmp({RegistrationInfo{k1}.name},'BackgroundThreshold')).value;
                %sl = RegistrationInfo{k1}(strcmp({RegistrationInfo{k1}.name},'stretchlim_0procent')).value;
                %sl(1,1:3) = 0;
            for k2 = 1:numel(CorrectedSlice)
                %histogram stretching               
                %IStrech{k2} = imadjust(CorrectedSlice{1},sl,[]);

                %histogram equalization excluding background
                bp = sum(CorrectedSlice{k2}(:) > bkt(1)) / numel(CorrectedSlice{k2}); %fraction of background pixels   
                hgram = [(1-bp)/254 * ones(1,254),bp]; %flat histogram unless large peak containing background
                IEqual{k2} = histeq(CorrectedSlice{1},hgram);
                %{
                figure;imagesc(CorrectedSlice{k2}(:,:,2))
                figure;imagesc(IStrech{k2}(:,:,2))
                figure;imagesc(IEqual{k2}(:,:,2))
                %}
            end
        else
            IEqual = CorrectedSlice;
        end
        %}
    
                
    
%% frame around annotated spots

    backrcolorsT = [];
    spotcolorsT = [];
    
    
    if image_TF
        subplotvar = ceil(sqrt(2*numel(spots)));
        figure;
    end
    
    
    for l = 1:numel(spots)
        Spot1 = IEqual{spots(l)}(f+spota(l,2),f+spota(l,1),1:3);
        Spot2 = imgaussfilt(Spot1,0.7); 
        GSpot2 = sng_RGB2Gray(Spot2,ColorToGrayVector,false); %gray transformed spot
        Mask2 = sng_SpotWavelet(GSpot2,ScaleLevels,ScaleBase,MPlevels,MPthreshold,false);
        SE1 = strel('cube',CubeSizeSpotErode); %spot structure element for erode
        SE2 = strel('cube',CubeSizeSpotDilate); %background structure element for dilate       
        SpotMask = imerode(Mask2,SE1);
        BackMask = ~imdilate(Mask2,SE2);
        
        SpotF = Spot1;SpotF(repmat(SpotMask,1,1,3) == 0) = 0;
        BackF = Spot1;BackF(repmat(BackMask,1,1,3) == 0) = 0;

        %{
        figure;imagesc(Spot1)
        figure;imagesc(Spot2)
        figure;imagesc(GSpot2)
        figure;imagesc(Mask2)
        figure;imagesc(SpotMask)
        figure;imagesc(BackMask)
        figure;imagesc(SpotF)
        figure;imagesc(BackF)
        %}

        spotcolorsarray = Spot1(repmat(SpotMask,1,1,3));
        spotcolors = reshape(spotcolorsarray,numel(spotcolorsarray)/3,3);

        backcolorsarray = Spot1(repmat(BackMask,1,1,3));
        backcolors = reshape(backrcolorsarray,numel(backrcolorsarray)/3,3);

        %used for previous color measure probability
        %spotcolorsT = [spotcolorsT;spotcolors];
        %backrcolorsT = [backrcolorsT;backrcolors];
        
        
        
        
        if image_TF
            subplot(subplotvar,subplotvar,(2*l)-1);imagesc(SpotF);axis off tight equal       
            subplot(subplotvar,subplotvar,(2*l));imagesc(BackF);axis off tight equal       
        end
        
    end
    toc
drawnow
spotcolorsTT = [spotcolorsTT;spotcolorsT];
backrcolorsTT = [backrcolorsTT;backrcolorsT];

end

S = double(spotcolorsTT);
B = double(backrcolorsTT);

%%
spotcolorsTT2 = unique(spotcolorsTT,'rows');
backrcolorsTT2 = unique(backrcolorsTT,'rows');



%
figure;scatter3(spotcolors(:,1),spotcolors(:,2),spotcolors(:,3));
hold on;scatter3(backrcolors(:,1),backrcolors(:,2),backrcolors(:,3));
axis equal;xlim([1,255]);ylim([1,255]);zlim([1,255]);xlabel('Red');ylabel('Green');zlabel('Blue')

figure;scatter3(spotcolorsTT2(:,1),spotcolorsTT2(:,2),spotcolorsTT2(:,3));
hold on;scatter3(backrcolorsTT2(:,1),backrcolorsTT2(:,2),backrcolorsTT2(:,3));
axis equal;xlim([1,255]);ylim([1,255]);zlim([1,255]);xlabel('Red');ylabel('Green');zlabel('Blue')
drawnow
%

%% feature segmentation.
%{
data = double([backrcolorsTT2;spotcolorsTT2]);
label = [(ones(size(backrcolorsTT2,1),1));2*(ones(size(spotcolorsTT2,1),1))];
A = prdataset(data(:,1:3),label)

W2 =  perlc(A)

figure;scatterd(A(:,1:2))
plotc({W2})

mat = eye(3)
mat(:,1:2) = struct(W2).data.rot

D = A(:,1:2)*W2

figure;scatterd(D)
%}


%% histogram spot/background contrast.


%{
bin = 255;


Shist=sng_chistcount(S',linspace(0,255,bin+1));      %create bins
Bhist=sng_chistcount(B',linspace(0,255,bin+1));      %

Snorm = Shist/size(S,1);
Bnorm = Bhist/size(B,1);

[minC1,C1,Q1] = sng_threshold2(Snorm,Bnorm);                %find threshold

sng_chistplot(Snorm,Bnorm,minC1);set(gcf,'numbertitle','off','name','RGB') 

%% a new way to transform to gray instead of taking the green values
bm = mean(B)
sm = mean(S)


b = (bm-sm)
bnorm = b/norm(b)

mx = dot([255,255,255],bnorm)
mn = dot([0,0,0],bnorm)

st = double(S)*bnorm'; %transformed spot
bt = double(B)*bnorm'; %transformed background

st = st* 255/mx;
bt = bt* 255/mx;

SThist = histcounts(st,linspace(0,255,bin+1));
BThist = histcounts(bt,linspace(0,255,bin+1));

STnorm = SThist/size(S,1);
BTnorm = BThist/size(B,1);

[minC2,C2,Q2] = sng_threshold2(STnorm,BTnorm);                %find threshold

sng_chistplot(STnorm,BTnorm,minC2);set(gcf,'numbertitle','off','name','RGB') 


Overlap = min(STnorm,BTnorm);
Combination = max(STnorm,BTnorm);
Jaccard = 1-sum(Overlap)/sum(Combination)
%}
%%


[cartvec bestJac polarvec index] = FindBestContrastVector(S,B,7,11,false,true)
bestvec = cartvec(:,end);

if image_TF
    Jc1 = sng_Rgb2BwContrast(S,B,bestvec,true);set(gca,'XLim',[75,160])
    Jc2 = sng_Rgb2BwContrast(S,B,[1;0;0],true);
    Jc3 = sng_Rgb2BwContrast(S,B,[0;1;0],true);
    Jc4 = sng_Rgb2BwContrast(S,B,[0;0;1],true);
end


%% Compute BW image based on transform vector
k10=10
CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k10).stackname,'.tif']);           
Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k10).ExtendedDeptOfField.IndexMatrix);

Img2 = sng_RGB2Gray(Icombined, bestvec,true);
%Img2 = sng_RGB2Gray(Icombined,[0;1;0],true);

 
%% create volume that contains which a color vector can compared with to determine if it is in the range of spots

%spot colorvectors scatter plot
%{
figure;scatter3(spotcolorsTT2(:,1),spotcolorsTT2(:,2),spotcolorsTT2(:,3),'filled',1);
axis equal;xlim([1,255]);ylim([1,255]);zlim([1,255]);xlabel('Red');ylabel('Green');zlabel('Blue')

%mask spot colorvectors
SpotMask = false(255,255,255);
for n = 1:size(spotcolorsTT2,1)
    SpotMask(spotcolorsTT2(n,1),spotcolorsTT2(n,2),spotcolorsTT2(n,3)) = true;
end
figure;hi = imagesc(SpotMask(:,:,1))
for m=1:255
    set(hi,'CData',SpotMask(:,:,m))
    drawnow
end
%}

%mask spot colorvectors with closed holes. Doesnt work for to less data, i.e. too
%much single pixels
%{
SpotMask2=imfill(SpotMask,'holes');
figure;hi = imagesc(SpotMask2(:,:,1))
for m=1:255
    set(hi,'CData',SpotMask2(:,:,m))
    drawnow
end
%}

%counts how much from each color vector, i.e. generates a valued Mask
%correction from 255 to 256 and added +1
%the spotdetectionlink function has also to be changed
SpotColorImg = zeros(256,256,256);
for n2 = 1:size(S,1)
    SpotColorImg(S(n2,1)+1,S(n2,2)+1,S(n2,3)+1) = SpotColorImg(S(n2,1)+1,S(n2,2)+1,S(n2,3)+1) + 1;
end    
%figure;hist(SpotImg(SpotImg ~= 0),1:8)

SpotMaskGauss = imgaussfilt3(double(SpotColorImg),5);


SVAP_index = find(SpotMaskGauss ~= 0);
%[subscript(:,1),subscript(:,2),subscript(:,3)] = (ind2sub([255,255,255],SVAP_index));
%SVAP_subscript = uint8(subscript);
SpotVectorArrayProbability = SpotMaskGauss(SVAP_index);


%{
writerObj1 = VideoWriter('blurred spots','MPEG-4');
writerObj1.FrameRate = 50;    % to perform realtime movie
open(writerObj1);
figure;hi = imagesc(SpotMaskGauss(:,:,1),[0 max(SpotColorImg(:))]);colormap(jet(1000));sng_imfix(2)
for m = 1:255
    set(hi,'CData',SpotMaskGauss(:,:,m))
    drawnow
    frame = getframe(gcf);
    writeVideo(writerObj1,frame)
end
close(writerObj1);    
%}

%sum(SpotMask(:))
%sum(SpotMask2(:))

%Regionstats = regionprops(SpotMask) %to much single pixels
%Regionstats = regionprops(boolean(SpotMaskGauss)) %to much single pixels



%{
[X,Y,Z] = meshgrid(1:255,1:255,1:255);
S2 = double(spotcolorsTT2);
kk = boundary(S2(:,1),S2(:,2),S2(:,3));
figure;trisurf(kk,S2(:,1),S2(:,2),S2(:,3),'Facecolor','red','FaceAlpha',0.1)


DT = delaunayTriangulation(S2(:,1),S2(:,2),S2(:,3));
K = convexHull(DT)
figure;trisurf(K,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3),'Facecolor','cyan','FaceAlpha',0.1)
%}



SpotTemplateVar.SpotContrastVector = bestvec;
SpotTemplateVar.SpotVectorArrayProbability = SpotVectorArrayProbability;
SpotTemplateVar.SVAP_index = SVAP_index;

if Save_TF
    save([Basepath,'/Template Spot/SpotTemplate.mat'],'SpotTemplateVar');
    save([TemplatePath3dpf,'/SpotTemplate.mat'],'SpotTemplateVar');
end

%% strech image the optain lower variance of spot color
%{



%}

end