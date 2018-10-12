function [CorrectedSlice,ppoutput,ImageSliceCor,FiltIm] = PreprocessionSNG(obj,ImageSlice)
%applys rgb and stack correction i.e. translate channels and slices using iat
%add different varaiable structure

% if ~exist('zfinput','var');
%     zfinput = struct
%     zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','onoff',true,'');    
%     zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','sigmalp',1,'moderate'); %sigma for bandpassfilter (lowpass)
%     zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','sigmahp',4,'moderate'); %for bandpassfilter (highpass)
%     zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','scaleC',1/4,'low'); %lowered initial scale to increase computation speed for correlation
%     zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','levelsC',2,'low'); %number of scaled levels correlation is performed
%     zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','iterationsC',10,'low'); %number of iterations iat is performed
%     zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','scaleS',1/4,'low'); %lowered initial scale to increase computation speed for correlation
%     zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','levelsS',3,'low'); %number of scaled levels correlation is performed
%     zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','iterationsS',20,'low'); %iterations iat is performed
% end
% sng_zfinputAssign(zfinput,'preprocession')

obj.SngInputParameters.assign('Preprocession')

n1 = numel(ImageSlice);
NCC_Img_before = cell(1,n1);
NCC_Img_after = cell(1,n1);
NCC_BP_before = cell(1,n1);
NCC_BP_after = cell(1,n1);
ECCWarp = cell(1,n1);
ImageSliceCor = cell(1,n1);
crop = cell(1,n1);

if nargout >= 4
    FiltIm = cell(1,n1);
end

%% does RGB correction on all images of one fish  
for k4=1:n1

    if onoff 
        %bandpass filter focussed on spotsize passing
        FilteredImage = sng_RGBselection(ImageSlice{k4}(:,:,1:3),sigmalp,sigmahp);
    else
        %no bandpass filtering
        FilteredImage = ImageSlice{k4}(:,:,1:3);
    end
    
    if nargout >= 4
        FiltIm{k4} = FilteredImage;
    end
    
    %{
    figure;imagesc((FilteredImage(:,:,1:3)));axis equal off tight
    figure;imagesc(uint8(ImageSlice{1}(:,:,1:3)));axis equal off tight
    %}
    %calculates the needed warp
    [ECCWarp{k4}, ~]= sng_RGB2(FilteredImage(:,:,1:3),scaleC,'translation',levelsC,iterationsC);
    %warp images
    [ImageSliceCor{k4},crop{k4}] = sng_RGB_IATwarp2(ImageSlice{k4}(:,:,1:3),ECCWarp{k4});
    %imwrite(uint8(Img2), fullFileName{l}, 'WriteMode', 'append', 'Compression','none');

    %normal correlation
    %{ 
    CC1 = corr2(FilteredImage(:,:,1),FilteredImage(:,:,2));
    CC2 = corr2(FilteredImage2(:,:,1),FilteredImage2(:,:,2));
    CC3 = corr2(ImageSlice{k4}(:,:,1),ImageSlice{k4}(:,:,3));
    CC4 = corr2(ImageSliceCor{k4}(:,:,1),ImageSliceCor{k4}(:,:,3)); 
    %}
    [FilteredImage2,~] = sng_RGB_IATwarp2(FilteredImage,ECCWarp{k4});

    %original image before and after normalised correlation (green-red and blue-red channel)
    NCC_Img_before{k4}(1) = sng_NCC(ImageSlice{k4}(:,:,1),ImageSlice{k4}(:,:,2));
    NCC_Img_before{k4}(2) = sng_NCC(ImageSlice{k4}(:,:,1),ImageSlice{k4}(:,:,3));
    
    NCC_Img_after{k4}(1) = sng_NCC(ImageSliceCor{k4}(:,:,1),ImageSliceCor{k4}(:,:,2));    
    NCC_Img_after{k4}(2) = sng_NCC(ImageSliceCor{k4}(:,:,1),ImageSliceCor{k4}(:,:,3));
    
    %bandpass filtered before and after normalised correlation  (green-red and blue-red channel)
    NCC_BP_before{k4}(1) = sng_NCC(FilteredImage(:,:,1),FilteredImage(:,:,2));
    NCC_BP_before{k4}(2) = sng_NCC(FilteredImage(:,:,1),FilteredImage(:,:,3));
    
    NCC_BP_after{k4}(1) = sng_NCC(FilteredImage2(:,:,1),FilteredImage2(:,:,2));
    NCC_BP_after{k4}(2) = sng_NCC(FilteredImage2(:,:,1),FilteredImage2(:,:,3));
    
end

%% allign different field dept images of one fish and save it as fullFileName1
[ECCWarp2,CorrectedSlice] = sng_AllignFishBatch3(ImageSliceCor,scaleS,'translation',levelsS,iterationsS)  

ppoutput.ColorWarp = ECCWarp;
ppoutput.SliceWarp = ECCWarp2;
ppoutput.NCC_Img_before = NCC_Img_before;
ppoutput.NCC_Img_after = NCC_Img_after;
ppoutput.NCC_BP_before = NCC_BP_before;
ppoutput.NCC_BP_after = NCC_BP_after;


%% RGB Correction figures and movies

%{

figure;imagesc(uint8(ImageSlice{1}(:,:,1:3)));axis equal off tight
set(gca,'position',[0 0 1 1],'units','normalized')


writerObj1 = VideoWriter('zoom in initial slice','MPEG-4');
writerObj1.FrameRate = 25;    % to perform realtime movie
open(writerObj1);
figure;imagesc(uint8(ImageSlice{1}(:,:,1:3)));axis equal off tight
set(gca,'position',[0 0 1 1],'units','normalized')
h=gca
siz = size((ImageSlice{1}(:,:,1:3)));
zoomn = [1,20]
zoompoint1 = siz(1:2)/2
zoompoint2 = [459,768]
step = 60
[zoomp,yzp,xzp] = sng_smoothzoom(zoomn,zoompoint1,zoompoint2,step,h)
for k = 1:step    
    sng_zoom(zoomp(k),[yzp(k),xzp(k)],siz(1:2),h)
    frame = getframe(gcf);
    writeVideo(writerObj1,frame);   
end
close(writerObj1);

siz = size((ImageSliceCor{1}(:,:,1:3)));
zoomn = [20,1]
zoompoint1 = [459,768] - crop(1:2)
zoompoint2 = siz(1:2)/2
step = 40

figure;imagesc(uint8(ImageSliceCor{1}));axis equal off tight
sng_zoom(zoomn(1),zoompoint1,siz(1:2),gca)
set(gca,'position',[0 0 1 1],'units','normalized')


writerObj1 = VideoWriter('zoom out corrected slice','MPEG-4');
writerObj1.FrameRate = 25;    % to perform realtime movie
open(writerObj1);
figure;imagesc(uint8(ImageSliceCor{1}));axis equal off tight
set(gca,'position',[0 0 1 1],'units','normalized')
h=gca
[zoomp,yzp,xzp] = sng_smoothzoom(zoomn,zoompoint1,zoompoint2,step,h)
for k = 1:step    
    sng_zoom(zoomp(k),[yzp(k),xzp(k)],siz(1:2),h);
    frame = getframe(gcf);
    writeVideo(writerObj1,frame);   
end
close(writerObj1);
%}
    
    
%% Stack Correction figures and movies
%{

%frame parameters
a = numel(ImageSlice) %scroll
up = 1:a
down = linspace(a,1,a)
freeze = ones(1,numel(ImageSlice)) %freeze
k = [up,a*freeze,down,freeze]
k = repmat(k,1,4)

writerObj1 = VideoWriter('stack scroll uncor','MPEG-4');
writerObj1.FrameRate = 6;    % to perform realtime movie
open(writerObj1);
figure;imagesc(uint8(ImageSliceCor{1}));axis equal off tight
set(gca,'position',[0 0 1 1],'units','normalized')
sz = size(ImageSliceCor{1})

for j = k    
    imagesc(uint8(ImageSliceCor{j}));axis equal off tight 
    t = text(sz(2)*0.95,sz(1)*0.05,num2str(j),'FontSize',20)
    frame = getframe(gcf);
    writeVideo(writerObj1,frame);   
end
close(writerObj1);



%frame parameters
a = numel(ImageSlice) %scroll
up = 1:a
down = linspace(a,1,a)
freeze = ones(1,numel(ImageSlice)) %freeze
k = [up,a*freeze,down,freeze]
k = repmat(k,1,4)

writerObj1 = VideoWriter('stack scroll cor','MPEG-4');
writerObj1.FrameRate = 6;    % to perform realtime movie
open(writerObj1);
figure;imagesc(uint8(CorrectedSlice{1}));axis equal off tight
set(gca,'position',[0 0 1 1],'units','normalized')

for j = k    
    imagesc(uint8(CorrectedSlice{j}));axis equal off tight
    t = text(sz(2)*0.95,sz(1)*0.05,num2str(j),'FontSize',20)
    frame = getframe(gcf);
    writeVideo(writerObj1,frame);   
end
close(writerObj1);
%}




%{

figure;imshow(FilteredImage)
figure;imshow(ImagesSliceCor{1})
figure;imshow(Images2)
figure;imshow(Icombined)

%}





