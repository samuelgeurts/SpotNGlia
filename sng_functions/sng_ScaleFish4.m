function [tform_2,ScaleCorrelationOutput,I50] = sng_ScaleFish4(I40,template,ScaleRange_rg,ScaleSteps_rg,Scale2_rg)
%this function persforms correlation between images and differnt scales of
%the template. The maxima in the correlationmatrix is determined and the
%scale selected.
% Version sng_ScaleFish4 usable for better control of input and output

%{ 
I40=double(I30);
%CropSize_rg = [1000,1360];
ScaleRange_rg = [0.6,1.3];
ScaleSteps_rg = 200;
Scale2_rg = 1/16;
%}

scales = linspace(ScaleRange_rg(1),ScaleRange_rg(2),ScaleSteps_rg);

Img2 = rgb2gray(uint8(I40));
template2 = rgb2gray(uint8(template));

[template_sc,Img_sc] = sng_scale4(Scale2_rg,template2,Img2);

%minimal overlap number of pixels for cross correlation for speeding up
mins= min(size(template_sc),size(Img_sc));
C = mins(1) * mins(2)/6;

%preallocation
maxk = zeros(2,numel(scales));
ck = zeros(2,numel(scales));
siz = zeros(2,numel(scales));

%when only the calculation scale as output

for j=1:numel(scales)  
    %templatevar = sng_scale2(template,scales(j));
    Imgvar = sng_scale4(scales(j),Img_sc);
    I100 = normxcorr2_general(Imgvar,template_sc,C);
    I100flip = normxcorr2_general(rot90(Imgvar,2),template_sc,C);
    %[maxk(1,j),ck(1,j)] = max(max(I100))
    %[maxk(2,j),ck(2,j)] = max(max(I100flip))
    [maxk(1,j),ck(1,j)] = max(I100(:));
    [maxk(2,j),ck(2,j)] = max(I100flip(:));
    siz(:,j) = size(I100);
end

%find de maxima of all correlation image maxima
%csflip is 1 for no flip and 2 for a flipped template
%csj could be 1 to numel(scales) dependent on the best match
[~, cs] = max(maxk(:));                %cs is de index with the highest correlation
[csflip,csj] = ind2sub(size(maxk),cs); 
scalebest = scales(csj);


%{
figure;plot(1:200,maxk(1,:))
hold on; plot(1:200,maxk(2,:))

figure;imagesc(Imgvar)
figure;imagesc(I100)
figure;imagesc(I100flip)    
    Imgvar = sng_scale4(scalebest,Img_sc);
    if csflip == 1; I101 = normxcorr2_general(Imgvar,template_sc,C);end;
    if csflip == 2; I101 = normxcorr2_general(rot90(Imgvar,2),template_sc,C);end;
    figure;imagesc(I101)
%}

%recalculates the best template match and correlation image

    I40(I40(:,:,1)==255) = 100;
    [sy,sx,~] = size(I40);
    
    Imgvar = sng_scale4(scalebest,Img_sc);
    if csflip == 1           
        Tflip = [1 0 0;0 1 0; 0 0 1];
        Corrbest = normxcorr2_general(Imgvar,template_sc,C);
    elseif csflip == 2 
        Tflip = [-1 0 0;0 -1 0; sx,sy 1];        
        Corrbest = normxcorr2_general(rot90(Imgvar,2),template_sc,C);
    end

%peak coordinates in correlationimage
[y,x] = ind2sub(siz(:,csj),ck(csflip,csj));
offset = [y,x] - size(Imgvar);
scaledoffset = offset / Scale2_rg;

Ttranslate = [1,0,0;0,1,0;scaledoffset(2),scaledoffset(1),1]; %transform back to pixel coordinates
Tscale = [scales(csj) 0 0; 0 scales(csj) 0;0 0 1];    
tform_2 = affine2d(Tflip*Tscale*Ttranslate);

I50 = imwarp(I40,tform_2,'FillValues',255,'OutputView',imref2d(size(template(:,:,1))));

%{
%to test performance indiviual transforms
I41 = imwarp(I40,affine2d(Tscale));figure;imagesc(uint8(I41))
I42 = imwarp(I41,affine2d(Tflip));figure;imagesc(uint8(I42))
I43 = imwarp(I42,affine2d(Ttranslate),'OutputView',imref2d(size(template(:,:,1))));figure;imagesc(uint8(I43))

Tscale = eye(3)
Ttranslate = eye(3)
Tflip = eye(3)
%}





%{
figure;imagesc(Corrbest)

figure;imagesc(I40)
figure;imagesc(uint8(I50))
figure;imshowpair(template,uint8(I50))

figure; surf(Corrbest)
%}

ScaleCorrelationOutput(1).name = 'MaxCorrelationValue';
ScaleCorrelationOutput(2).name = 'Tflip';
ScaleCorrelationOutput(3).name = 'Tscale';
ScaleCorrelationOutput(4).name = 'Ttranslate';
ScaleCorrelationOutput(5).name = 'BestScaleCorrelationImage';
ScaleCorrelationOutput(6).name = 'MaxCorrelationPlot';
ScaleCorrelationOutput(7).name = 'FlipIsTwo';


ScaleCorrelationOutput(1).value = scalebest;
ScaleCorrelationOutput(2).value = Tflip;
ScaleCorrelationOutput(3).value = Tscale;
ScaleCorrelationOutput(4).value = Ttranslate;
ScaleCorrelationOutput(5).value = Corrbest;
ScaleCorrelationOutput(6).value = maxk;
ScaleCorrelationOutput(7).value = csflip;

end

