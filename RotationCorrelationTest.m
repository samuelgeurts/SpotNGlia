close all
for k1 = stacknumbers
%
    CC = imread([FolderPath,'/',stackinfo(k1).stackname,'.tif']);


    [xcom1,ycom1] = sng_CenterOfMassColor(imcomplement(CC),1);
    shift = (fliplr(size(CC(:,:,1))/2)-[xcom1,ycom1])
%    
    
    CC2 = imtranslate(CC,shift,'OutputView','same','FillValues',255);
    %figure;imagesc(CC2)

%    
    %sh = [abs(ceil(shift)),0]
    %ref1 = imref2d(size(CC))    
    %ref3 = imref2d(size(CC)-2*sh,[ref1.XWorldLimits(1)+sh(1),ref1.XWorldLimits(2)-sh(1)],[ref1.YWorldLimits(1)+sh(2),ref1.YWorldLimits(2)-sh(2)])    
    %CC2 = imwarp(CC,affine2d([1 0 0; 0 1 0; shift(1) shift(2) 1]),'FillValues',0,'Interp','cubic','OutputView',ref3,'SmoothEdges',true);
   
    %[s1,s2,~] = size(CC);
    %CC2 = imcrop(CC, [abs(shift),-abs(shift)] - [shift, shift] + [0 0 s2,s1]);
    %figure;imagesc(CC2)
   
%    
    cm = sng_MaxCircleMask(CC2);
    
    
    CC3=CC2;
    CC3(~cm)=0;
    
    %figure;imagesc(CC)
    %figure;imagesc(CC2)
    figure;imagesc(CC3)

    k1
end
sng_imageslider3

%% Test the performance of the scale-flip-translate correlation method


[Reg] = {stackinfo.Registration}



ref_temp = imref2d(CompleteTemplate.Size);

for k4 = 1:numel(stackinfo)

CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k4).stackname,'.tif']);           
Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k4).ExtendedDeptOfField.IndexMatrix);
%Icombined = imread([ExtendedDeptOfFieldPath,'/',stackinfo(k1).stackname,'-ExtendedDeptOfField.tif']);
flip = stackinfo(k4).Registration(16).value
tform_1 = stackinfo(k4).Registration(17).value;
tform_2 = stackinfo(k4).Registration(18).value;
tform_3 = stackinfo(k4).Registration(19).value;
tform_4 = stackinfo(k4).Registration(20).value;

tform_12 = affine2d(tform_1.T*tform_2.T);
tform_1234 = stackinfo(k4).Registration(21).value;
    
Ialligned = imwarp(Icombined,tform_12,'FillValues',255,'OutputView',ref_temp);
figure;imagesc(uint8(Ialligned))

end

tform_1.T
tform_2.T
tform_1.T*tform_2.T


