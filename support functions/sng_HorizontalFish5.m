function [tform_rc2,RotationCorrelationOutput,varargout] = sng_HorizontalFish5(Img,varargin)
%this function rotates a fish so it is orientated horizontal
%the code is in 2 main steps, where first a global rotation angel is
%calculateds and finetuned afterwards
% 
% version sng_HorizontalFish3
%   added figures in code between brace   
%   07-03-2017 add image and rotation movie between brace
% version sng_HorizontalFish4 20170616
%   change for zfinput wih direct variables
% version sng_HorizontalFish5
%   changed tform output, add transformation to pixel coordinates with
%   centered image over the y-axis with size FinalCenteredYCrop


%{
    Img = I20;
    AngleSteps1 = 100;
    Scale1 = 1/16;    

    AngleSteps2 = 50;
    AngleRange2 = 0.02;
    FinalCenteredYCrop = 1000;
%}

%Example
%    %'fast'
%    angles1 = linspace(2*pi/50,2*pi,50);
%    scale1 = 1/32;

%    
%    %'precise'
%    angles1 = linspace(2*pi/100,2*pi,100);
%    scale1 = 1/16;
%    +Imcrop parameters
%    angles2 = linspace(0.99*pi,1.01*pi,50)
%
%    [tform,Img6,ref6] = sng_HorizontalFish4(Img,AngleSteps1,Scale1,AngleSteps2,AngleRange2)
 


AngleSteps1 = varargin{1};
Scale1 = varargin{2};
AngleRange2 = varargin{3};
AngleSteps2 = varargin{4};
FinalCenteredYCrop = varargin{5};

angles1 = linspace(2*pi/AngleSteps1,2*pi,AngleSteps1);
    

%% first iteration global angle calculation
Imgc = imcomplement(uint8(Img));
%calculate center of mass.
[xcom1,ycom1] = sng_CenterOfMassColor(Imgc,1);
%scale image for global calculation
[Img2,xcom2,ycom2] = sng_scale(Imgc,xcom1,ycom1,Scale1);
%apply normalized rotation correlation
[theta1,maxangle1,peaksep1,ncorr1] = sng_NormalizedRotationCorrelation3(Img2,angles1,[xcom2,ycom2]);

%{
figure;imagesc(Img)
figure;imagesc(Imgc)
figure;imagesc(Img2)
[theta1,maxangledeg1,peaksep1,ncorr1,I40,I50,I61] = sng_NormalizedRotationCorrelation3(Img2,angles1,[xcom2,ycom2],1:10:100)
[theta,angles,peaksep,ncorr,I40,I50,I61] = sng_NormalizedRotationCorrelation2(Img2,rad2deg(angles1),[xcom2,ycom2],1:50);
figure;imagesc(I40)
figure;imagesc(I50)
for k = 1:numel(I61)
    figure;imagesc(uint8(I61{k}))
end
figure;plot(angles1,ncorr);xlim([0,2*pi])

%}

%rotate
Ttrans= [1,0,0;0,1,0;-xcom1,-ycom1,1];
Trotate = [cos(theta1),-sin(theta1),0;sin(theta1),cos(theta1),0;0,0,1]; 
tform_rc1 = affine2d(Ttrans*Trotate);
%tform = affine2d(Trotate);

%% a more precise rotation calculation

if ~isempty(varargin{3})

    angles2 = pi + 2*pi*linspace(- AngleRange2,AngleRange2,AngleSteps2);
    
    %!!! misschien crop veranderen naar een ingebouwde imwarp crop
    %I60 = imwarp(I50,tform,'bilinear','OutputView',Rout, 'SmoothEdges',true,'FillValues',0);
       
    [Img3,ref3] = imwarp(uint8(Imgc),tform_rc1,'FillValues',0);
    %figure;imagesc(Img3)
    %crops image s.t. iamges become more symmetric
    [Img4,~] = sng_fishcrop3(Img3,ref3,[size(Img3,2),size(Img3,1)/2]);  %kan meer exact
    [xcom4,ycom4] = sng_CenterOfMassColor(Img4,1);

    %more precise angle calculation
    [theta2,maxangle2,~,ncorr2] = sng_NormalizedRotationCorrelation3(Img4,angles2,[xcom4,ycom4]);
    theta3 = theta1 + theta2 + pi;

    %final rotation
    Trotate = [cos(theta3),-sin(theta3),0;sin(theta3),cos(theta3),0;0,0,1]; 
    tform_rc1 = affine2d(Ttrans*Trotate);
end


%{
figure;imagesc(Img)
figure;imagesc(Imgc)
figure;imagesc(Img2)
figure;imagesc(Img3);sng_imfix
figure;imagesc(uint8(Img4))

figure;plot(ncor)
figure;plot(ncor2)
%}

%% transform such that the center lays on the y axis
        
%transformation with (0.0) as the center coordinate
[I30,ref30] = imwarp(Img,tform_rc1,'FillValues',255);
%reference object with cropped Y-axis such that the center is the real center, x-reference stays the same
ref31 = imref2d([FinalCenteredYCrop,size(I30,2)],ref30.XWorldLimits,[0.5-FinalCenteredYCrop/2,0.5+FinalCenteredYCrop/2]);
%cropped image
%I32 = imwarp(Img,tform_rc1,'FillValues',255,'OutputView',ref31);

tform1ToPixel = affine2d([1,0,0;0,1,0;-ref31.XWorldLimits(1),-ref31.YWorldLimits(1),1]);

tform_rc2 = affine2d(tform_rc1.T*tform1ToPixel.T);

if nargout >=3
    [I34,ref34] = imwarp(Img,tform_rc2,'FillValues',255,'OutputView',imref2d([FinalCenteredYCrop,size(I30,2)]));        
    varargout{1} = I34;
    varargout{2} = ref34;
end

%{
figure;imagesc(ref30.XWorldLimits,ref30.YWorldLimits,uint8(I30));
figure;imagesc(ref32.XWorldLimits,ref32.YWorldLimits,uint8(I32));
figure;imagesc(ref34.XWorldLimits,ref34.YWorldLimits,uint8(I34));
%}




%% output variables

RotationCorrelationOutput(1).name = 'rotation angle1';
RotationCorrelationOutput(2).name = 'raw angle1';
RotationCorrelationOutput(3).name = 'quality 1-peak2/peak1';
RotationCorrelationOutput(4).name = 'correlation1';
RotationCorrelationOutput(5).name = 'angles1';
RotationCorrelationOutput(6).name = 'rotation angle2';
RotationCorrelationOutput(7).name = 'raw angle2';
RotationCorrelationOutput(8).name = 'correlation2';
RotationCorrelationOutput(9).name = 'angles2';

RotationCorrelationOutput(1).value = theta1;
RotationCorrelationOutput(2).value = maxangle1;
RotationCorrelationOutput(3).value = peaksep1;
RotationCorrelationOutput(4).value = ncorr1;
RotationCorrelationOutput(5).value = angles1;
RotationCorrelationOutput(6).value = theta3;
RotationCorrelationOutput(7).value = maxangle2;
RotationCorrelationOutput(8).value = ncorr2;
RotationCorrelationOutput(9).value = angles2;



%% show fasthorizontal correlation with movie

%{
[theta,param,ncorr,I22,I23,I24] = sng_NormalizedRotationCorrelation2(Img2,angles,[xcom2,ycom2],[1:50]);

figure;imagesc(Img);axis equal off tight
set(gca,'position',[0 0 1 1],'units','normalized')
figure;imagesc(Imgc);axis equal off tight
set(gca,'position',[0 0 1 1],'units','normalized')
figure;imagesc(Img2);axis equal off tight
set(gca,'position',[0 0 1 1],'units','normalized')
figure;imagesc(I22);axis off tight equal;
set(gca,'position',[0 0 1 1],'units','normalized')
figure;imagesc(I23);axis off tight equal;
set(gca,'position',[0 0 1 1],'units','normalized')

figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
set(gcf,'color','white')
k=50

writerObj1 = VideoWriter('rotation correlation','MPEG-4');
writerObj1.FrameRate = 10;    % to perform realtime movie
open(writerObj1);

subplot(2,2,1);imagesc(I22);axis off tight equal;
subplot(2,2,2);imagesc(I23);axis off tight equal;
for k=1:50
    subplot(2,2,2);imagesc(I24{k});axis off tight equal
    subplot(2,1,2);plot(angles(1:k),ncorr(1:k));
    xlim([0,360]);ylim([0,1])
    xlabel('angle')
    frame = getframe(gcf);
    writeVideo(writerObj1,frame);    
end
close(writerObj1);    
%}


end