function Img2 = sng_RGBselection(Img,sigma1,sigma2)


if ~exist('sigma1')
    sigma1 = 1
end

if ~exist('sigma2')
    sigma2 = 4
end


%bandpassfilter RGB image to focus on spots
%crops image to remove boundary and scale bar
%calculates center of mass according bandpass filtered image
%
%{
Img = ImageSlice{k}(:,:,1:3);
%}
%Img = imread(Imgloc1{1}{1});
%Img = Img(:,:,1:3);

%gaussian blur and difference image
Img = double(Img(:,:,1:3));
%Imglow1 = imgaussfilt(Img,sigma1);
%Imglow2 = imgaussfilt(Img,sigma2);
%for older matlab versions:
    h1 = fspecial('gaussian',2*ceil(2*sigma1)+1,sigma1);
    Imglow1 = imfilter(Img,h1,'replicate');
    h2 = fspecial('gaussian',2*ceil(2*sigma2)+1,sigma2);
    Imglow2 = imfilter(Img,h2,'replicate');
Imgband =  (Imglow1 - Imglow2);

%geen balkje en geen randen
Img1 = Imgband(20:round(size(Img,1)*(11/12)),20:size(Img,2)-20,1:3);
Img2a = Img1;
Img2b = Img1;

Img2a(Img1 <= 0) = 0;
Img2b(Img1 >= 0) = 0;

Img2 = cat(1,Img2a,-Img2b);


%figure;imagesc(uint8(Img2*10));


% %select area around center of mass detail image
% %does not work for some fishes
% [comx,comy] = sng_CenterOfMassColor(Img1,1);
% x1 = round(comx-300);if x1 < 1; x1 = 1;end;
% x2 = round(comx+300);if x2 > size(Img1,2); x2 = size(Img1,2);end;
% y1 = round(comy-300);if y1 < 1; y1 = 1;end;
% y2 = round(comy+300);if y2 > size(Img1,1); y2 = size(Img1,1);end;
% Img2 = Img1(y1:y2,x1:x2,:);

end

%{
figure;imagesc(uint8(Img(:,:,1:3)));


figure;imagesc(0.15*abs(Img1(:,:,1:3)))
figure;imagesc(0.15*Img2a(:,:,1:3))
figure;imagesc(-0.15*Img2b(:,:,1:3))

l = Img1-min(Img1(:));

figure;imagesc(QImg);
figure;imagesc(Img2(:,:,1));colormap('gray')
figure;imagesc(Img2(:,:,2));colormap('gray')
figure;imagesc(Img2(:,:,3));colormap('gray')
sng_figureslide
%}

% [Img2a,ECCWarp] = sng_RGB(Img1b,1/2,'translation',2,20);
% %[Img2b,tform] = sng_RGB_imreg(Img1b,1,'translation');
% figure;imagesc(10*Img2a)
% figure;imagesc(10*Img2b)
% 
% QImg = Img(:,:,1:3);
% for j=2:3
% QImg(:,:,j) = imwarp(QImg(:,:,j),tform{j},'OutputView',imref2d(size(QImg)));
% end
% 

