function [Img2,ref2] = sng_fishcrop3(Img,ref,WantedSize)
% this function crops or extends a figure to the input coordinate range


% smaller than WantedSize, the function does nothing.
%
% The xcoord of the center of the square is the middle of the imput figure 
% The ycoord of the center of the square is accoording the spatial reference
% object
%
% This function doenst work for a reference object with a PixelExtentInWorld
% not 1.
%
% example:
% [I76,ref76] = sng_fishcrop3(I75,ref75,[1,1500;1,1000]);

%{
 Img = I75(:,:,1);
 ref = ref75
%}

%{
%I start to change the function but didnt finish it.

cl = [1,1500;1,1000]

Img2 = zeros((cl(2,2) - cl(2,1) + 1),(cl(1,2) - cl(1,1) + 1));

figure;imagesc(ref75.XWorldLimits,ref75.YWorldLimits,uint8(I75));

xwl = round(ref75.XWorldLimits)
ywl = round(ref75.YWorldLimits)

xframe1 = max([xwl(1);cl(1,1)])
xframe2 = min([xwl(2);cl(1,2)])
yframe1 = max([ywl(1);cl(2,1)])
yframe2 = min([ywl(2);cl(2,2)])

Img2(yframe1:yframe2,xframe1:xframe2) = Img(yframe1-ywl+1:yframe2-ywl,xframe1-xwl+1:xframe2-xwl);

figure;imagesc(Img2)
%}







ImageSize = size(Img);
SizeDifference = WantedSize(1:2) - ImageSize(1:2);
x = 1;
y = 1;
height = min(WantedSize(1),ImageSize(1));
width = min(WantedSize(2),ImageSize(2));


%for crop in Y direction dependend on world coordinate center (rows)
if SizeDifference(1) < 0
    %round because we only crop full pixels
    %calculates start postionen of y depending on the world coordinates
    y = -round(ref.YWorldLimits(1)+WantedSize(1)/2);

    %to prevent a larger crop than the wanted size
    if y > round(-SizeDifference(1))
        y = round(-SizeDifference(1));
    end
    %no crop when the wanted image is larger
    if y < 1; 
        y = 1;
    end
end
    
%for crop in X direction which crop equal on both sites 
%for an odd crop, it could differ one pixel
if SizeDifference(2) < 0  
    x = round(-SizeDifference(2)/2);
end

%Img2 = imcrop(Img,[x,y,width - 1,height - 1]);
Img2 = Img(y:y + height - 1,x:x + width - 1,1:ImageSize(3));



%set spatial reference object
ref2 = imref2d(size(Img2));

ref2.XWorldLimits(1) = ref.XWorldLimits(1) + x - 1;
ref2.XWorldLimits(2) = ref2.XWorldLimits(1) + width;

ref2.YWorldLimits(1) = ref.YWorldLimits(1) + y - 1;
ref2.YWorldLimits(2) = ref2.YWorldLimits(1) + height;

end



