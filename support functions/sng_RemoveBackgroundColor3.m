function [Img2,BackgroundMask,thr] = sng_RemoveBackgroundColor3(Img,method,varargin)
%removes background using


% 
% example:
% Img = I10;
% I20 = sng_RemoveBackgroundColor3(I10,'TriangleSmooth',1,'cuboid');
% 
% 
% I20 = sng_RemoveBackgroundColor3(I10,SmoothingMethod,Variable,ColorCombiningMethod);
% 
% SmoothingMethod:    'TriangleSmooth'
%   Variable            [naturalnumber]
%                       amount of smoothing of the histogram
%                    
% Method:             'cuboid'
%                     In rgb space, the thersholded color area is cuboid
%                     'ellipsoid'
%                     In rgb space, the threshollded color area is ellipsoid 
%                     with the center at [255,255,255]                     
% 




%{

method = 'TriangleSmooth';
Img = I10;
var1 = 1
thr_method = 'cuboid'

%}



var1=[];
Img2 = zeros(size(Img));
thr = zeros(3,1);

%variable used in grayscale RemoveBackground
if nargin>=3
    var1 = varargin{1};
end

%variable to determine cubic of eliptical histogram threshold 
if nargin>=4
    thr_method = varargin{2};
end

%three times graythreshold
if strcmp(thr_method,'separate') || strcmp(thr_method,'cuboid')
    for j = 1:3    
        [Img2(:,:,j),thr(j)] = sng_RemoveBackground(Img(:,:,j),method,var1);    
    end
end

if strcmp(thr_method,'cuboid')
     A = Img2(:,:,1)==255; 
     B = Img2(:,:,2)==255; 
     C = Img2(:,:,3)==255;
     BackgroundMask = A|B|C;
end

if strcmp(thr_method,'ellipsoid')
    for j = 1:3
        [~,thr(j)] = sng_RemoveBackground(Img(:,:,j),method,var1);    
    end
    dist = ((single(Img(:,:,1)).^2/(thr(1))^2)+...
        (single(Img(:,:,2)).^2/(thr(2))^2)+...
        (single(Img(:,:,3)).^2/(thr(3))^2));
    BackgroundMask=dist>=3;
end

Img(cat(3,BackgroundMask,BackgroundMask,BackgroundMask))=255;
Img2=Img;

%      Img(1,1,1:3) = [255 255 255];
%      Img(2,1,1:3) = [180 180 180];    
%      Img(3,1,1:3) = [ 0 0 180];   
%      Img(4,1,1:3) = [ 182 180 178];   

     
end