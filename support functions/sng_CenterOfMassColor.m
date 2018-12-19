function [xc,yc] = sng_CenterOfMassColor(Img,varargin)
%computes the center off mass
%input is a single image or with coordinates in array

Img = single(Img);
if nargin <= 2
    if varargin{1} == 1
        Img2 = (Img(:,:,1)+Img(:,:,3)+Img(:,:,3))/3;
    else
        Img2 = sqrt(Img(:,:,1).^2+Img(:,:,3).^2+Img(:,:,3).^2);
    end
end

totalintensity = sum(Img2(:));


%input is image

[X,Y] = meshgrid(1:size(Img2,2),1:size(Img2,1));
xc = sum(X(:).* Img2(:))/totalintensity;
yc = sum(Y(:).* Img2(:))/totalintensity; 
        
        
        
        
        
            
end










