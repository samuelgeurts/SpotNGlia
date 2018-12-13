function [varargout] = sng_scale4(scale,varargin)
%resize if input image with scale
%multiply with scale if input is 1d array or matrix with one dimension
%equal to two, which is the case for coordinates

%example:
%scale = 1/8;
%coordinate = [234,555;546,456]
%Img1 = zeros(100)
%Img2 = ones(200)
%[Img1 Img2 coordinate] = sng_scale3(1/8,Img1 Img2 coordinate)

varargout = cell(nargin-1,1);

for j = 1:nargin-1    
    if (size(varargin{j},1) <= 3 ||...
            size(varargin{j},2) <= 3);
        varargout{j} = varargin{j} * scale;
    else
        varargout{j} = imresize(varargin{j},scale);
    end
end
