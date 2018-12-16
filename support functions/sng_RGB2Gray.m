function Img2 = sng_RGB2Gray(Img,vector,ynimg)
%transforms a color image to a grayscale image by multiplying with a vector
%use the vector with the highest separability to improve contrast

%{
    k10=10
    CorrectedSlice = sng_openimstack2([PreprocessionPath,'/',stackinfo(k10).stackname,'.tif']);           
    Icombined = sng_SliceCombine(CorrectedSlice,stackinfo(k10).ExtendedDeptOfField.IndexMatrix);
    Img = Icombined;
    vector = [-0.431770623113389;0.847397560890843;-0.309016994374947] 
%}    

r = double(Img(:,:,1));
g = double(Img(:,:,2));
b = double(Img(:,:,3));

BWarray = [r(:),g(:),b(:)] * vector;
Img2 = reshape(BWarray,size(Img(:,:,1)));

if (nargin >= 3) && (ynimg == true)
    figure;imagesc(Img);
    figure;imagesc(Img2);colormap gray
end
