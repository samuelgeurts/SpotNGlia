function [Img2] = sng_boxaroundcenter2(Img1,coords1,boxsize,exteriorColor)
%extends images and cut out a 1000x1000 image around a given center
    %extend with 200 pixels on every side

%{    
    Img1 = CompleteTemplate.MidBrainDistanceMap;
    coords1 = cxy
    
    Img1 = objt.BandMidBrain;
    coords1 = objt.CenterMidBrain;
    
%}    
    siz = size(Img1);
    if ndims(Img1) == 2
        siz(3) = 1;
    end

    if ~exist('exteriorColor')
        exteriorColor = 255;
    end
        
    ext = exteriorColor*ones(siz(1)+400,siz(2)+400,siz(3));
    ext(201:siz(1)+200,201:siz(2)+200,:) = Img1;

    
    %{
    figure;imagesc(uint8(ext))
    hold on
    scatter(cxy2(1),cxy2(2))
    %}


    %square 1000x1000 images with skeleton based center in the middle
    cxy2 = round(coords1 + 200);
    rangex = cxy2(1)-499:cxy2(1)+500;
    rangey = cxy2(2)-499:cxy2(2)+500;
    Img2 = ext(rangey,rangex,:);
    
    %{
    figure;imagesc(uint8(Img2))
    %}
    
   
    %{
    Img3 = ones(1000,1000)*255;
    c1 = round(coords1)
    r1 = max([c1 - 499;1,1])
    r2 = min([c1 + 500;siz(1),siz(2)])
    
    c1-499   substract the first part
    
    Img3(r1(1):r2(1),r1(2):r2(2)) = Img1(r1(1):r2(1),r1(2):r2(2))
    
    figure;imagesc(uint8(Img3))
    %}
end

