function I60 = sng_imrotate(I50,angle)
%instead of imrotate, we use the more general transformation function
%imwarp where sng_imrotate is based on because imrotate shows some 
%unknown artifacts for rotating 90 degrees
 

tform = affine2d([cos(angle) -sin(angle) 0; sin(angle) cos(angle) 0; 0 0 1]);
    RA = imref2d(size(I50));
    Rout = images.spatialref.internal.applyGeometricTransformToSpatialRef(RA,tform);
        %crop to original size
        Rout.ImageSize = RA.ImageSize;
        xTrans = mean(Rout.XWorldLimits) - mean(RA.XWorldLimits);
        yTrans = mean(Rout.YWorldLimits) - mean(RA.YWorldLimits);
        Rout.XWorldLimits = RA.XWorldLimits+xTrans;
        Rout.YWorldLimits = RA.YWorldLimits+yTrans;    
    I60 = imwarp(I50,tform,'bilinear','OutputView',Rout, 'SmoothEdges',true,'FillValues',0);

end