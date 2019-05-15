function croppedImage = sng_CropExtendAroundCoord(inputImage, coords, boxsize, method, color)

%extends images triple the size and crops to size 'boxsize' around a given center 'coords'

%20181105 Version sng_cropAroundCoord comes from the function
%sng_boxaroundcenter2. As it is a very specific function, it is
%changed to a static function under the SNGMidBrinDetection class
%20190507 Version sng_CropExtendAroundCoord

%{
inputImage = Image(:,:,5:7);
coords = obj.BrainSegmentationObject(1).midbrainCenter)
%}


imageSize = size(inputImage);
extendedCoords = coords + fliplr(imageSize(1:2));

imageExtention = ((boxsize - 1) / 2);
croppedCoords = floor((imageExtention));

columns = extendedCoords(1) - floor(imageExtention(1)):extendedCoords(1) + ceil(imageExtention(1));
rows = extendedCoords(2) - floor(imageExtention(2)):extendedCoords(2) + ceil(imageExtention(2));





if ndims(inputImage) == 2
    imageSize(3) = 1;
end

if ~exist('coords', 'var') || isempty(coords)
    coords = floor(imageSize(1:2)/2);
end

if ~exist('boxsize', 'var') || isempty(boxsize)
    boxsize = [1000, 1000]; 
end

if ~exist('method', 'var') || isempty(method)
    method = 'mode';
end
%inputImage, coords1, boxsize, exteriorColor


switch method
    case 'mirror'
        extendedImage = mirrorExtend(inputImage);
    case 'mode'
        extendedImage = modeExtend(inputImage);
    case 'single'
            
            %if only a single color value is given and the image is rgb
            %'method' is extended to a 1:3 array 
            if length(color) < imageSize(3)
                color(1:imageSize(3)) = color(1);
            end
            
        extendedImage = singleExtend(inputImage);
     otherwise
        error('extention method unknown')
end

croppedImage = crop(extendedImage);




    function extendedImage = mirrorExtend(inputImage)
        
        %mirror extention
        lr = fliplr(inputImage);
        ud = flipud(inputImage);
        lrud = fliplr(ud);
        extendedImage = [lrud, ud, lrud; lr, inputImage, lr; lrud, ud, lrud];
        
    end
    function extendedImage = modeExtend(inputImage)
        
        %select border pixels for all color channels
        mask = true(imageSize(1:2));
        mask(2:end-1,2:end-1) = false;
        
        for ichannel = 1:imageSize(3)
            channel = inputImage(:,:,ichannel);
            borderPixels(:,ichannel) = channel(mask);
        end

        mostFrequentValue = mode(borderPixels);
        
        extendedImage = zeros([3*imageSize(1:2),imageSize(3)],class(mostFrequentValue)); %preallocation
             
        for ichannel = 1:imageSize(3)
            extendedImage(:,:,ichannel) = (mostFrequentValue(ichannel));
        end
        
        extendedImage(imageSize(1)+1:imageSize(1)+imageSize(1), imageSize(2)+1:imageSize(2)+imageSize(2),:) = inputImage;
    end
    function extendedImage = singleExtend(inputImage)
        
        extendedImage = zeros([3*imageSize(1:2),imageSize(3)]);
        for ichannel = 1:imageSize(3)
            extendedImage(:,:,ichannel) = color(ichannel) * ones(3*imageSize(1:2));
        end
        
        
        extendedImage(imageSize(1)+1:imageSize(1)+imageSize(1), imageSize(2)+1:imageSize(2)+imageSize(2),:) = inputImage;

        
    end
    function croppedImage = crop(extendedImage)
        croppedImage = extendedImage(rows, columns, :);
    end

    function show
       figure(1);
       imagesc(inputImage);
       hold on; 
       scatter(coords(1),coords(2),'Marker','x','MarkerEdgeColor','red');
       hold off
       axis equal off

       figure(2);
       imagesc(extendedImage);
       hold on; 
       scatter(extendedCoords(1),extendedCoords(2),'Marker','x','MarkerEdgeColor','red');
       hold off
       axis equal off

       figure(3);
       imagesc(croppedImage);
       hold on; 
       scatter(croppedCoords(1),croppedCoords(2),'Marker','x','MarkerEdgeColor','red');
       hold off
       axis equal off
    end

end

