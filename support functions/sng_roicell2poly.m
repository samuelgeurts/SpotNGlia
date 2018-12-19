
function [poly,slice] = sng_roicell2poly(RoiBrain,n)

%selects polygon coordinates (and slices) from imported RoiCell
%check for largest area and place them in order

%Example
%   %put all largest area values in cell structure
%   [poly,slice] = sng_roi2poly(RoiBrain)
%   %give largest area values
%   [poly,slice] = sng_roi2poly(RoiBrain,1)
%   %give n largest area values in cell structure
%   [poly,slice] = sng_roicell2poly(RoiBrain,2)


%preallocate
r = cell(1,numel(RoiBrain));
a = zeros(1,numel(RoiBrain));



for k = 1:numel(RoiBrain)
    r{k} = RoiBrain{k}.mnCoordinates;
    a(k) = polyarea(r{k}(:,1),r{k}(:,2));
end
[~,index] = sort(a,'descend');

if exist('n','var') && (n == 1)
    poly = RoiBrain{index(1)}.mnCoordinates;
    if nargout >= 2
        slice = RoiBrain{index(1)}.vnSlices;
    end
elseif exist(n,'var') && (n >= 2)
    poly = cell(1,n);
    for k2 = 1:n
        poly{k2} = RoiBrain{index(k2)}.mnCoordinates;
        if nargout >= 2
            slice{k2} = RoiBrain{index(k2)}.vnSlices;
        end
    end
else
    poly = cell(1,numel(RoiBrain));
    for k2 = 1:numel(RoiBrain)
        poly{k2} = RoiBrain{index(k2)}.mnCoordinates;
        if nargout >= 2
            slice{k2} = RoiBrain{index(k2)}.vnSlices;
        end
    end
end
    
        
    


    %{
    figure;imshow(uint8(Icombined))
    hold on;plot(mbr(:,1),mbr(:,2))
    %}    
    
    %transform coordinates to template 
    %[mbr2{k1}(:,1),mbr2{k1}(:,2)] = transformPointsForward(tform_1234,mbr(:,1),mbr(:,2));
    %mbr2{k1} = double(mbr2{k1});