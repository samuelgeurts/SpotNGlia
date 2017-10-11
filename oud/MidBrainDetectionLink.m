function [IBrain,Parameters] = MidBrainDetectionLink(Ialligned,Parameters_allign,CompleteTemplate);
%Detect midbrain
%works in combination with other "Link" functions
%based on MidBraindetection 2 and GeneratePolarTransforms2


%% generate useable edge filter from edge template
%add more templates later on
filterwidth = 5

EdgeTemp = CompleteTemplate.EdgeTemplate1
EdgeTempMean = mean(EdgeTemp,2)
EdgeFilter = repmat(EdgeTempMean,1,filterwidth)

%{
figure;imshow(uint8(EdgeTemp))
figure;imshow(uint8(EdgeTempMean))
figure;imshow(uint8(EdgeFilter))
%}

cxy = Parameters_allign.TemplateCenterMidBrain

%extend with 200 pixels on every side
siz = size(Ialligned)
ext = 255*ones(siz(1)+400,siz(2)+400,3);
ext(201:siz(1)+200,201:siz(2)+200,:) = Ialligned;

%{
figure;imagesc(uint8(ext))
hold on
scatter(cxy2(1),cxy2(2))
%}


%square 1000x1000 images with skeleton based center in the middle
cxy2 = round(cxy + 200) 
rangex = cxy2(1)-499:cxy2(1)+500;
rangey = cxy2(2)-499:cxy2(2)+500;
Isquare = ext(rangey,rangex,:);

%{
figure;imagesc(uint8(Isquare));axis equal tight off
%}

%polar coordinates
Ipolar = sng_Im2Polar3(Isquare);
Ipolar = permute(Ipolar,[2,1,3]);

%{
view(gca,[180 90]);
figure;imshow(uint8(Ipolar));
%}

%TODO: normalized or not normalized correlation?

%correlation with Edgetemplates
ICorr(:,:,1) = normxcorr2_general(EdgeFilter(:,:,1),Ipolar(:,:,1),numel(EdgeFilter(:,:,1))/2);
ICorr(:,:,2) = normxcorr2_general(EdgeFilter(:,:,2),Ipolar(:,:,2),numel(EdgeFilter(:,:,2))/2);
ICorr(:,:,3) = normxcorr2_general(EdgeFilter(:,:,3),Ipolar(:,:,3),numel(EdgeFilter(:,:,3))/2);

%crop the correlated image to the original size
stf = floor(size(EdgeFilter(:,:,1))/2); %start coordinates to add
stc = ceil(size(EdgeFilter(:,:,1))/2); %end coordinates to add
si = size(Ipolar); %image size
crR = stc(1):si(1)+stf(1); %row coordinates
crC = stc(2):si(2)+stf(2); %column coordinates
ICorr = ICorr(stc(1):si(1)+stf(1),stc(2):si(2)+stf(2),:);

%{
figure;imshow(ICorr)
%}

%remove values belos zero
ICorr2 = ICorr;
ICorr2 (ICorr2 < 0) = 0;


%{
figure;imshow(ICorr2)
%}

%ISum = sum(ICorr2,3);
INorm = sqrt(sum(ICorr2.^2,3));

%{
figure;imshow(ISum)
figure;imshow(INorm)
%}

%invert for lowest cost path
INorm = -(INorm-1);

%{
figure;imshow(INorm)
%}

%TODO: use more suffisticated method to exclude area in shortest path
%       finding, umultiply shortest path with based on brains blurred 
%       distancemap
%exclude some area from pathfinding

INorm(1:100,:) = 1;
INorm(400:end,:) = 1;

%{
figure;imshow(INorm)
%}


%% Find shortest Path

s = size(INorm)
%compute edges
[edges] = sng_SquareGridDiGraph(s);
%compute directer graph
DGraph = digraph([edges(:,1)],[edges(:,2)],INorm(edges(:,2)));


%compute shortest path start and end at equal row
rows = 100:10:400
nrows = numel(rows)
path = zeros(nrows,1000);


%TODO: check if area calculation is right

for l = 1:nrows
    l
    a = sub2ind(s,rows(l),1);
    b = sub2ind(s,rows(l),s(2));   
    [path,d(l)]= shortestpath(DGraph,a,b,'Method','acyclic');
    [I,J] = ind2sub(s,path);
    
    Area(l) = sum(pi*(I.^2)/1000);
    Hor(l) = I(end) + I(end/2);
    Ver(l) = I(round((3*s(2))/4)) + I(round((s(2))/4));

%     hold on
%     plot(J,I,'r')
%     drawnow
end

%{
    figure;plot(Area)
    figure;plot(Hor)
    figure;plot(Ver)
    figure;plot(d)
%}

%% find lowestcost path from all lowest cost past with equal end and
%start row
%TODO:  apply area/horizontal/vertical distance brain constrained
%       find max value of hor and ver
[mx,mr] = min(d)
a = sub2ind(s,rows(mr),1);
b = sub2ind(s,rows(mr),s(2)); 
[path,d(rows(mr))]= shortestpath(DGraph,a,b,'Method','acyclic');
[I,J] = ind2sub(s,path);


%% create mask in cartesian space
path2 = zeros(s);
path2(path) = 1;       
%make everything below the path one 
for kk = 1:size(path2,2);
    path2(1:I(kk),kk) = 1;
end

path2 = permute(path2,[2,1]); %restore original permutation
Mask =  sng_Polar2Im(path2);
    
%{
    figure;imagesc(path2)
    figure;imagesc(uint8(Mask))
%}
    
frame = false(siz(1:2)+400);
frame(rangey,rangex) = Mask;
IBrain = frame(201:end-200,201:end-200);

%{
    figure;imagesc(uint8(Ialligned))
    figure;imagesc(uint8(IBrain))
    figure;plot(J,I)
    
    IA = double(Ialligned);
    IA(:,:,3) = IA(:,:,3) + 200 * double(IBrain);
    figure;imagesc(uint8(IA))



%}

Parameters.EdgeFilterWidth = filterwidth
%Parameter.BrainArea
%Parameter.BrainHorizontalsection
%Parameter.BrainVerticalsection
Parameters.ShortestPath = [J',I'];
Parameters.ShortestPathValue = mx












