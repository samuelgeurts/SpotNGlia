function [Ibrain,Parameters] = MidBrainDetectionLink3(Ialligned,CompleteTemplate,zfinput);
%Detect midbrain
%works in combination with other "Link" functions
%based on MidBraindetection 2 and GeneratePolarTransforms2
%
%Version MidBrainDetectionLink2
%   applied on template aligned fishes
%   add skeletonband for path selection
%   put 1000x1000 selection in other function (sng_boxaroundcenter)
%
%Version MidBrainDetectionLink3
%   applicable for zfinput and broutput


if ~exist('zfinput','var');
    zfinput = struct;
    zfinput = sng_zfinput(zfinput,0,'BrainSegmentation','','Method','TriangleSmooth','?');
    zfinput = sng_zfinput(zfinput,0,'BrainSegmentation','','Method2',4,'?');2
end
sng_zfinputAssign(zfinput,'BrainSegmentation')



%% generate useable edge filter from edge template
%add more templates later on
filterwidth = 21; %has to be an odd numberl

EdgeTemp = CompleteTemplate.EdgeTemplate1;
EdgeTempMean = mean(EdgeTemp,2);
EdgeFilter = repmat(EdgeTempMean,1,filterwidth);


%{
figure;
imagesc(uint8(EdgeTemp));axis off tight equal
set(gca,'position',[0 0 1 1],'units','normalized')
truesize(gcf,12*sng_size(get(get(gca,'Children'),'Cdata'),[1,2]))

figure;imagesc(uint8(EdgeTemp));axis off equal tight
figure;imagesc(uint8(EdgeTempMean));sng_imfix

figure;
imagesc(uint8(EdgeFilter));axis off tight equal
set(gca,'position',[0 0 1 1],'units','normalized')
truesize(gcf,12*sng_size(get(get(gca,'Children'),'Cdata'),[1,2]))

%}

%% polar transform

cxy = CompleteTemplate.CenterMidBrain;

[Isquare] = sng_boxaroundcenter(Ialligned,cxy);

    %extend with 200 pixels on every side
    siz = size(Ialligned);
    cxy2 = round(cxy + 200); 
    rangex = cxy2(1)-499:cxy2(1)+500;
    rangey = cxy2(2)-499:cxy2(2)+500;

%{
sng_show(Isquare)
%}

%polar coordinates
Ipolar = sng_Im2Polar3(Isquare);
Ipolar = permute(Ipolar,[2,1,3]);

%{
%view(gca,[180 90]);
sng_show(Ipolar)
%}

%TODO: normalized or not normalized correlation?

%correlation with Edgetemplates
ICorr(:,:,1) = normxcorr2_general(EdgeFilter(:,:,1),Ipolar(:,:,1),numel(EdgeFilter(:,:,1))/2);
ICorr(:,:,2) = normxcorr2_general(EdgeFilter(:,:,2),Ipolar(:,:,2),numel(EdgeFilter(:,:,2))/2);
ICorr(:,:,3) = normxcorr2_general(EdgeFilter(:,:,3),Ipolar(:,:,3),numel(EdgeFilter(:,:,3))/2);

%{
figure;imagesc((ICorr+1)/2);sng_imfix
%}

%crop the correlated image to the original size
stf = floor(size(EdgeFilter(:,:,1))/2); %start coordinates to add
stc = ceil(size(EdgeFilter(:,:,1))/2); %end coordinates to add
si = size(Ipolar); %image size
crR = stc(1):si(1)+stf(1); %row coordinates
crC = stc(2):si(2)+stf(2); %column coordinates
ICorr2 = ICorr(stc(1):si(1)+stf(1),stc(2):si(2)+stf(2),:);

%remove values below zero or lift above zero
%ICorr2 (ICorr2 < 0) = 0;
ICorr2 =  (ICorr2+1)/2; %function lift


%{
figure;imagesc(ICorr2);sng_imfix
%}

%combine channels
%ISum = sum(ICorr2,3);
INorm = sqrt(sum(ICorr2.^2,3));

%{
figure;imshow(ISum)
figure;imagesc(INorm);colormap gray;sng_imfix
figure;imagesc(INorm);sng_imfix

%}

%DONE: use more suffisticated method to exclude area in shortest path
%finding, multiply shortest path with based on brains blurred distancemap
%to exclude some area from pathfinding
%TODO: add band to generateskeleton function output instead of computing from distancemap
IsquareMB = sng_boxaroundcenter(CompleteTemplate.MidBrainDistanceMap,cxy);

%{
figure;imagesc(boolean(IsquareMB));sng_imfix;colormap gray
%}

IDistmap = sng_Im2Polar3(IsquareMB);
IDistmap = permute(IDistmap,[2,1,3]);
IDistmap(IDistmap ~= 0) = 1; %set every nonzero value 2 one (band) 

%gaussion blur on IDistmap, probably doesnt work as the main plan is
%blurred
IDistmap2 = imgaussfilt(IDistmap,[40,5]);

%IDistmap3 contains added gausian filter at te boudary i.e. the values of
%IDistmap which are 1 stays 1, i.e. the plane is preserved
blurvar = 6;
IDistmap3 = IDistmap;
sz=size(IDistmap3)
center = 501;
for k2 = 1:sz(2)
a = find(IDistmap3(:,k2) == 1,1,'first');
b = find(IDistmap3(:,k2) == 1,1,'last');
%y = gaussmf([-500:1:500],[(b-a)/blurvar 0])';
%for versions without fuzzy logic toolbox
y = exp(-((-500:1:500)-0).^2/(2*((b-a)/blurvar).^2))';

IDistmap3(1:a,k2) = y(center-a+1:center);
IDistmap3(b:sz(1),k2) = y(center:center+sz(1)-b);
end
%a bit of gaussian blur added
IDistmap4 = imgaussfilt(IDistmap3,[5,5]);
%figure;plot(y) % a gaussian

INorm2 = INorm .* IDistmap4;

%{
figure;imagesc(IDistmap);sng_imfix;colormap gray
figure;imagesc(IDistmap2);sng_imfix;colormap gray
figure;imagesc(IDistmap3);sng_imfix;colormap gray
figure;imagesc(IDistmap4);sng_imfix;colormap gray

figure;imagesc(INorm2);sng_imfix;
%}

%invert for lowest cost path and normalize
INorm3 = -INorm2/max(INorm2(:)) + 1;

%{
figure;imagesc(INorm3);sng_imfix;
%}


%% Find shortest Path

s = size(INorm3)
%compute edges
[edges] = sng_SquareGridDiGraph(si(1:2));
%compute directer graph
DGraph = digraph([edges(:,1)],[edges(:,2)],INorm3(edges(:,2)));


%compute shortest path start and end at equal row
rows = 100:10:400
nrows = numel(rows)
path = zeros(nrows,1000);


%TODO: check if area calculation is right
%figure;imagesc(INorm3);sng_imfix

for l = 1:nrows
    l
    a = sub2ind(s,rows(l),1);
    b = sub2ind(s,rows(l),s(2));   
    [path,lengthd(l)]= shortestpath(DGraph,a,b,'Method','acyclic');
    [I,J] = ind2sub(s,path);
    
    Area(l) = sum(pi*(I.^2)/1000);
    Hor(l) = I(end) + I(end/2);
    Ver(l) = I(round((3*s(2))/4)) + I(round((s(2))/4));

     %hold on
     %plot(J,I,'r','LineWidth',2)
     %drawnow
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
[mx,mr] = min(lengthd)
a = sub2ind(s,rows(mr),1);
b = sub2ind(s,rows(mr),s(2)); 
[path,~]= shortestpath(DGraph,a,b,'Method','acyclic');
[I,J] = ind2sub(s,path);

%{
figure;imagesc(INorm3);sng_imfix;hold on;plot(J,I,'r','LineWidth',2)
figure;imagesc(INorm2);sng_imfix;hold on;plot(J,I,'r','LineWidth',2)
%}



%% create mask in cartesian space
path2 = zeros(s);
path2(path) = 1;       
%make everything below the path one 
for kk = 1:size(path2,2);
    path2(1:I(kk),kk) = 1;
end

%{
    figure;imagesc(path2);sng_imfix
%}


path3 = permute(path2,[2,1]); %restore original permutation
Mask =  sng_Polar2Im(path3);
    
%{
    figure;imagesc(path3)
    figure;imagesc(uint8(Mask));sng_imfix

%}
    
frame = false(siz(1:2)+400);
frame(rangey,rangex) = Mask;
Ibrain = frame(201:end-200,201:end-200);

%{
    figure;imagesc(uint8(Ialligned))

    figure;imagesc(uint8(Ibrain))
    figure;plot(J,I)

    IA = double(Ialligned);
    IA(:,:,3) = IA(:,:,3) + 110 * double(IBrain);
    IA(:,:,2) = IA(:,:,2) + 50 * double(IBrain);
    figure;imagesc(uint8(IA));sng_imfix

bc = bwboundaries(Ibrain)
figure;imagesc(uint8(Ialligned));sng_imfix
hold on;plot(bc{1}(:,2),bc{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',2)

%}


bwbb = bwboundaries(Ibrain);

%%

Parameters.EdgeFilterWidth = filterwidth;
%Parameter.BrainArea
%Parameter.BrainHorizontalsection
%Parameter.BrainVerticalsection
Parameters.ShortestPath = [J',I'];
Parameters.ShortestPathValue = mx;

Parameters.BrainEdge = bwbb{1}








%% the fish movie (not finished, look for mov_polartransform3)
%{
siz = size(Ialligned)
cxy2 = round(cxy + 200) 
rangex = cxy2(1)-499:cxy2(1)+500;
rangey = cxy2(2)-499:cxy2(2)+500;


figure;imagesc(uint8(Isquare));axis equal tight off


tf = 200 %totalframes
scale = 1
viewrot = 90
cutrot = 90

Img = Isquare;
Img = imrotate(Img,cutrot);

h= figure(1);

R=Img(:,:,1);
G=Img(:,:,2);
B=Img(:,:,3);

s = size(Img);
[X,Y] = meshgrid(1:s(2),1:s(1));

%make center zero
X2 = X - cxy2(2);
Y2 = Y - cxy2(1);

%image to polar
Radius = sqrt(X2.^2+Y2.^2);
Theta = atan2(Y2,X2); %branchcut on the left -pi -> pi
r = scale * 500;%readius of circle

%make area outsite radius white
R2=double(R);R2(Radius >= r) = 256;
G2=double(G);G2(Radius >= r) = 256;
B2=double(B);B2(Radius >= r) = 256;

writerObj1 = VideoWriter('polar tranform2','MPEG-4');
writerObj1.FrameRate = 30;    % to perform realtime movie
open(writerObj1);

%the transformation parameters. s is for the angle transformation. 
s = linspace(1,0.01,tf).^2; %final s should be small but nog 0
% s2 isfor the radius transformation based on the values of s
s2 = (linspace(0,10,tf).^2)./s
%used for the x-limit
sxmin = linspace(1,0,tf).^2; 


HX = Radius .* cos(Theta);
HY = Radius .* sin(Theta);
figure(1);
ScatH = scatter(HX(:),HY(:),100,double(cat(2,R2(:),G2(:),B2(:)))/255,'.','LineWidth',2);
view(gca,[viewrot 90]);
axis tight off equal

for k = 1:numel(s)

HX = (Radius+s2(k)) .* cos(Theta*s(k));
HY = (Radius+s2(k)) .* sin(Theta*s(k));

set(ScatH,'Xdata',HX(:))
set(ScatH,'Ydata',HY(:))

%xlimit such that the initial center of the image is smootly shifted to the
%right. if k = 1, sxmin(1)=1 and xmin = 2*HX(end/2,end/2) - max(HX(R2 ~= 255)));
% if k = 100, sxmin(1)=0 and xmin = HX(end/2,end/2), quadrating makes it smoother
xmin = HX(floor(end/2),floor(end/2)) + sxmin(k)*(HX(floor(end/2),floor(end/2))- max(HX(R2 ~= 256)));
xlim([xmin,max(HX(R2 ~= 256))]);
ylim([min(HY(R2 ~= 256)),max(HY(R2 ~= 256))]);

axis square off


%xlim([xmin,HX(floor(end/2),end)]);
%min(HX(R2~=256))


%pause (0,1)

frame = getframe(gcf);
writeVideo(writerObj1,frame); 
end

close(writerObj1);    

%% from fish to circel in fish
k=1
HX = (Radius+s2(k)) .* cos(Theta*s(k));
HY = (Radius+s2(k)) .* sin(Theta*s(k));

%inital fish
figure;scatter(HX(:),HY(:),100,double(cat(2,R(:),G(:),B(:)))/255,'.','LineWidth',2);
axis equal tight off;
view(gca,[viewrot 90]);

%with circle
figure;scatter(HX(:),HY(:),100,double(cat(2,R(:),G(:),B(:)))/255,'.','LineWidth',2);
axis equal tight off
pos = [-500,-500, 1000, 1000]*scale;
rectangle('Position',pos,'Curvature',[1 1]);
view(gca,[viewrot 90]);

%remove out circle
figure;scatter(HX(:),HY(:),100,double(cat(2,R2(:),G2(:),B2(:)))/255,'.','LineWidth',2);
axis equal tight off
pos = [-500,-500, 1000, 1000]*scale;
rectangle('Position',pos,'Curvature',[1 1]);
view(gca,[viewrot 90]);

%scale movie
figure;scatter(HX(:),HY(:),100,double(cat(2,R2(:),G2(:),B2(:)))/255,'.','LineWidth',2);
axis equal tight off
pos = [-500,-500, 1000, 1000]*scale;
rectangle('Position',pos,'Curvature',[1 1]);
view(gca,[viewrot 90]);



writerObj2 = VideoWriter('polar tranform zoom','MPEG-4');
writerObj2.FrameRate = 30;    % to perform realtime movie
open(writerObj2);

z = linspace(750*scale,500*scale,50);
xlim([-z(1),z(1)])
for j = 1:numel(z)
xlim([-z(j),z(j)])
ylim([-500*scale,500*scale])


frame = getframe(gcf);
writeVideo(writerObj2,frame); 
end

close(writerObj2);   

%}

end
