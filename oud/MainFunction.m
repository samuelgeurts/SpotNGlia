function [parameter] = MainFunction(ImageSlice)

%preprocession
%extended dept of field
%registration



%

%% Main function applied on new data single fishes

%% Load images (and spot roi)
%FolderPath = uigetdir('/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Data Nynke')
FolderPath = uigetdir('/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Data Nynke/Template brain/3 dpf/original')
[~,FolderName] = fileparts(FolderPath)
%Savepath = ['/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Results/',FolderName]
[ImagePath,dirinfo] = sng_FilesFromMap(FolderPath,'tif');


%P.Folderpath = FolderPath
%P.FolderName = FolderName
%P.ImagePath = ImagePath
%P.dirinfo = dirinfo


for k=1:numel(ImagePath)
    %num = sscanf(dirinfo(k).name, '%*s %*s %d'); %werkt nu niet
    %TODO: retrieve number of image to name new images    
    ImageSlice{k} = imread(ImagePath{k}); 
end

RoiPath = sng_FilesFromMap(FolderPath,'roi');

%% template
CompleteTemplate = LoadTemplateLink(3) %input is age
P.UsedTemplate = CompleteTemplate

%}



%% Check Images
%TODO check images for proper resolution, enlighting, color


% RGB and slice correction
[CorrectedSlice,Parameters_preproc] = PreprocessionLink(ImageSlice)


%{
for ki=1:numel(CorrectedSlice)
    figure;imshow(uint8(CorrectedSlice{ki}));
end
figure;imshow(uint8(Images.Combined));
%}        
% Combine Slices, extended dept of field
[Icombined,Parameters_edof] = ExtendedDeptofFieldLink(CorrectedSlice);
%{
figure;imagesc(uint8(Icombined));sng_imfix
%}


%% Allignment

[Ialligned,Parameters_allign] = AllignmentLink2(Icombined,CompleteTemplate);

%{
figure;imagesc(uint8(Ialligned))
%}

%% Brain segmentation
[Ibrain,Parameters_brain] = MidBrainDetectionLink2(Ialligned,Parameters_allign,CompleteTemplate);

%{
figure;imagesc(uint8(Ialligned))
figure;imagesc(uint8(Ibrain))
bc = bwboundaries(Ibrain)
figure;imagesc(uint8(Ialligned))
hold on;plot(bc{1}(:,2),bc{1}(:,1))
%}    

%% Spot Detection
%TODO:  apply on original images which are not dept of field combined as 
%       some spots appear due to dof artifacts. Combine images afterwards
%TODO:  filter spots based on color
%TODO:  optimize many parameters
%       minspotsize,
%TODO:  add more parameters

[Ispots,Parameters_spots] = SpotDetectionLink(Ibrain,Ialligned);

%{
figure;imshow(Ispots)
%}


%
figure;imagesc(uint8(Ialligned));sng_imfix
hold on;
for ka = 1:Parameters_spots.KThreshold
    spts(ka,:) = getfield(Parameters_spots.SpotProperties,{ka},'Centroid',{1:2})
end
sq = scatter(spts(:,1),spts(:,2),300)
sq.LineWidth = 0.2
sq.MarkerEdgeColor = 'white'
sq.Marker = 'o'
bc = bwboundaries(Ibrain)
plot(bc{1}(:,2),bc{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
%}

%

%% compare with Nynke spot-roi's
%function SpotCompareLink(Ispots,Icombined)

[sROI] = ReadImageJROI(RoiPath{1})
sc = sROI.mnCoordinates %spot coordinates
tf = Parameters_allign.tform_complete %complete fish transformation
[xsc,ysc] = transformPointsForward(tf,sc(:,1),sc(:,2));
bc = bwboundaries(Ibrain)

%
figure;imshow(uint8(Icombined))
hold on
scatter(sc(:,1),sc(:,2))

figure;imagesc(uint8(Ialligned))
hold on
scatter(xsc,ysc)
%}


figure;imagesc(uint8(Ialligned - Ispots))
hold on
sq = scatter(xsc,ysc,50)
sq.LineWidth = 2
sq.MarkerEdgeColor = 'white'
sq.Marker = 'o'
plot(bc{1}(:,2),bc{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
axis equal tight off
set(gca,'position',[0 0 1 1],'units','normalized')
truesize

figure;imagesc(uint8(Ispots));axis off tight equal
hold on;
plot(bc{1}(:,2),bc{1}(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
sq = scatter(xsc,ysc,50)
sq.LineWidth = 2
sq.MarkerEdgeColor = 'white'
sq.Marker = '+'
set(sq,'MarkerEdgeColor','white','LineWidth',1.5)
axis equal tight off
set(gca,'position',[0 0 1 1],'units','normalized')
truesize



%}
    

end
    