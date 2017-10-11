%clear all;close all

%PARAMETERS
zfinput = struct('stage',[],'substage',[],'name',[]','value',[],'sensitivity',[])

zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','levels',1,'low'); %number of scaled levels correlation is performed
zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','iterations',10,'low'); %number of iterations iat is performed
zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','scale',1/16,'low'); %lowered initial scale to increase computation speed for correlation
zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','threshold',0.96,'severe'); %threshold for correlation selection
zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','threshold2',2,'severe'); %threshold for warp modulus


zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','sigmalp',1,'moderate'); %sigma for bandpassfilter (lowpass)
zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','sigmahp',4,'moderate'); %for bandpassfilter (highpass)
zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','scaleC',1/4,'low'); %lowered initial scale to increase computation speed for correlation
zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','levelsC',2,'low'); %number of scaled levels correlation is performed
zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','iterationsC',10,'low'); %number of iterations iat is performed
zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','scaleS',1/4,'low'); %lowered initial scale to increase computation speed for correlation
zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','levelsS',3,'low'); %number of scaled levels correlation is performed
zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','iterationsS',20,'low'); %iterations iat is performed


%% Load images (and spot roi)

%[Names,FolderPath] = uigetfile('/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Data Nynke/','MultiSelect','on')

%FolderPath = uigetdir('/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Data Nynke/')
%FolderPath = uigetdir('C:/Users/260018/Documents/dropbox/');
FolderPath = uigetdir('\\bioinf-filesrv2\cluster15\Ham\Research_group')





[imageinfo] = ImageInfoLink(FolderPath,zfinput);

[stackinfo] = StackInfoLink(imageinfo);


for k1 = 1:numel(batchinfo)
    %% select slices   
    selec = find([imageinfo.fishnumber] == 1);
    ImageSlice = cell(1,numel(selec));
    for k2 = 1:numel(selec)
        ImageSlice{k2} = imread([FolderPath,'/',imageinfo(selec(k2)).name]);
        ImageSlice{k2} = ImageSlice{k2}(:,:,1:3);
    end
        
    %sng_SaveCell2TiffStack(ImageSlice,[FolderPath,'/Stack/',batchinfo(k1).batchname,'.tif'])    
    
    [CorrectedSlice,ppoutput] = PreprocessionLink(ImageSlice,zfinput);
    
    sng_SaveCell2TiffStack(CorrectedSlice,[FolderPath,'/Stack_pp/',batchinfo(k1).batchname,'.tif'])    
    batchinfo(k1).Preprocession = ppoutput;
end

save([FolderPath,'/Stack_pp/imageinfo.mat'],'imageinfo')
save([FolderPath,'/Stack_pp/batchinfo.mat'],'stackinfo')
save([FolderPath,'/Stack_pp/inputparameters.mat'],'zfinput')






%RoiPath = sng_FilesFromMap(FolderPath,'roi');
