clear all;close all

%PARAMETERS
zfinput = struct('stage',[],'substage',[],'name',[]','value',[],'sensitivity',[])
zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','levels',1,'low'); %number of scaled levels correlation is performed
zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','iterations',10,'low'); %number of iterations iat is performed
zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','scale',1/16,'low'); %lowered initial scale to increase computation speed for correlation
zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','threshold',0.96,'severe'); %threshold for correlation selection
zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','threshold2',2,'severe'); %threshold for warp modulus

%% Load images (and spot roi)

%FolderPath = uigetdir('/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Data Nynke/')
%FolderPath = uigetdir('C:/Users/260018/Documents/dropbox/');
FolderPath = uigetdir('\\bioinf-filesrv2\cluster15\Ham\Research_group\Data voor Samuel\Input','select folder with tif files');

[imageinfo] = ImageInfoLink2(FolderPath,zfinput);
[stackinfo] = StackInfoLink2(imageinfo);


%% save image, stack and imput info

%create savepath
%[FolderPath1,FolderName] = fileparts(FolderPath)
%SavePath = [fileparts(FolderPath1),'/Output/',FolderName]
%mkdir(SavePath)
%SavePath = uigetdir('/Volumes/Seagate Expansion Drive/ZebraFish/Input/')

save([FolderPath,'/imageinfo.mat'],'imageinfo')
save([FolderPath,'/stackinfo.mat'],'stackinfo')
save([FolderPath,'/zfinput.mat'],'zfinput')

