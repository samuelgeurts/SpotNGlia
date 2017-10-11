%% setting parameters and image check

%clear all;close all

savefig1_TF = false;
savefig2_TF = false;


if ~exist('FolderPath')
    disp('select loadpath')
    %Load images (and spot roi)
    %[Names,FolderPath] = uigetfile('/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Data Nynke/','MultiSelect','on')
    %FolderPath = uigetdir('/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Data Nynke/')
    %FolderPath = uigetdir('C:/Users/260018/Documents/dropbox/');
    %FolderPath = uigetdir('\\bioinf-filesrv2\cluster15\Ham\Research_group\Data voor Samuel\Input')
    FolderPath = uigetdir('/Volumes/Seagate Expansion Drive/ZebraFish/Input/','select folder with original tif files (input)');
end

load([FolderPath,'/imageinfo.mat'],'imageinfo')
load([FolderPath,'/stackinfo.mat'],'stackinfo')
load([FolderPath,'/zfinput.mat'],'zfinput')

%parameters
zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','onoff',false,''); %sigma for bandpassfilter (lowpass)    
zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','sigmalp',1,'moderate'); %sigma for bandpassfilter (lowpass)
zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','sigmahp',4,'moderate'); %for bandpassfilter (highpass)
zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','scaleC',1/4,'low'); %lowered initial scale to increase computation speed for correlation
zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','levelsC',2,'low'); %number of scaled levels correlation is performed
zfinput = sng_zfinput(zfinput,0,'preprocession','RGB correction','iterationsC',10,'low'); %number of iterations iat is performed
zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','scaleS',1/4,'low'); %lowered initial scale to increase computation speed for correlation
zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','levelsS',3,'low'); %number of scaled levels correlation is performed
zfinput = sng_zfinput(zfinput,0,'preprocession','slice stack correction','iterationsS',20,'low'); %iterations iat is performed

open imageinfo
open stackinfo


button = questdlg({'Check slicecombination in imageinfo and stackinfo. Apply corrections in "imageinfo.CorNextStack". Set value "2" for removing images from stack. Did you apply corrections?'})

if button == 'Yes'

    [stackinfo] = StackInfoLink2(imageinfo);

    disp('select savepath')
    SavePath = uigetdir('/Volumes/Seagate Expansion Drive/ZebraFish/Input/','select location to save preprocessed images (output)')

    %% preprocession
    for k1 = 1:numel(stackinfo)
    disp(k1)
        %select slices
        ImageSlice = cell(1,stackinfo(k1).stacksize);%preallocate for every new file
        for k2 = 1:stackinfo(k1).stacksize
            ImageSlice{k2} = imread([FolderPath,'/',stackinfo(k1).imagenames{k2}]);
        end
        
        if savefig1_TF
            sng_SaveCell2TiffStack(ImageSlice,[SavePath,'/stack/',stackinfo(k1).stackname,'.tif'])
        end
            
        [CorrectedSlice,ppoutput] = PreprocessionLink(ImageSlice,zfinput);

        if savefig2_TF
            sng_SaveCell2TiffStack(CorrectedSlice,[SavePath,'/',stackinfo(k1).stackname,'.tif'])
        end
        
        stackinfo(k1).Preprocession = ppoutput;
    end

    save([SavePath,'/stackinfo.mat'],'stackinfo')
    save([SavePath,'/imageinfo.mat'],'imageinfo')
    save([SavePath,'/input.mat'],'zfinput')
end
