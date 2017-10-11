%generate artificial imageas based on red channel with random translations
%as the blue en green channel. Als generate an artificial image stack based
%on a single image


if ~exist('FolderPath')
    disp('select loadpath')
    FolderPath = uigetdir('/Volumes/Seagate Expansion Drive/ZebraFish/Input/');
end

load([FolderPath,'/imageinfo.mat'],'imageinfo')
load([FolderPath,'/stackinfo.mat'],'stackinfo')
load([FolderPath,'/input.mat'],'zfinput')
%select only ok images
imageinfo = imageinfo([imageinfo.nextstack]~=2);

clearvars -except imageinfo stackinfo FolderPath zfinput

pp = zfinput(find(strcmp({zfinput.stage},'preprocession'))) 
%fsv = [fs.value];
for k3 = 1:numel(pp)
    eval([pp(k3).name,'_pp=',num2str(pp(k3).value)])
end

SliceAlTestTable = table;
cell = [2;3;4;5];



%%





for k20 = 1:numel(imageinfo)
    
    mod = rand(4,1) * 60;
    ang = rand(4,1) * 2 * pi - pi; 
    x = sin(ang) .* mod;
    y = cos(ang) .* mod;
    IniWarp = [x,y];
    
    Image = imread([FolderPath,'/',imageinfo(k20).name]);

    Imagestack{1} = Image(61:end-60,61:end-60,1:3);
    for k21 = 2:5
        Img_slice = imtranslate(Image(:,:,1:3),IniWarp(k21-1,:),'FillValues',0,'Method','cubic');
        Imagestack{k21} = Img_slice(61:end-60,61:end-60,1:3);  
        %{
        figure;imagesc(imagestack{k21})
        %}
    end    

    [ECCWarp2,CorrectedSlice] = sng_AllignFishBatch3(Imagestack,scaleS_pp,'translation',levelsS_pp,iterationsS_pp);  

    %fill subtable
    [file{1:4,1}] = deal(imageinfo(k20).name);
    CalWarp = [ECCWarp2{:}]';

    T = table(file,cell,IniWarp,CalWarp);        
    SliceAlTestTable = [SliceAlTestTable;T];
end        

SliceAlTestTable.Difference = abs(SliceAlTestTable.IniWarp-SliceAlTestTable.CalWarp)


save([FolderPath,'/SliceAlignmentTest2.mat'],'SliceAlTestTable')





%%      
    %{
        figure;imagesc(Image(:,:,1:3))
        figure;imagesc(Image(:,:,1))        
        figure;imagesc(Img1)
        figure;imagesc(Img2)
        figure;imagesc(Img3)
        figure;imagesc(ImgFilt1)
        figure;imagesc(ImgFilt2)
        figure;imagesc(ImgFilt3)

        [ImgRGBcor1,crop] = sng_RGB_IATwarp2(Img1c,ECCWarp1);
        [ImgRGBcorF1,crop] = sng_RGB_IATwarp2(Img1c,ECCWarpFilt1);
        [ImgRGBcor2,crop] = sng_RGB_IATwarp2(Img2c,ECCWarp2);
        [ImgRGBcorF2,crop] = sng_RGB_IATwarp2(Img2c,ECCWarpFilt2);
        [ImgRGBcor3,crop] = sng_RGB_IATwarp2(Img3c,ECCWarp3);
        [ImgRGBcorF3,crop] = sng_RGB_IATwarp2(Img3c,ECCWarpFilt3);

        figure;imagesc(ImgRGBcor1)
        figure;imagesc(ImgRGBcorF1)
        figure;imagesc(ImgRGBcor2)
        figure;imagesc(ImgRGBcorF2)
        figure;imagesc(ImgRGBcor3)
        figure;imagesc(ImgRGBcorF3)
    %}  





