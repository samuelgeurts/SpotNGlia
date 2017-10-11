%generate artificial imageas based on red channel with random translations
%as the blue en green channel. Als generate an artificial image stack based
%on a single image


% original image path

OriginalPath = uigetdir('/Volumes/Seagate Expansion Drive/ZebraFish/Input/','select loadpath');
PreprocessionPath = uigetdir('/Volumes/Seagate Expansion Drive/ZebraFish/Input/','select loadpath');




load([OriginalPath,'/imageinfo.mat'],'imageinfo')
load([OriginalPath,'/stackinfo.mat'],'stackinfo')
imageinfo = imageinfo([imageinfo.CorNextStack]~=2);

load([PreprocessionPath,'/zfinput.mat'],'zfinput')
pp = zfinput(find(strcmp({zfinput.stage},'preprocession'))) 





%%

%clearvars -except imageinfo stackinfo FolderPath zfinput



%fsv = [fs.value];
for k3 = 1:numel(pp)
    eval([pp(k3).name,'_pp=',num2str(pp(k3).value)])
end

for k20 = 1:10%numel(imageinfo)

        Image = imread([OriginalPath,'/',imageinfo(k20).name]);

        mod = rand(6,1) * 7;
        ang = rand(6,1) * 2 * pi - pi;
        x = sin(ang).*mod;
        y = cos(ang).*mod;
        IniWarp{k20} = [x,y];
         
        %work for older Matlab 2013 version (a little bit different as SmoothEdges cant be used)
        Ref1 = imref2d(size(Image(:,:,1)));
        Img11 = imwarp(Image(:,:,1),affine2d([1 0 0; 0 1 0; x(1) y(1) 1]),'FillValues',0,'Interp','cubic','OutputView',Ref1);      
        Img12 = imwarp(Image(:,:,1),affine2d([1 0 0; 0 1 0; x(2) y(2) 1]),'FillValues',0,'Interp','cubic','OutputView',Ref1); 
        Img21 = imwarp(Image(:,:,2),affine2d([1 0 0; 0 1 0; x(3) y(3) 1]),'FillValues',0,'Interp','cubic','OutputView',Ref1); 
        Img22 = imwarp(Image(:,:,2),affine2d([1 0 0; 0 1 0; x(4) y(4) 1]),'FillValues',0,'Interp','cubic','OutputView',Ref1); 
        Img31 = imwarp(Image(:,:,3),affine2d([1 0 0; 0 1 0; x(5) y(5) 1]),'FillValues',0,'Interp','cubic','OutputView',Ref1); 
        Img32 = imwarp(Image(:,:,3),affine2d([1 0 0; 0 1 0; x(6) y(6) 1]),'FillValues',0,'Interp','cubic','OutputView',Ref1); 
       
        %Img11 = imtranslate(Image(:,:,1),[x(1),y(1)],'FillValues',0,'Method','cubic');              
        %Img12 = imtranslate(Image(:,:,1),[x(2),y(2)],'FillValues',0,'Method','cubic');
        %Img21 = imtranslate(Image(:,:,2),[x(3),y(3)],'FillValues',0,'Method','cubic');
        %Img22 = imtranslate(Image(:,:,2),[x(4),y(4)],'FillValues',0,'Method','cubic');
        %Img31 = imtranslate(Image(:,:,3),[x(5),y(5)],'FillValues',0,'Method','cubic');
        %Img32 = imtranslate(Image(:,:,3),[x(6),y(6)],'FillValues',0,'Method','cubic');
        
        %{
            figure;imagesc(Img11)
        %}
        
        Img1 = cat(3,Image(:,:,1),Img11,Img12);
        Img2 = cat(3,Image(:,:,2),Img21,Img22);
        Img3 = cat(3,Image(:,:,3),Img31,Img32);

        Img1c = Img1(8:end-7,8:end-7,:);
        Img2c = Img2(8:end-7,8:end-7,:);
        Img3c = Img3(8:end-7,8:end-7,:);

        [ECCWarp1, ~] = sng_RGB2(Img1c,scaleC_pp,'translation',levelsC_pp,iterationsC_pp);
        [ECCWarp2, ~] = sng_RGB2(Img2c,scaleC_pp,'translation',levelsC_pp,iterationsC_pp);
        [ECCWarp3, ~] = sng_RGB2(Img3c,scaleC_pp,'translation',levelsC_pp,iterationsC_pp);
        
        ImgFilt1 = sng_RGBselection(Img1c,sigmalp_pp,sigmahp_pp);
        ImgFilt2 = sng_RGBselection(Img2c,sigmalp_pp,sigmahp_pp);
        ImgFilt3 = sng_RGBselection(Img3c,sigmalp_pp,sigmahp_pp);
        
        [ECCWarpFilt1, ~] = sng_RGB2(ImgFilt1,scaleC_pp,'translation',levelsC_pp,iterationsC_pp);
        [ECCWarpFilt2, ~] = sng_RGB2(ImgFilt2,scaleC_pp,'translation',levelsC_pp,iterationsC_pp);
        [ECCWarpFilt3, ~] = sng_RGB2(ImgFilt3,scaleC_pp,'translation',levelsC_pp,iterationsC_pp);            
        
        %fill subtable
        CalWarp{k20} = [ECCWarp1{1}';ECCWarp1{2}';ECCWarp2{1}';ECCWarp2{2}';ECCWarp3{1}';ECCWarp3{2}'];
        CalWarpFilt{k20} = [ECCWarpFilt1{1}';ECCWarpFilt1{2}';ECCWarpFilt2{1}';ECCWarpFilt2{2}';ECCWarpFilt3{1}';ECCWarpFilt3{2}'];
        
end

%%
file = repmat({imageinfo(:).name},[6,1]);file = {file{:}}';
base = repmat({'red';'red';'green';'green';'blue';'blue'},[numel(imageinfo),1]);
InitialWarp = vertcat(IniWarp{:});
ComputedWarp = vertcat(CalWarp{:});
ComputedWarpFilt = vertcat(CalWarpFilt{:});
Difference = abs(InitialWarp-ComputedWarp);        
DifferenceFilt = abs(InitialWarp-ComputedWarpFilt);

RgbAlTestTable = table(InitialWarp,ComputedWarp,ComputedWarpFilt,Difference,DifferenceFilt)     


RgbAlTestTable = table(file,base,InitialWarp,ComputedWarp,ComputedWarpFilt,Difference,DifferenceFilt);       

save([OriginalPath,'/RgbAlignmentTest.mat'],'RgbAlTestTable')




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




        
