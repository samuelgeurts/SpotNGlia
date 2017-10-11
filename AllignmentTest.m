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

clearvars -except imageinfo stackinfo FolderPath zfinput

pp = zfinput(find(strcmp({zfinput.stage},'preprocession'))) 
%fsv = [fs.value];
for k3 = 1:numel(pp)
    eval([pp(k3).name,'_pp=',num2str(pp(k3).value)])
end

T2 = table;
base = {'red';'red';'green';'green';'blue';'blue'};
channel = {'gr-red';'bl-red';'red-gr';'bl-gr';'re-bl';'gr-bl'};



%%
for k20 = 1:numel(imageinfo)

        Image = imread([FolderPath,'/',imageinfo(k20).name]);

        mod = rand(6,1) * 7;
        ang = rand(6,1) * 2 * pi - pi;
        x = sin(ang).*mod;
        y = cos(ang).*mod;
        IniWarp = [x,y];

        Img11 = imtranslate(Image(:,:,1),[x(1),y(1)],'FillValues',0,'Method','cubic');              
        Img12 = imtranslate(Image(:,:,1),[x(2),y(2)],'FillValues',0,'Method','cubic');
        Img21 = imtranslate(Image(:,:,2),[x(3),y(3)],'FillValues',0,'Method','cubic');
        Img22 = imtranslate(Image(:,:,2),[x(4),y(4)],'FillValues',0,'Method','cubic');
        Img31 = imtranslate(Image(:,:,3),[x(5),y(5)],'FillValues',0,'Method','cubic');
        Img32 = imtranslate(Image(:,:,3),[x(6),y(6)],'FillValues',0,'Method','cubic');

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
        [file{1:6,1}] = deal(imageinfo(k20).name);
        CalWarp = [ECCWarp1{1}';ECCWarp1{2}';ECCWarp2{1}';ECCWarp2{2}';ECCWarp3{1}';ECCWarp3{2}'];
        CalWarpFilt = [ECCWarpFilt1{1}';ECCWarpFilt1{2}';ECCWarpFilt2{1}';ECCWarpFilt2{2}';ECCWarpFilt3{1}';ECCWarpFilt3{2}'];
        
        T = table(file,base,channel,IniWarp,CalWarp,CalWarpFilt);        
        T2 = [T2;T];
end        
        
save([Folderpath,'/RgbAllginmentTest.mat'],'T2')



dif1 = abs(T2.IniWarp-T2.CalWarp)
mean(dif1(:))
std(dif1(:))
        
dif2 = abs(T2.IniWarp-T2.CalWarpFilt)
mean(dif2(:))
std(dif2(:))
        

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


%{


%%
        imagestack{1} = Image(61:end-60,61:end-60,:);
        for k21 = 2:5
            mod3 = rand(1) * 60;
            ang3 = rand(1) * 2 * pi - pi; 
            x3 = sin(ang3)*mod3;
            y3 = cos(ang3)*mod3;
            Img_slice = imtranslate(Image(:,:,1:3),[x3,y3],'FillValues',0,'Method','cubic');
            imagestack{k21} = Img_slice(61:end-60,61:end-60,:);  
            %{
            figure;imagesc(imagestack{k21})
            %}
        end
%%      
        %{
        figure;imagesc(imagestack{1})
        %}



n10 = numel(stackinfo)
n = 1;m=1
for k10 = 1:n10
    for k11 = 1:stackinfo(k10).stacksize
        
        
        
        
        
        
    FilteredImage = sng_RGBselection(ImageSlice{k4}(:,:,1:3),sigmalp_pp,sigmahp_pp);
            
    [ECCWarp{k4}, ~]= sng_RGB2(FilteredImage(:,:,1:3),scaleC_pp,'translation',levelsC_pp,iterationsC_pp);

    [ImageSliceCor{k4},crop{k4}] = sng_RGB_IATwarp2(ImageSlice{k4}(:,:,1:3),ECCWarp{k4});

    
    [FilteredImage2,~] = sng_RGB_IATwarp2(FilteredImage,ECCWarp{k4});
        
        
        WGR(n,:) = stackinfo(k10).Preprocession.ColorWarp{k11}{1};
        WBR(n,:) = stackinfo(k10).Preprocession.ColorWarp{k11}{2};
        %if k11 ~= 1
        %    WS(n,:) = stackinfo(k10).Preprocession.SliceWarp{k11};
        %end
        n = n + 1;
    end
    
    
    for k12 = 2:stackinfo(k10).stacksize
        WS(m,:) = stackinfo(k10).Preprocession.SliceWarp{k12};
        m = m + 1
    end
end

%}