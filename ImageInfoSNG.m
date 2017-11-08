function [imageinfo] = ImageInfoSNG(imageinfo,zfinput)
%this function measures the correlation af adjacent fish-images in a folder en
%give the boundaries of batches of fishes that are similar

%Version ImageInfoLink2 20170616
%  add CorNextStack for corrections on nextstack , added warning row for
%  ambigious stackcuts.
% Version ImageInfoLink3 20171026
% remove questionbox and variable to zfinput

%% parameters


if ~exist('zfinput','var')
    zfinput = struct;
    zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','sorting','Date','high');   
    zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','levels',1,'low'); %number of scaled levels correlation is performed
    zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','iterations',10,'low'); %number of iterations iat is performed
    zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','scale',1/16,'low'); %lowered initial scale to increase computation speed for correlation
    zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','threshold',0.96,'severe'); %threshold for correlation selection
    zfinput = sng_zfinput(zfinput,0,'imageinfo','stackselection','threshold2',2,'severe'); %threshold for warp modulus
end

sng_zfinputAssign(zfinput,'imageinfo')



%% imageinfo output variable

if strcmp(sorting,'Date')
    [~, newindex] = sort({imageinfo.date});
    imageinfo = imageinfo(newindex);
end

%{
using the histogram of the images does not lead to very reliable results.
for k = 1:numel(ImagePath)
    ImageSlice = imread(ImagePath{k});
    H(k,:) = histcounts(ImageSlice,0:255);
end

for k = 1:numel(ImagePath)-1
    diff(k) = sum(((H(k,:) - H(k+1,:)).^2))
end
figure;plot(diff)
%}

%% correlation

IS = cell(numel(imageinfo),1);

par.transform = 'translation';
par.levels = levels;
par.iterations = iterations; %iterations per level

ni = numel(imageinfo);

bit = zeros(ni,1);
histmax = zeros(ni,3);
warp = zeros(ni,2);
corcoef = zeros(ni,1);
%fishnumber = zeros(ni,1);
modulusw = zeros(ni,1);
anglew = zeros(ni,1);
warning = cell(ni,1);
nextstack = zeros(ni,1);


for k = 1:ni
    slice = imread([imageinfo(k).folder,'/',imageinfo(k).name]);
    
    %sometimes images are 16 bit
    slice =im2uint8(slice);
    
  
    % store bit info
    info = imfinfo([imageinfo(k).folder,'/',imageinfo(k).name]);
    bit(k) = info.BitsPerSample(1);

    % store histogram peak info to determine if it is oversaturated
    %[~,histmax(k,1)] = max(histcounts(slice(:,:,1),linspace(0,255,256)));
    %[~,histmax(k,2)] = max(histcounts(slice(:,:,2),linspace(0,255,256)));
    %[~,histmax(k,3)] = max(histcounts(slice(:,:,3),linspace(0,255,256)));
    %for older matlab versions:
    s1 = slice(:,:,1);
    s2 = slice(:,:,2);
    s3 = slice(:,:,3);
    
    [~,histmax(k,1)] = max(hist(s1(:),linspace(0,255,256)));
    [~,histmax(k,2)] = max(hist(s2(:),linspace(0,255,256)));
    [~,histmax(k,3)] = max(hist(s3(:),linspace(0,255,256)));

    histmax(k,:) = histmax(k,:) - 1; %as the values go from 0-255 
        
    % compute correlation and warp and add it to imageinfo
    IS{k} = imcomplement(imresize(slice,scale));
    if k>=2
        [warp(k,:),corcoef(k)] = ...
            iat_ecc(IS{k-1}(:,:,1:3), IS{k}(:,:,1:3), par);
        % using the correlation between fishes works but alligning them
        % first works better.
        %l = corrcoef(double(IS{k-1}(:,:,1:3)),double(IS{k}(:,:,1:3)))
        %cc(k) = l(1,2)
        %n(k) = norm(ECCWarp{k}')
    end
    fprintf('\b ');disp(num2str(k));
    
    modulusw(k) = sqrt(warp(k,1).^2+warp(k,2).^2);%modulus warp
    anglew(k) = angle(warp(k,:)*[1;1i]);%angle warp
    

    
   

    %compute fish number based on given threshold and add the number to imageinfo    
    if corcoef(k) <= threshold || isnan(corcoef(k)) || modulusw(k) >= threshold2
        nextstack(k) = 1;
    else
        nextstack(k) = 0;
    end
    
    %warning for ambigious stack-cuts
    if corcoef(k) >= 0.96 && corcoef(k) <= 0.99 && modulusw(k) <=3 
        warning{k} = 'Warn!';
    elseif corcoef(k) >= 0.97 && modulusw(k) >=1 
        warning{k} = 'Warn!';
    elseif bit(k) ~= 8
        warning{k} = 'Wrong Res!';
    elseif (anglew(k) >= -0.18 || anglew(k) <= -0.53) && nextstack(k) == 0
        warning{k} = 'Angle warn';
    elseif (anglew(k) <= -0.18 && anglew(k) >= -0.53) && nextstack(k) == 1
        warning{k} = 'Angle warn';
    else
        warning{k} = '';
    end    
    
  
end


%add variables to imageinfo struct

tmp = num2cell(warp,2); [imageinfo.warp] = tmp{:};
tmp = num2cell(modulusw); [imageinfo.moduluswarp] = tmp{:};
tmp = num2cell(anglew); [imageinfo.anglewarp] = tmp{:};
tmp = num2cell(corcoef); [imageinfo.corcoef] = tmp{:};
tmp = num2cell(nextstack); [imageinfo.nextstack] = tmp{:};
tmp = {imageinfo.nextstack}';[imageinfo.CorNextStack] = tmp{:};
[imageinfo.warning] = warning{:};
tmp = num2cell(bit); [imageinfo.bit] = tmp{:};
tmp = num2cell(histmax,2); [imageinfo.histmax] = tmp{:};
%tmp = num2cell(fishnumber); [imageinfo.fishnumber] = tmp{:};
















%{
figure;plot(q);
figure;histogram(q,linspace(0.9,1,200));
%it can be seen that a good boundary is a correlation coefficeint of 0.96
for images with 1/16 scaling
%}


%{
figure;plot(n)
histogram(n,linspace(0.9,1,200))
figure;plot(cc)
histogram(cc,linspace(0.9,1,200))
%}


%{
%check fish 50 and 51 which gives NaN output do to bad hessian matrix
%it means that when NaN is the output, the fishes differs a lot so a new
%stack is found
f2=IS{50}
g2=IS{51}
figure;imagesc(f2(:,:,1:3))
figure;imagesc(g2(:,:,1:3))
k=51
[a,b] = iat_ecc(IS{k-1}(:,:,1:3), IS{k}(:,:,1:3), par);

%}


%TRIED as measures for discrimination fishes
%histogram correlation of images
%correlation of images
%alligning-distance
%correlation coefficient of alligend fishes works the best

end
