function [tform,Img2,CorrCoef,supportECC] = sng_AllignFish2Template2(template,Img,scale,transform,levels,iterations,initialization)
%allign one fish with a template 

%Version: sng_AllignFish2Template2 
%       Add final correlation parameter between template and fish    

%example
%[tform4] = sng_AllignFish2Template(template,I60,1/4,'translation',3,40);
%{
Img=I55;
scale = 1/4
transform = 'translation'
levels = 3
iterations = 40
initialization = [-600/4;0]

transform = 'affine'
initialization = [eye(2),[-600/4;0]]

%}

%algorithm parameters
par.initwarp = initialization;
par.transform = transform;
par.levels = levels;
par.iterations = iterations; %iterations per level

tempsc = imresize(template,scale); %resize to speed up
Imgsc = imresize(Img,scale);


[ECCWarp, CorrCoef]= iat_ecc(Imgsc, tempsc, par);
%ECCWarp= iat_LucasKanade(Imgsc, tempsc, par);


%{

[ECCWarp, CorrCoef]= iat_ecc(tempsc,Imgsc, par);
[wimageECC, supportECC] = iat_inverse_warping(tempsc, ECCWarp, par.transform, 1:size(Imgsc,2),1:size(Imgsc,1));


[wimageECC, supportECC] = iat_inverse_warping(Imgsc, ECCWarp, par.transform, 1:size(tempsc,2),1:size(tempsc,1));
figure;imagesc(uint8(wimageECC))
figure;imshowpair(tempsc,wimageECC)
%}

%correct for scaling. Separated because of parfor loop
if strcmp(par.transform,'translation')
    ECCWarp = ECCWarp*(1/scale);
    T = [1 0 0;0 1 0;-ECCWarp(1),-ECCWarp(2),1];
end

if strcmp(par.transform,'affine')
    ECCWarp(1:2,3) = ECCWarp(1:2,3) * (1/scale);
    T = inv([ECCWarp;[0,0,1]]');
    T(1,3) = 0;
    T(2,3) = 0;
    T(3,3) = 1
end
%{
figure;imagesc(Img)
figure;imagesc(Imgsc)
figure;imagesc(tempsc)



%}


if nargout >= 2
    [M,N,Oh] = size(template);
    [Img2, supportECC] = iat_inverse_warping(Img, ECCWarp, par.transform, 1:N, 1:M);
    %setbackground to white
    Img2 = Img2 + repmat(255*(1-supportECC),1,1,3); 
end
    
%{
xtemplatecoords_wc = (ref60.XWorldLimits(1) + templatecoords{1}{1}(:,1) + 0.5);
ytemplatecoords_wc = (ref60.XWorldLimits(1) + templatecoords{1}{1}(:,2) + 0.5);


figure; imagesc(template)
figure; imagesc(Img)
figure; imshowpair(template,I65)
figure; imagesc(supportECC)
figure; imagesc(uint8(Img2))
figure;imshowpair(template,Img2)

figure;imagesc(uint8(Img2));
hold on
plot(templatecoords{1}{1}(:,1),templatecoords{1}{1}(:,2))
plot(templatecoords{1}{2}(:,1),templatecoords{1}{2}(:,2))
axis equal
hold off

imagesc(ref70.XWorldLimits,ref70.YWorldLimits,uint8(I70));colormap(gray);
plot(xcoord_all{k},ycoord_all{k})
plot(xtemplatecoords_match,ytemplatecoords_match)







%}

tform = affine2d(T);

% %i dont know if this is right
% if strcmp(par.transform,'affine')
%     ECCWarp(1,2) = ECCWarp(1,2)*(1/scale);
%     ECCWarp(2,1) = ECCWarp(2,1)*(1/scale);
% end    
    
end
