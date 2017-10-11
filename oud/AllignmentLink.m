function [I80,Parameters] = AllignmentLink(Icombined,CompleteTemplate)
% this function is a subfunction of "Main Function" which alligns a fish
% with a given template in a couple of steps.
% The function is strongly based on "HorizontalFish7"


template = CompleteTemplate.Template;
templatecoords = CompleteTemplate.EyeRegion;

scales = linspace(0.6,1.3,200); %for scale and flip finding
I10 = Icombined;
        
%{
[xcom1,ycom1] = sng_CenterOfMassColor(imcomplement(I10),1);
shift = (fliplr(size(I10(:,:,1))/2)-[xcom1,ycom1]);
I11= imtranslate(I10,shift,'OutputView','same');
figure;imagesc(I11)
%}

%horizontal allignment
%remove background to compute center of mass in horizontalfish well
I20 = sng_RemoveBackgroundColor3(I10,'TriangleSmooth',1,'cuboid');
tform1 = sng_HorizontalFish2(I20,'precise'); %translation + rotation
[I30,ref30] = imwarp(I10,tform1,'FillValues',255);

%scale and leftrigh fitting with template
%cropping is needed to perform well for ScaleFish3 as it is assumed that fishes are centered.
%this crop is not interesting as the image I50 after scaling does not use the cropped image
[I40,ref40] = sng_fishcrop3(I30,ref30,[1000,1360]);
%crop1 = [ref40.XWorldLimits-ref30.XWorldLimits;ref40.YWorldLimits-ref30.YWorldLimits]

tform2 = sng_ScaleFish3(double(template),double(I40),scales,1/16);
    tform_12 = affine2d(tform1.T*tform2.T);
[I50,ref50] = imwarp(I10,tform_12,'FillValues',255);


%allign fish to template. cropping in advance speeds up allignment
%but is it also neccesary?
[I60,ref60] = sng_fishcrop3(I50,ref50,[1024,1360]);

[tform3,I65,CorCoef] = sng_AllignFish2Template2(template,I60,1/4,'translation',3,40);
    tform_123 = affine2d(tform1.T*tform2.T*tform3.T);

[I70,ref70] = imwarp(I10,tform_123,'FillValues',255);
[I80,ref80] = sng_fishcrop3(I70,ref70,[1000,1500]); 


%{
figure;imagesc(I10)
figure;imagesc(ref30.XWorldLimits,ref30.YWorldLimits,uint8(I30));
figure;imagesc(ref40.XWorldLimits,ref40.YWorldLimits,uint8(I40));
figure;imagesc(ref50.XWorldLimits,ref50.YWorldLimits,uint8(I50));
figure;imagesc(ref60.XWorldLimits,ref60.YWorldLimits,uint8(I60));
figure;imagesc(uint8(I65));
figure;imagesc(ref70.XWorldLimits,ref70.YWorldLimits,uint8(I70));
figure;imagesc(uint8(I70))
figure;imagesc(ref75.XWorldLimits,ref75.YWorldLimits,uint8(I75));
figure;imagesc(uint8(I71))
figure;imagesc(uint8(I75))

figure;imagesc(ref80.XWorldLimits,ref80.YWorldLimits,uint8(I80));
figure;imagesc(uint8(template));
figure;imagesc(uint8(CompleteTemplate.MidBrainDistanceMap));
%}



%world coordinatestemplate of template eye after template matching projected on fish (eyes)
xteye1_wc = (ref60.XWorldLimits(1) + templatecoords{1}{1}(:,1) + 0.5);
yteye1_wc = (ref60.YWorldLimits(1) + templatecoords{1}{1}(:,2) + 0.5);
xteye2_wc = (ref60.XWorldLimits(1) + templatecoords{1}{2}(:,1) + 0.5);
yteye2_wc = (ref60.YWorldLimits(1) + templatecoords{1}{2}(:,2) + 0.5);
%{
figure;
imagesc(ref80.XWorldLimits,ref80.YWorldLimits,uint8(I80));colormap(gray);
hold on; plot(xtemplatecoords1_wc,ytemplatecoords1_wc)
hold on; plot(xtemplatecoords2_wc,ytemplatecoords2_wc)
%}
%pixel coordinates of template eyes in final horizontal fish
xteye1_pix = xteye1_wc - ref80.XWorldLimits(1);
yteye1_pix = yteye1_wc - ref80.YWorldLimits(1);
xteye2_pix = xteye2_wc - ref80.XWorldLimits(1);
yteye2_pix = yteye2_wc - ref80.YWorldLimits(1);

%{
figure;
imagesc(uint8(I80));colormap(gray);
hold on; plot(xteye1_pix,yteye1_pix)
hold on; plot(xteye2_pix,yteye2_pix)

%pixel coordinates of template eyes in initial (preprocessed) fish
[xteye1_i,yteye1_i] = transformPointsInverse(tform_123,xteye1_wc,yteye1_wc);
[xteye2_i,yteye2_i] = transformPointsInverse(tform_123,xteye2_wc,yteye2_wc);

figure;
imagesc(uint8(I10));colormap(gray);
hold on; plot(xteye1_i,yteye1_i)
hold on; plot(xteye2_i,yteye2_i)
%}

%add the folowing to template coordinates to acchieve pixel coordinates
temp2fish = [ref60.XWorldLimits(1) + 0.5 - ref80.XWorldLimits(1),...
    ref60.YWorldLimits(1) + 0.5 - ref80.YWorldLimits(1)];
%computes center of mean midbrain projected on fish
cxy = CompleteTemplate.CenterMidBrain + temp2fish
%computes mean midbrain projected on fish 
mmbx = CompleteTemplate.MeanMidBrain(:,1) + temp2fish(1)
mmby = CompleteTemplate.MeanMidBrain(:,2) + temp2fish(2)





%transform to pixel coordinates current fish




Parameters.tform_hor = tform1;
Parameters.tform_scale = tform2;
Parameters.tform_trans = tform3;
Parameters.tform_complete = tform_123;
Parameters.ref70 = ref70;
Parameters.ref80 = ref80;
Parameters.TemplateEyes1 = [xteye1_pix,yteye1_pix];
Parameters.TemplateEyes2 = [xteye2_pix,yteye2_pix];
Parameters.temp2fish = temp2fish;
Parameters.TemplateCenterMidBrain = cxy;
Parameters.TemplateMeanMidBrain = [mmbx,mmby]
Parameters.TemplateCorrelation = CorCoef;




%{
figure
imshowpair(template,uint8(I65))
%}

end


%% Stuff from Horizontalfish7
%{
    %% coordinate calculations    

    
    %roi coordinates per fish head
    if ~isempty(fullpaths2) 
        [xcoord_30,ycoord_30] = transformPointsForward(tform1,coords_head{k}{1}(:,1),coords_head{k}{1}(:,2));
        [xcoord_50,ycoord_50] = transformPointsForward(tform_12,coords_head{k}{1}(:,1),coords_head{k}{1}(:,2));
        [xcoord_70,ycoord_70] = transformPointsForward(tform_123,coords_head{k}{1}(:,1),coords_head{k}{1}(:,2));
    end
    %roi coordinates per fish brain
    if ~isempty(fullpaths3) 
        [xcoord_30b,ycoord_30b] = transformPointsForward(tform1,coords_brain{k}{1}(:,1),coords_brain{k}{1}(:,2));
        [xcoord_50b,ycoord_50b] = transformPointsForward(tform_12,coords_brain{k}{1}(:,1),coords_brain{k}{1}(:,2));
        [xcoord_70b,ycoord_70b] = transformPointsForward(tform_123,coords_brain{k}{1}(:,1),coords_brain{k}{1}(:,2));
    end
    if ~isempty(fullpaths4) 
    %roi coordinates per fish eyes 1
        [xcoord_30c,ycoord_30c] = transformPointsForward(tform1,coords_eyes{k}{1}(:,1),coords_eyes{k}{1}(:,2));
        [xcoord_50c,ycoord_50c] = transformPointsForward(tform_12,coords_eyes{k}{1}(:,1),coords_eyes{k}{1}(:,2));
        [xcoord_70c,ycoord_70c] = transformPointsForward(tform_123,coords_eyes{k}{1}(:,1),coords_eyes{k}{1}(:,2));    
        %roi coordinates per fish eyes 2
        [xcoord_30d,ycoord_30d] = transformPointsForward(tform1,coords_eyes{k}{2}(:,1),coords_eyes{k}{2}(:,2));
        [xcoord_50d,ycoord_50d] = transformPointsForward(tform_12,coords_eyes{k}{2}(:,1),coords_eyes{k}{2}(:,2));
        [xcoord_70d,ycoord_70d] = transformPointsForward(tform_123,coords_eyes{k}{2}(:,1),coords_eyes{k}{2}(:,2));    
    end   
 
    
    %template roi coordinates after template matching projected on fish (eyes)
    xtemplatecoords1_wc = (ref60.XWorldLimits(1) + templatecoords{1}{1}(:,1) + 0.5);
    ytemplatecoords1_wc = (ref60.YWorldLimits(1) + templatecoords{1}{1}(:,2) + 0.5);
    xtemplatecoords2_wc = (ref60.XWorldLimits(1) + templatecoords{1}{2}(:,1) + 0.5);
    ytemplatecoords2_wc = (ref60.YWorldLimits(1) + templatecoords{1}{2}(:,2) + 0.5);
    %to calculate the template coordinates for I60
    %[xtemplatecoords_match,ytemplatecoords_match] = transformPointsInverse(tform_12,xtemplatecoords_wc,ytemplatecoords_wc);    
    
    %roi coordinates of fish to template (head)
    if ~isempty(fullpaths2) 
        xtemplatecoord_60 =  xcoord_50  - ref60.XWorldLimits(1) - 0.5;
        ytemplatecoord_60 =  ycoord_50  - ref60.YWorldLimits(1) - 0.5;
        [xtemplatecoord_65{k},ytemplatecoord_65{k}] = transformPointsForward(tform3,xtemplatecoord_60,ytemplatecoord_60);
    end
    %roi coordinates of fish to template (brain)
    if ~isempty(fullpaths3) 
        xtemplatecoord_60b = xcoord_50b - ref60.XWorldLimits(1) - 0.5;
        ytemplatecoord_60b = ycoord_50b - ref60.YWorldLimits(1) - 0.5;
        [xtemplatecoord_65b{k},ytemplatecoord_65b{k}] = transformPointsForward(tform3,xtemplatecoord_60b,ytemplatecoord_60b);
    end
        
%% save images and coordinates
    if strcmp(SaveImagesAndVariables,'yes')
        fullsavepath = [savepath1,filenames(k).name];
        imwrite(uint8(I80),fullsavepath,...
                'WriteMode', 'overwrite', 'Compression','none');

        baseFileName = sscanf(filenames(k).name, '%c',7);

        if ~isempty(fullpaths2) 
            headcoords = [xcoord_70,ycoord_70];
            headcoordsalligned = [xtemplatecoord_65{k},ytemplatecoord_65{k}];
        else
            headcoords = [];
            headcoordsalligned = [];
        end
        
        if ~isempty(fullpaths3) 
            braincoords = [xcoord_70b,ycoord_70b];
            braincoordsalligned = [xtemplatecoord_65b{k},ytemplatecoord_65b{k}];
        else
            braincoords = [];
            braincoordsalligned = [];
        end
        
        if ~isempty(fullpaths4) 
            eyecoords1 = [xcoord_70c,ycoord_70c];
            eyecoords2 = [xcoord_70d,ycoord_70d];
        else
            eyecoords1 = [];
            eyecoords2 = [];
        end
        
        tempeyecoords1 = [xtemplatecoords1_wc,ytemplatecoords1_wc];
        tempeyecoords2 = [xtemplatecoords2_wc,ytemplatecoords2_wc];

        

        %saves variables (and coordinates)
        if strcmp(TransformCoordinates,'no');
            save([savepath1,baseFileName,'.mat'],'ref60','ref80','tform1','tform_12','tform_123');
        elseif strcmp(TransformCoordinates,'yes');
            save([savepath1,baseFileName,'.mat'],'ref60','ref80','braincoords','headcoords','eyecoords1','eyecoords2'...
            ,'tform1','tform_12','tform_123','tempeyecoords1','tempeyecoords2','braincoordsalligned','headcoordsalligned')
        end   
    end
    

%% Images    

    if strcmp(ShowImages,'yes')
        for j = 8
            figure
            sng_FigureLocation(1.7,'downright')
            hold on

            if j == 1;
                imagesc(I10);
                if ~isempty(fullpaths2) 
                    plot(coords_head{k}{1}(:,1),coords_head{k}{1}(:,2))
                end
                if ~isempty(fullpaths3) 
                    plot(coords_brain{k}{1}(:,1),coords_brain{k}{1}(:,2))
                end
                if ~isempty(fullpaths4) 
                    plot(coords_eyes{k}{1}(:,1),coords_eyes{k}{1}(:,2))
                    plot(coords_eyes{k}{2}(:,1),coords_eyes{k}{2}(:,2))
                end
            elseif j == 2
                imagesc(ref30.XWorldLimits,ref30.YWorldLimits,uint8(I30));colormap(gray);
                if ~isempty(fullpaths2) 
                    plot(xcoord_30,ycoord_30)
                end
                if ~isempty(fullpaths3)
                    plot(xcoord_30b,ycoord_30b)
                end
                if ~isempty(fullpaths4) 
                    plot(xcoord_30c,ycoord_30c)  
                    plot(xcoord_30d,ycoord_30d)
                end
            elseif j == 3;
                imagesc(ref40.XWorldLimits,ref40.YWorldLimits,uint8(I40));colormap(gray);
                if ~isempty(fullpaths2)
                    plot(xcoord_30,ycoord_30)
                end
                if ~isempty(fullpaths3)
                    plot(xcoord_30b,ycoord_30b)
                end
                if ~isempty(fullpaths4) 
                    plot(xcoord_30c,ycoord_30c)
                    plot(xcoord_30d,ycoord_30d)
                end
            elseif j == 4;
                imagesc(ref50.XWorldLimits,ref50.YWorldLimits,uint8(I50));colormap(gray);
                if ~isempty(fullpaths2)
                    plot(xcoord_50,ycoord_50)
                end
                if ~isempty(fullpaths3)
                    plot(xcoord_50b,ycoord_50b)
                end
                if ~isempty(fullpaths4) 
                    plot(xcoord_50c,ycoord_50c)
                    plot(xcoord_50d,ycoord_50d)
                end
            elseif j == 5;
                imagesc(ref60.XWorldLimits,ref60.YWorldLimits,uint8(I60));colormap(gray);
                if ~isempty(fullpaths2)
                    plot(xcoord_50,ycoord_50)
                end
                if ~isempty(fullpaths3)
                    plot(xcoord_50b,ycoord_50b)
                end
                if ~isempty(fullpaths4) 
                    plot(xcoord_50c,ycoord_50c)
                    plot(xcoord_50d,ycoord_50d)
                end
                %plot(xtemplatecoords_match,ytemplatecoords_match)
            elseif j == 6;
                imagesc(ref70.XWorldLimits,ref70.YWorldLimits,uint8(I70));colormap(gray);
                if ~isempty(fullpaths2)
                    plot(xcoord_70,ycoord_70)
                end
                if ~isempty(fullpaths3)
                    plot(xcoord_70b,ycoord_70b)
                end
                if ~isempty(fullpaths4) 
                    plot(xcoord_70c,ycoord_70c)
                    plot(xcoord_70d,ycoord_70d)
                end
                plot(xtemplatecoords1_wc,ytemplatecoords1_wc)
                plot(xtemplatecoords2_wc,ytemplatecoords2_wc)
            elseif j == 7;
                imagesc(ref80.XWorldLimits,ref80.YWorldLimits,uint8(I80));colormap(gray);
                if ~isempty(fullpaths2)
                    plot(xcoord_70,ycoord_70)
                end
                if ~isempty(fullpaths3)
                    plot(xcoord_70b,ycoord_70b)
                end
                if ~isempty(fullpaths4) 
                    plot(xcoord_70c,ycoord_70c)
                    plot(xcoord_70d,ycoord_70d)
                end
                plot(xtemplatecoords1_wc,ytemplatecoords1_wc)
                plot(xtemplatecoords2_wc,ytemplatecoords2_wc)
            elseif j == 8;
                imshowpair(template,uint8(I65))
            end
            
            axis ij equal %off tight
            hold off   
        end
        drawnow
    end

 %%       
end

if strcmp(ShowTemplate,'yes') && ~isempty(fullpaths3)
    figure; imshowpair(template,uint8(I65))
    figure
    sng_FigureLocation(1.7,'downright')
    hold on
    imagesc(uint8(template));colormap(gray);
    for l = 1:numel(fullpaths)    
        %plot(xtemplatecoord_65{l},ytemplatecoord_65{l})
        plot(xtemplatecoord_65b{l},ytemplatecoord_65b{l}) 
    end
    hold off
end

if strcmp(ShowTemplate,'yes') || strcmp(ShowImages,'yes')
    sng_figureslide3
end

%{
figure;hold on
    imagesc(uint8(I65));colormap(gray);
    plot(xtemplatecoord_65b{l},ytemplatecoord_65b{l}) 
hold off
%}    




%}


