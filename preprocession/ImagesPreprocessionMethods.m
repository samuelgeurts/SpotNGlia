%images for report chapter preprocessiong methods
PathsComplete('or','bp','pp')

%export figures to desktop
export1_TF = 0

%generate images
colorchannelcorrection_TF = 0
stackcorrection_TF = 0
extendeddepthoffield_TF = 1


Im2{1} = imread([OriginalPath,'/set2-Image053.tif']);
Im2{2} = imread([OriginalPath,'/set2-Image054.tif']);
Im2{3} = imread([OriginalPath,'/set2-Image055.tif']);
Im2{4} = imread([OriginalPath,'/set2-Image056.tif']);
Im2{5} = imread([OriginalPath,'/set2-Image057.tif']);


CorrectedSlice = sng_openimstack2([PreprocessionPath,'/','set2-Image(053-057).tif']);       
%Icombined = sng_openimstack2([PreprocessionPath,'/','set2-Image(053-057)-ExtendedDeptOfField.tif']);       



%%


%input variables        
load([PreprocessionPath,'/zfinput.mat'],'zfinput')    
%compute color correction and stack correction




%% Color Channel Correction
if colorchannelcorrection_TF

    zfinput = sng_zfinput(zfinput,0,'preprocession','bandpass filtering','onoff',true,'');
    [CorrectedSlice,ppoutput,ColorCorrect,FiltIm] = PreprocessionLink(Im2,zfinput);
    
%compare channels - unfinished

%bandpass filtered - unfinished
    figure;imagesc(0.15*FiltIm{1}(1:end/2,:,:)) 
    figure;imagesc(0.15*FiltIm{1}(end/2+2:end,:,:)) 

%spot compare color correction
    psiz3= [40,50];
    fsx = 7; %widht in cm 
    fsy = fsx*(4/5); %height in cm

    rx2 = 782
    ry2 = 623
    h4 = figure;
    imagesc(uint8(CorrectedSlice{1}(ry2:ry2+psiz3(1)-1,rx2:rx2+psiz3(2)-1,1:3)))
    set(gca, 'Units', 'normalized', 'Position', [0 0 1 1],'Visible','off');
    sng_figcm(fsx,fsy)
    %
    rx2 = 835
    ry2 = 623
    h5 = figure;
    imagesc(uint8(Im2{1}(ry2:ry2+psiz3(1)-1,rx2:rx2+psiz3(2)-1,1:3)))
    %
    set(gca, 'Units', 'normalized', 'Position', [0 0 1 1],'Visible','off');
    sng_figcm(fsx,fsy)
    %
    if export1_TF
        export_fig(h4 ,['/Users/samuelgeurts/Desktop/','spotrgb',num2str(1)], '-png', '-r600', '-nocrop');
        export_fig(h5 ,['/Users/samuelgeurts/Desktop/','spotrgb',num2str(2)], '-png', '-r600', '-nocrop');
    end    

end
%% stack only color corrected
if stackcorrection_TF
    

%uncorrected stack
    wt = 6;
    for k = 1:5
        sz = size(ColorCorrect{k}); %to lock acpect ratio
        rt = sz(2)/sz(1)   
        h2{k} = figure;
        imagesc(ColorCorrect{k})
        set(gca, 'Units', 'normalized', 'Position', [0 0 1 1],'Visible','off');
        sng_figcm(wt,wt/rt)
        if export1_TF
            export_fig(h2{k} ,['/Users/samuelgeurts/Desktop/','fig',num2str(k)], '-png', '-r600', '-nocrop');
        end    
    end
    
%zoom
    psiz2 = [120,120];
    ht2 = 3;
    rx2 = 780;
    ry2 = 610
    for k = 1:5   
        h2{k} = figure;
        imagesc(ColorCorrect{k}(ry2:ry2+psiz2(1)-1,rx2:rx2+psiz2(2)-1,1:3))
        set(gca, 'Units', 'normalized', 'Position', [0 0 1 1],'Visible','off');
        sng_figcm(ht2,ht2)
        %sng_figcm(5,4,113.6)
        if export1_TF
            export_fig(h2{k} ,['/Users/samuelgeurts/Desktop/','figzoom',num2str(k)], '-png', '-r600', '-nocrop');
        end    
    end

%Stack Correction example Corrected
    wt = 6;
    for k = 1:5       
        sz = size(ColorCorrect{k});
        rt = sz(2)/sz(1)
        h2{k} = figure;
        imagesc(uint8(CorrectedSlice{k}))
        set(gca, 'Units', 'normalized', 'Position', [0 0 1 1],'Visible','off');
        sng_figcm(wt,wt/rt)
        %sng_figcm(5,4,113.6)
        if export1_TF
            export_fig(h2{k} ,['/Users/samuelgeurts/Desktop/','figcor',num2str(k)], '-png', '-r600', '-nocrop');
        end    
    end

%zoom 
    rx2 = 732
    ry2 = 610
    for k = 1:5
        h2{k} = figure;
        imagesc(uint8(CorrectedSlice{k}(ry2:ry2+psiz2(1)-1,rx2:rx2+psiz2(2)-1,1:3)))
        set(gca, 'Units', 'normalized', 'Position', [0 0 1 1],'Visible','off');
        sng_figcm(ht2,ht2)
        if export1_TF
            export_fig(h2{k} ,['/Users/samuelgeurts/Desktop/','figcorzoom',num2str(k)], '-png', '-r600', '-nocrop');
        end    
    end

end
%% Extended Dept of Field
if extendeddepthoffield_TF

    zfinput = sng_zfinput(zfinput,0,'ExtendedDeptOfField','edof','variancedisksize',7,'moderate'); %sigma for bandpassfilter (lowpass)    
    [Icombined,edoutput] = ExtendedDeptofFieldLink2(CorrectedSlice,zfinput);

% combined fish
    wt = 6;
    sz = size(Icombined);
    rt = sz(2)/sz(1)   
    h3 = figure;
    imagesc(Icombined)
    set(gca, 'Units', 'normalized', 'Position', [0 0 1 1],'Visible','off');
    sng_figcm(wt,wt/rt)
    if export1_TF
        export_fig(h3,[Basepath,'/','edof',num2str(1)], '-png', '-r600', '-nocrop');
    end    
       
    rx2 = 732
    ry2 = 610
    ht2 = 3;
    psiz2 = [120,120]

    h3b = figure;
    imagesc(uint8(Icombined(ry2:ry2+psiz2(1)-1,rx2:rx2+psiz2(2)-1,1:3)))
    set(gca, 'Units', 'normalized', 'Position', [0 0 1 1],'Visible','off');
    sng_figcm(ht2,ht2)
    if export1_TF
        export_fig(h3b,[Basepath,'/','edof',num2str(2)], '-png', '-r600', '-nocrop');
    end    
end



%}
