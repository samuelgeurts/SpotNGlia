%creates images of fishes for chapter preproccesion testing of the report 


export1_TF = 0


% Load Paths
PreprocessionPathNF = [Basepath,'/1_preprocessed_nofilt']



%stackpp = sng_openimstack2([PreprocessionPath,'/set3-Image(022-025).tif'])

%% Comparison of Color Spot Correction with and without filter

Im{1} = imread([OriginalPath,'/set1-Image004.tif']);
Im{2} = imread([PreprocessionPath,'/set1-Image(004-008).tif']);
Im{3} = imread([PreprocessionPathNF,'/set1-Image(004-008).tif']);

Im{4} = imread([OriginalPath,'/set3-Image022.tif']);
Im{5} = imread([PreprocessionPath,'/set3-Image(022-025).tif']);
Im{6} = imread([PreprocessionPathNF,'/set3-Image(022-025).tif']);

Im{7} = imread([OriginalPath,'/set5-Image0000.tif']);
Im{8} = imread([PreprocessionPath,'/set5-Image(0000-0002).tif']);
Im{9} = imread([PreprocessionPathNF,'/set5-Image(0000-0002).tif']);

%{
    for k = 7:9
    figure;imagesc(Im{k}(:,:,1:3))
    end
%}

%pixelsize of image
psiz = [40,50]
%startpoint of image, top left corner
rx = [730 681 679 685 647 647 860 842 842]
ry = [528 527 527 522 522 520 577 576 576]

for k = 1:9
    h{k} = figure;
    imagesc(Im{k}(ry(k):ry(k)+psiz(1)-1,rx(k):rx(k)+psiz(2)-1,1:3))
    set(gca, 'Units', 'normalized', 'Position', [0 0 1 1],'Visible','off');
    sng_figcm(5,4,113.6)
    
    if export1_TF
        export_fig(h{k} ,['/Users/samuelgeurts/Desktop/','fig',num2str(k)], '-png', '-r600', '-nocrop');
    end    
end


