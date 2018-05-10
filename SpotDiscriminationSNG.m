
function [SpotsDetected,Spotpar, SpotN] = SpotDiscriminationSNG(Spotpar,cmbr, Reginfo, wp, ZFParameters)


%{
INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters','RegistrationInfo');
Reginfo = INFO.RegistrationInfo{1}
Spotpar = INFO.SpotParameters{1}
cmbr = INFO.BrainSegmentationInfo(fn).BrainEdge
wp = TEMPLATE.Classifier
%}

%{
fn=3
INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters','RegistrationInfo');
Reginfo = INFO.RegistrationInfo{fn}
Spotpar = SpotParameters{fn}
cmbr = obj.checkup(fn).Midbrain
wp = TEMPLATE.Classifier,...




%}

%% retrieve features from SpotParameters(Spotpar) and RegistrationInfo(Reginfo)

    %spot color features
    SpotMeanhsv = reshape([Spotpar.ColorMeanhsv], 3, numel(Spotpar))';
    SpotMean = reshape([Spotpar.ColorMean], 3, numel(Spotpar))';
    %spot features
    Area = [Spotpar.Area]';
    %contrast features
    Azim = [Spotpar.azimuth]';
    Elev = [Spotpar.elevation]';
    R = [Spotpar.r]';
    %fish image feature
    mmc(1:numel(Spotpar),1) = mean(Reginfo(strcmp({Reginfo.name},'MaxFishColor')).value);

    %create feature matrix with discriminative features
    
    DS =[... 
    [SpotMeanhsv(:, 1)], ...
    [SpotMeanhsv(:, 2)], ...
    [SpotMean(:,2)], ...
    [Area], ...
    [Azim], ...
    [Elev], ...
    [R],...
    [mmc]...
    ];

    pra = prdataset(DS(:, 1:8));
    pra = pra * scalem(pra, 'variance'); %scaling ;

    classified = pra * wp;
    classprobability = getdata(classified);
    temp = num2cell(classprobability(:,1)); [Spotpar.ClassProbability] = temp{:};


%% compute selection insite brain

    [rc] = reshape([Spotpar.Centroid], 2, numel(Spotpar))';
    [in, ~] = inpolygon(rc(:, 1), rc(:, 2), cmbr(:, 2), cmbr(:, 1));
    temp = num2cell(in); [Spotpar.Insite] = temp{:};

    
%% threshold spots    

SpotsDetected = Spotpar([Spotpar.Insite] == 1 & ...
    [Spotpar.ClassProbability] >= 0.3);

SpotN = numel(SpotsDetected);



%{

%% other discriminations
cm = [Regions1.ColorProbability] >= MinProbability;
temp = num2cell(cm); [Regions1.MinProbability] = temp{:};


%% determine if insite brain polygon
[rc] = reshape([Regions1.Centroid], 2, numel(Regions1))';
[in, ~] = inpolygon(rc(:, 1), rc(:, 2), cmbr(:, 2), cmbr(:, 1));
temp = num2cell(in); [Regions1.Insite] = temp{:};

%{
    figure;
    imagesc(Ialigned{1})
    hold on
    plot(cmbr(:,2),cmbr(:,1))
    plot(rc(in,1),rc(in,2),'r+')
    %plot(rc(on,1),rc(on,2),'g+')
    plot(rc(~in,1),rc(~in,2),'bo')
 
    figure;hist([Regions1.Area],30)
%}

%% spotsize

%spots larger than MinSpotSize
lt = [Regions1.Area] >= MinSpotSize;
temp = num2cell(lt); [Regions1.LargerThan] = temp{:};

%spots smaller than MaxSpotSize
st = [Regions1.Area] <= MaxSpotSize;
temp = num2cell(st); [Regions1.SmallerThan] = temp{:};

%%
%{
    Ifilt = ismember(L,find(...
        [Regions1.Insite] == 1 &...
        [Regions1.LargerThan] == 1 &...
        [Regions1.SmallerThan] == 1 &...
        [Regions1.MinProbability] == 1));
%}

%nspots = sum([Regions1.Insite] == 1 &...
%    [Regions1.LargerThan] == 1 &...
%    [Regions1.SmallerThan] == 1 &...
%    [Regions1.MinProbability] == 1);


SpotsDetected = Regions1([Regions1.Insite] == 1 & ...
    [Regions1.LargerThan] == 1 & ...
    [Regions1.SmallerThan] == 1 & ...
    [Regions1.MinProbability] == 1);


%{
    figure;imagesc(Ifilt);axis off tight equal
        hold on;plot(cmbr(:,2),cmbr(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
    colormap gray
    set(gca,'position',[0 0 1 1],'units','normalized')
    truesize
    figure;scatter3(colormean(:,1),colormean(:,2),colormean(:,3))
    %}
    
    %Ispots = Ialigned;
    %Ispots(repmat(Ifilt,1,1,3) == 0) = 0;
    
    %{
        figure;imagesc(uint8(Ispots));axis off tight equal
        hold on;plot(cmbr(:,2),cmbr(:,1),'Color',[0,176/255,240/255],'LineWidth',2)
 
        colormap gray
        set(gca,'position',[0 0 1 1],'units','normalized')
        truesize
 
       Ispotshsv =  rgb2hsv(Ispots);
       figure;imagesc(Ispotshsv);axis off tight equal
 
 
       Ihsv1 = Ispotshsv(:,:,1);
       Ihsv2 = Ispotshsv(:,:,2);
       Ihsv3 = Ispotshsv(:,:,3);
        figure;hist(Ihsv1(Ihsv1~=0))
        figure;hist(Ihsv2(Ihsv2~=0))
        figure;hist(Ihsv3(Ihsv3~=0))
 
        Ehsi=sng_rgb2hsi2(colormean')'
        pca1 = pca(colormean','NumComponents',1)
 
        SpotCom = cell2mat({Regions2.Centroid}')
 
        figure;imagesc(uint8(Ispots));axis off tight equal
        hold on;
        scatter(SpotCom(pca1<=0.07,1),SpotCom(pca1<=0.07,2),500,'blue')
        scatter(SpotCom(pca1>=0.25,1),SpotCom(pca1>=0.25,2),500,'yellow')
    %}
    
%}