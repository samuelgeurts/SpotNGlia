
INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters','RegistrationInfo');


for k1 = 1:numel(INFO.SpotParameters)

    %retrieve featrues from SpotParameters and RegistrationInfo
    
    Spotpar = INFO.SpotParameters{1}
    %spot color features
    SpotMeanhsv = reshape([Spotpar.ColorMeanhsv], 3, numel(Spotpar))';
    SpotMean = reshape([Spotpar.ColorMean], 3, numel(Spotpar))';
    %spot features
    Area = [Spotpar.Area]'
    %contrast features
    Azim = [Spotpar.azimuth]'
    Elev = [Spotpar.elevation]'
    R = [Spotpar.r]'
    %fish image feature
    mmc(1:numel(Spotpar),1) = mean(INFO.RegistrationInfo{k1}(strcmp({INFO.RegistrationInfo{k1}.name},'MaxFishColor')).value);

    %create feature matrix
    
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

lab0 = zeros(numel(Spotpar),1)

a = prdataset(DS(:, 1:8),lab0);
a = a * scalem(a, 'variance'); %scaling 

classified = a * wp
getdata(classified)





