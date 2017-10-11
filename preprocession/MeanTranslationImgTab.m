%Displays tables of mean and std of computed translations with and without
%filter color correction and images. Also stack correcting tables en
%images.


fsx = 8.5; %widht in cm 
fsy = 8.5; %height in cm
Color1 = [0 0 0.7];
Color2 = [0 0.7 0 ];
Color3 = [0 166/255 214/255];

PathsComplete('bp','pp')
PreprocessionPathNF = [Basepath,'/1_preprocessed_nofilt']


for k0 = 1:2
    
    if k0 == 1
 
        load([PreprocessionPathNF,'/stackinfo.mat'],'stackinfo')
    elseif k0 == 2
        load([PreprocessionPath,'/stackinfo.mat'],'stackinfo')
    end
    %applied in report on preprocced with filter and without filter

    %% measure the warp
    %clearvars -except imageinfo stackinfo FolderPath

    n10 = numel(stackinfo)
    n = 1;m=1
    for k10 = 1:n10
        for k11 = 1:stackinfo(k10).stacksize        
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

    %WGR = Warp Green-Red Channel
    meanWGR = mean(WGR);
    stdWGR = std(WGR);
    maxWGR = max(WGR);
    minWGR = min(WGR);
    maxstdWGR = meanWGR + 3 * stdWGR
    minstdWGR = meanWGR - 3 * stdWGR
    modWGR = sqrt(WGR(:,1).^2+WGR(:,2).^2);
    rangemodWGR = mean(modWGR) + 3 * std(modWGR) * [-1 1]
    angWGR = atan2(WGR(:,1),WGR(:,2))
    angrangeWGR = mean(angWGR) + 3 * std(angWGR) * [-1 1]

    %WBR = Warp Blue-Red Channel
    meanWBR = mean(WBR);
    stdWBR = std(WBR);
    maxWBR = max(WBR);
    minWBR = min(WBR);
    maxstdWBR = meanWBR + 3 * stdWBR;
    minstdWBR = meanWBR - 3 * stdWBR;
    modWBR = sqrt(WBR(:,1).^2+WBR(:,2).^2);
    rangemodWBR = mean(modWBR) + 3 * std(modWBR) * [-1 1];
    angWBR = atan2(WBR(:,1),WBR(:,2));
    angrangeWBR = mean(angWBR) + 3 * std(angWBR) * [-1 1];

    %WS = Warp stack shift
    meanWS = mean(WS);
    stdWS = std(WS);
    maxWS = max(WS);
    minWS = min(WS);
    maxstdWS = meanWS + 3 * stdWS;
    minstdWS = meanWS - 3 * stdWS;
    modWS = sqrt(WS(:,1).^2+WS(:,2).^2);
    rangemodWS = mean(modWS) + 3 * std(modWS) * [-1 1];
    angWS = atan2(WS(:,1),WS(:,2));
    angrangeWS = mean(angWS) + 3 * std(angWS) * [-1 1];

    %Table variables
    Correction = {'green-red';'blue-red';'slice'};
    ModulusMean = [mean(modWGR);mean(modWBR);mean(modWS)];
    StdModulus = [std(modWGR);std(modWBR);std(modWS)];
    Std3Modulus = StdModulus*3;
    AngleMean = [mean(angWGR);mean(angWBR);mean(angWS)];
    StdAngle = [std(angWGR);std(angWBR);std(angWS)];
    Std3Angle = StdAngle * 3;

    T = table(ModulusMean,StdModulus,Std3Modulus,AngleMean,StdAngle,Std3Angle,'RowNames',Correction)

    %% scatterplot of color shift

    h1 = figure('PaperUnits','centimeters','Color',[1 1 1]);
    sng_figcm(fsx,fsy);
    hold on
    line([-10,4],[0,0],'Color',[0 0 0],'LineWidth',0.5)
    line([0,0],[-4,10],'Color',[0 0 0],'LineWidth',0.5)
    s1 = scatter(WBR(:,1),WBR(:,2),'MarkerEdgeColor',Color1,'SizeData',10,'LineWidth',0.5,'Marker','o')
    s2 = scatter(WGR(:,1),WGR(:,2),'MarkerEdgeColor',Color2,'SizeData',10,'LineWidth',0.5,'Marker','o')

    xlabel('x-displacement [pix]','FontSize',8,'FontName','arial');
    ylabel('y-displacement [pix]','FontSize',8,'FontName','arial');
    set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
    set(gca,'XLim',[-10,4])
    set(gca,'YLim',[-4,10])
    set(gca,'Units','centimeters','Position',[1.2 1.2 fsx-1.7 fsy-1.7])

    %[l] = legend([s1,s2],'green to red','blue to red','Location','northwest')
    %this legend scales well with the ScaleFigure function
    [l] = legendflex([s1,s2],{'green to red','blue to red'},'ref',gca,...
        'anchor',{'nw','nw'},'buffer', [10 -10],'padding',[-5 -5 10])

    %sng_figcm(fsx,fsy,113.6)
    %set(0, 'currentfigure', h1)

    export_fig(h1 ,['/Users/samuelgeurts/Desktop/','scat',num2str(1)], '-png', '-r600', '-nocrop');

end


%% Stack correction image mean shift processed on bpfilteredrgb preproc

%fit with a line going through zero
lin2zero = fittype({'x'})
[fit1,gof,fitinfo] = fit(WS(:,1),WS(:,2),lin2zero)


std(fitinfo.residuals)

%fsx = 20,fsy = 20

h1 = figure('PaperUnits','centimeters','Color',[1 1 1]);
sng_figcm(fsx,fsy);

hold on
line([-50,20],[0,0],'Color',[0 0 0],'LineWidth',0.5)
line([0,0],[-20,50],'Color',[0 0 0],'LineWidth',0.5)
s1 = scatter(WS(:,1),WS(:,2),'MarkerEdgeColor',Color3,'SizeData',10,'LineWidth',0.5,'Marker','o')

f1 = plot(fit1)

xlabel('x-displacement [pix]','FontSize',8,'FontName','arial');
ylabel('y-displacement [pix]','FontSize',8,'FontName','arial');
set(gca,'FontName','arial','FontSize',8,'XGrid','on','YGrid','on');
set(gca,'XLim',[-50,20])
set(gca,'YLim',[-20,50])
set(gca,'Units','centimeters','Position',[1.2 1.2 fsx-1.7 fsy-1.7])

%[l] = legend([s1,s2],'green to red','blue to red','Location','northwest')
%this legend scales well with the ScaleFigure function
[l] = legendflex([s1,f1],{'stack correction','trend line: y = -0.35x'},'ref',gca,...
    'anchor',{'nw','nw'},'buffer', [10 -10],'padding',[-5 -5 10])

lg1 = get(l,'children')
set(lg1(2),'XData',[8 20])

%sng_figcm(fsx,fsy,113.6)
%set(0, 'currentfigure', h1)

export_fig(h1 ,['/Users/samuelgeurts/Desktop/','scatstack',num2str(1)], '-png', '-r600', '-nocrop');
