clear all

sng_MultiDir(who,'PreprocessionPath')

drawnow

load([PreprocessionPath,'/imageinfo.mat'],'imageinfo')
load([PreprocessionPath,'/stackinfo.mat'],'stackinfo')
load([PreprocessionPath,'/input.mat'],'zfinput')


%% measure the warp
clearvars -except imageinfo stackinfo FolderPath

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

%WGR = Warp Green-Red Channelw
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

%scatterplot of shift
figure;scatter(WBR(:,1),WBR(:,2))
hold on;scatter(WGR(:,1),WGR(:,2))
line([-10,4],[0,0])
line([0,0],[-4,10])
legend






