
%minspotsize:   0-20-60
%maxspotsize    200-300-470
%MinProbability 0.03-0.06-0.12



obj = obj.NewPath(10)
load([obj.SavePath, '/', 'SpotOptList', '.mat'], 'SpotOptList')

%%
clear zfinputlist
%zfinputlist.ColorToGrayVectorL = {[-0.380803509067153;0.847820881126686;-0.369037181064066]};

objb.Objects{2}.BrainSegmentation
objb.Objects{2}.CheckBrain




zfinputlist.ColorToGrayVectorL = {[0; 1; 0]};
zfinputlist.ScaleBaseL = {0.5};
zfinputlist.KthresholdL = {0};
zfinputlist.MPlevelsL = {5:7};
zfinputlist.MPthresholdL = {200};
zfinputlist.MinSpotSizeL = [20];
zfinputlist.MaxSpotSizeL = [200];
zfinputlist.MinProbabilityL = [0:0.25:7]/100;



parfor k1 = 1:10
k1
obj = objb.Objects{k1}.SpotOptimizationRMSD([],zfinputlist,'Correction');

%%
%load([obj.SavePath, '/', 'SpotOptListRMSD', '.mat'], 'SpotOptListRMSD')
%open SpotOptListRMSD
end
%figure;plot([SpotOptList.MeanF1score])
    
    
%{       
%% remove duplicates    
[~,ia] = unique([[SpotOptList.MeanF1score]',[SpotOptList.MeanPrecision]',[SpotOptList.MeanRecall]'],'rows')    
SpotOptList2 = SpotOptList(sort(ia))
   
    
%% sort spotoptlist based on Color2GrayVector -> MPlevels -> MPthreshold
[~,ia,types] = unique([SpotOptList2.ColorToGrayVectorL]','rows')

SpotOptList3 = struct

for k1 = 1:numel(ia)
    SpotOptTemp = SpotOptList2(types == k1)
    
    temp = cell(0)
    for k2 = 1:numel(SpotOptTemp)
        temp{k2} = num2str(SpotOptTemp(k2).MPlevelsL)
    end
    
    [~,ib,types2] = unique(temp)
    for k3 = 1:numel(ib)
       SpotOptTemp2 = SpotOptTemp(types2 == k3)
       

         [~,indices] = sort([SpotOptTemp2.MPthresholdL])

         if isempty(fields(SpotOptList3))
            SpotOptList3 = SpotOptTemp2(indices)
         else
             n = numel(SpotOptList3);
             n2 = numel(SpotOptTemp2);
             
             SpotOptList3(n+1:n+n2) = SpotOptTemp2(indices)
         end
    end
end
    %}
%figure;plot([SpotOptList3.MeanF1score])

%%

%save([obj.SavePath, '/', 'SpotOptList', '.mat'], 'SpotOptList')

%%



    