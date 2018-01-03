 
name = 'SpotOptListSlice'
%name = 'SpotOptList'

obj = obj.NewPath(10)
load([obj.SavePath, '/', name, '.mat'], 'SpotOptList')

%%
clear zfinputlist
%zfinputlist.ColorToGrayVectorL = {[-0.380803509067153;0.847820881126686;-0.369037181064066]};

%
zfinputlist.SpotDistLimit = 10;
zfinputlist.ColorToGrayVectorL = {[0; 1; 0]};
zfinputlist.ScaleBaseL = {0.5};
zfinputlist.KthresholdL = {0};
zfinputlist.MPlevelsL = {4:6};
zfinputlist.MPthresholdL = {80 60 40 20};
%zfinputlist.MinSpotSizeL = {8.4084}
%zfinputlist.MaxSpotSizeL = {458}
%zfinputlist.MinProbabilityL = {0.066}
%zfinputlist.MinSpotRange = 1
%zfinputlist.MaxSpotRange = 1
%zfinputlist.MinProbabilityRange = 1
%

zfinputlist.MinSpotRange = [0:0.03:0.3]
zfinputlist.MaxSpotRange = [0.5:0.5:1.2]
zfinputlist.MinProbabilityRange = [0:0.02:0.2]

objTemp = obj.SpotOptimization(1:50,zfinputlist)

%%

load([obj.SavePath, '/', name, '.mat'], 'SpotOptList')
figure;plot([SpotOptList.MeanF1score])

    
    %[SpotOptList.MeanF1score]
    
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



    