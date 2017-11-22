load([obj.SavePath, '/', 'SpotOptList', '.mat'], 'SpotOptList')
 
obj = obj.NewPath(10)

%%
CompleteTemplate = LoadTemplateSNG([obj.SourcePath, '/', 'Template 3 dpf']);

clear zfinputlist
%zfinputlist.ColorToGrayVectorL = {[-0.380803509067153;0.847820881126686;-0.369037181064066]};

zfinputlist.ColorToGrayVectorL = {[0; 1; 0]};
zfinputlist.ScaleBaseL = {0.6};
zfinputlist.KthresholdL = {0};
%zfinputlist.MPlevelsL = {4:6;5:7;6:8;4:7;5:8};
zfinputlist.MPlevelsL = {5:7};

%zfinputlist.MPthresholdL = {0.5 1 2 4 8 16 32 64};
zfinputlist.MPthresholdL = {255 280 330 355};

%{
zfinputlist.MinSpotSizeL = {8.4084}
zfinputlist.MaxSpotSizeL = {458}
zfinputlist.MinProbabilityL = {0.0660}
zfinputlist.MinSpotRange = 1
zfinputlist.MaxSpotRange = 1
zfinputlist.MinProbabilityRange = 1
%}
zfinputlist.MinSpotRange = [0:0.05:0.2]
zfinputlist.MaxSpotRange = [0.4:0.05:1.2]
zfinputlist.MinProbabilityRange = [0:0.05:0.4]



        
obj.SpotOptimization(1:40,zfinputlist)

%%

    load([obj.SavePath, '/', 'SpotOptList', '.mat'], 'SpotOptList')


    figure;plot([SpotOptList.MeanF1score])

    
    [SpotOptList.MeanF1score]
    
       
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
    
figure;plot([SpotOptList3.MeanF1score])

%%

save([obj.SavePath, '/', 'SpotOptList', '.mat'], 'SpotOptList')

%%



    