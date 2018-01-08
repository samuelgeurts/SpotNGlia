
%minspotsize:   0-20-60
%maxspotsize    200-300-470
%MinProbability 0.03-0.06-0.12

name = 'SpotOptList2'

obj = obj.NewPath(10)
load([obj.SavePath, '/', name, '.mat'], 'SpotOptList2')

%
clear inp
%inp.ColorToGrayVectorL = {[-0.380803509067153;0.847820881126686;-0.369037181064066]};
%objb.Objects{2}.BrainSegmentation
%objb.Objects{2}.CheckBrain

%{
    inp.ColorToGrayVectorL = {[0; 1; 0]};
    inp.ScaleBaseL = {0.5};
    inp.KthresholdL = {0};
    inp.MPlevelsL = {5:7};
    inp.MPthresholdL = {200};
    inp.MinSpotSizeL = 8.4;
    inp.MaxSpotSizeL = 458;
    inp.MinProbabilityL = 0.066;
    inp.ComputeOnSlice = 1
    inp.SpotDistLimit = 10;
    inp.brain_sw = 'Annotation' 
    inp.spot_sw = 'Annotation'
%}
inp.ColorToGrayVectorL = {[0; 1; 0]};
inp.ScaleBaseL = {0.5};
inp.KthresholdL = {0};
inp.MPlevelsL = {5:7};
inp.MPthresholdL = {300};
inp.MinSpotSizeL = [5 10 15 20 25 30 45 50];
inp.MaxSpotSizeL = [350 400 450 500];
inp.MinProbabilityL = [8:1:20]/100;
inp.ComputeOnSlice = 1
inp.SpotDistLimit = 10;
inp.brain_sw = 'Annotation' %
inp.spot_sw = 'Annotation', %'Counts'

obj = obj.SpotOptimization2([],inp,[],'SpotOptList2');

inp.MPlevelsL = {4:6,5:7};
inp.MPthresholdL = {100 150 200 250 300 350};

obj = obj.SpotOptimization2([],inp,[],'SpotOptList2');






load([obj.SavePath, '/', name, '.mat'], 'SpotOptList2')
%{

parfor k1 = 1:10
k1

obj = objb.Objects{k1}.SpotOptimization2([],inp,'Annotation');

%obj2 = objb.Objects{k1}.SpotOptimizationRMSD([],inp,'Annotation');

%%
%load([obj.SavePath, '/', 'SpotOptListRMSD', '.mat'], 'SpotOptListRMSD')
%open SpotOptListRMSD
end

%}




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



    