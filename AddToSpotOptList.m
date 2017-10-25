function SpotOptList = AddToSpotOptList(SpotOptList,obj,linen)
%puts variables important for valiation in a struct SpotOpt
%saves SpotOpt to the next field in SpotOptList or to a specified field number "linen"

sng_zfinputAssign(obj.zfinput,'SpotDetection')
    
    SpotOpt.MPlevels = obj.zfinput(find(strcmp({obj.zfinput.name},'MPlevels'))).value;
    SpotOpt.MPthreshold = obj.zfinput(find(strcmp({obj.zfinput.name},'MPthreshold'))).value;
    SpotOpt.MinSpotSize = obj.zfinput(find(strcmp({obj.zfinput.name},'MinSpotSize'))).value;
    SpotOpt.MaxSpotSize = obj.zfinput(find(strcmp({obj.zfinput.name},'MaxSpotSize'))).value;
    SpotOpt.MinProbability = obj.zfinput(find(strcmp({obj.zfinput.name},'MinProbability'))).value;

    SpotOpt.MeanPrecision = mean([obj.SpotInfo.Precision]);
    SpotOpt.MeanRecall = mean([obj.SpotInfo.Recall]);
    SpotOpt.MeanF1score = mean([obj.SpotInfo.F1score]);
    SpotOpt.MeanAbsDifference = mean(abs([obj.SpotInfo.AbsDifference]));
    SpotOpt.MeanRelDifference = mean(abs([obj.SpotInfo.RelDifference]));

    SpotOpt.StdPrecision = std([obj.SpotInfo.Precision]);
    SpotOpt.StdRecall = std([obj.SpotInfo.Recall]);
    SpotOpt.StdF1score = std([obj.SpotInfo.F1score]);
    SpotOpt.StdAbsDifference = std(abs([obj.SpotInfo.AbsDifference]));
    SpotOpt.StdRelDifference = std(abs([obj.SpotInfo.RelDifference]));

    SpotOpt.Precision = [obj.SpotInfo.Precision];
    SpotOpt.Recall = [obj.SpotInfo.Recall];
    SpotOpt.F1score = [obj.SpotInfo.F1score];
    SpotOpt.AbsDifference = [obj.SpotInfo.AbsDifference];
    SpotOpt.RelDifference = [obj.SpotInfo.RelDifference];

    SpotOpt.date = date;
    

    if isempty(SpotOptList)
        SpotOptList = SpotOpt;
    elseif exist('linen','var') & linen ~= 0 
        SpotOptList(linen) = SpotOpt;               
    else
        SpotOptList(numel(SpotOptList)+1) = SpotOpt;
    end
   
    
end
