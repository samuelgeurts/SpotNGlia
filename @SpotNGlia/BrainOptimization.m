    function BrainOptimization(obj, fishnumbers)
            
            %CompleteTemplate = LoadTemplateLink3([obj.SourcePath, '/', 'Template 3 dpf']);
            
            
            if exist([obj.SavePath, '/', 'BrainOptList', '.mat'], 'file')
                load([obj.SavePath, '/', 'BrainOptList', '.mat'], 'BrainOptList')
            else
                BrainOptList = [];
                save([obj.SavePath, '/', 'BrainOptList', '.mat'], 'BrainOptList')
            end
            
            obj = LoadTemplate(obj);
                        
            if ~exist('fishnumbers', 'var')
                fishnumbers = 1:numel(obj.StackInfo);
            elseif max(fishnumbers) > numel(obj.StackInfo)
                error('at least one fish does not exist in RegistrationInfo')
            end
            
            nfishes = numel(fishnumbers);
            
            %{
            ColorToGrayVectorL = {[0; 1; 0]}
            ScaleBaseL = {0.5}
            KthresholdL = {0}
            MPlevelsL = {5:7};
            MPthresholdL = {256};
            %MinSpotSizeL = {}
            %MaxSpotSizeL = {}
            %MinProbabilityL = {}
            
            %stores every combination of indices given by the size of the
            %variables above
            [a, b, c, d, e] = ndgrid(1:numel(ColorToGrayVectorL), ...
                1:numel(ScaleBaseL), ...
                1:numel(KthresholdL), ...
                1:numel(MPlevelsL), ...
                1:numel(MPthresholdL));
            
            %create different variable-sets to test select only spotinfo parameters
            for k2 = 1:numel(a)
                ZFParametersTemp{k2} = obj.ZFParameters(strcmp({obj.ZFParameters.stage}, 'SpotDetection'));
                ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'RgbToGray', 'ColorToGrayVector', ColorToGrayVectorL{a(k2)}, ''); %color selection
                ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'Wavelet', 'ScaleBase', ScaleBaseL{b(k2)}, ''); %color selection
                ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'Wavelet', 'Kthreshold', KthresholdL{c(k2)}, ''); %color selection
                ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'MultiProduct', 'MPlevels', MPlevelsL{d(k2)}, ''); %color selection
                ZFParametersTemp{k2} = sng_zfinput(ZFParametersTemp{k2}, 0, 'SpotDetection', 'MultiProduct', 'MPthreshold', MPthresholdL{e(k2)}, ''); %color selection
                %ZFParametersTemp = sng_zfinput(ZFParametersTemp,0,'SpotDetection','SpotSelection','MinSpotSize',MinSpotSizeL{g},''); %color selection
                %ZFParametersTemp = sng_zfinput(ZFParametersTemp,0,'SpotDetection','SpotSelection','MaxSpotSize',MaxSpotSizeL{h},''); %color selection
                %ZFParametersTemp = sng_zfinput(ZFParametersTemp,0,'SpotDetection','SpotSelection','MinProbability',MinProbabilityL{i},''); %color selection
            end
            
            %compute Spots
            %}
            
            for k2 = 1%:numel(ZFParametersTemp)
                %compute spot parameters for certain ZFParemeter input
                for k1 = 1:nfishes
                    fn = fishnumbers(k1);
                    %fprintf([num2str(a(k2)), ',', num2str(b(k2)), ',', num2str(c(k2)), ',', num2str(d(k2)), ',', num2str(e(k2)), ',', num2str(k1), '\n'])
                    ambr = [obj.Annotations(fn).MidBrain];
                    AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
                    [~, BrainSegmentationInfo(k1)] = MidBrainDetectionSNG(AlignedFish, obj.CompleteTemplate, obj.ZFParameters);
                k1
                end
                

                objTemp = obj
                objTemp.BrainInfo = []
                objTemp.BrainStats = []
                objTemp.SpotInfo = []
                objTemp.SpotStats = []               
                objTemp.SpotBrainInfo = []
                objTemp.SpotBrainStats = []             
                
                
                
                objTemp = objTemp.BrainVal(BrainSegmentationInfo,1:numel(BrainSegmentationInfo));
                
                objTemp.ShowBoxPlot
                %{
                disp(num2str(mean([obj.SpotInfo.F1score])));
                
                %disp([num2str(l2),' ',num2str(l1),' ',num2str(mean([obj.SpotInfo.F1score]))]);
                
                SpotOpt.ColorToGrayVectorL = ColorToGrayVectorL{a(k2)};
                SpotOpt.ScaleBaseL = ScaleBaseL{b(k2)};
                SpotOpt.KthresholdL = KthresholdL{c(k2)};
                SpotOpt.MPlevelsL = MPlevelsL{d(k2)};
                SpotOpt.MPthresholdL = MPthresholdL{e(k2)};
                SpotOpt.MinSpotSize = MinSpotSize;
                SpotOpt.MaxSpotSize = MaxSpotSize;
                SpotOpt.MinProbability = MinProbability;
                
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
                %}
                BrainOpt.date = date;
                
                if isempty(BrainOptList)
                    BrainOptList = SpotOpt;
                    %elseif exist('linen','var') & linen ~= 0
                    %    BrainOptList(linen) = SpotOpt;
                else
                    BrainOptList(numel(BrainOptList)+1) = SpotOpt;
                end
            end
            save([obj.SavePath, '/', 'BrainOptList', '.mat'], 'BrainOptList');
            
        end
  