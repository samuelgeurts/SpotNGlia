function obj = TscoreVal(obj)
    
            %load([obj.SavePath,'/',obj.InfoName,'.mat'],'SpotParameters')
            load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected')
            
            
            nfishes = numel(obj.StackInfo);
            
            
%             LinkDistance = cell(nfishes, 1);
%             CorrectSpots = cell(nfishes, 1);
%             nCorrect = zeros(nfishes, 1);
%             FalsePosSpots = cell(nfishes, 1);
%             nFalsePos = zeros(nfishes, 1);
%             FalseNegSpots = cell(nfishes, 1);
%             nFalseNeg = zeros(nfishes, 1);
%             Precision = zeros(nfishes, 1);
%             Recall = zeros(nfishes, 1);
%             F1score = zeros(nfishes, 1);
%             AbsDifference = zeros(nfishes, 1);
%             RelDifference = zeros(nfishes, 1);
            
            %SpotCom = cell(nfishes,1);
            
            %value added to prevent for substraction by zero
            a = 0.001;
            
            for k1 = 1:nfishes
                
                Spotpar = SpotsDetected{k1};
                %ambr = obj.BrainInfo(k1).AnnotatedMidBrainRegistrated;
                
                SpotCom = reshape([Spotpar.Centroid], 2, numel(Spotpar))';
                ambs = obj.Annotations(k1).Spots;
                
                
                [Correct, FalsePos, FalseNeg, link] = sng_CoordinateMatching ...
                    (SpotCom, ambs, 10);
                
                LinkDistance{k1} = link;
                
                if exist('Correct', 'var')
                    CorrectSpots{k1} = Correct;
                    nCorrect(k1) = size(Correct, 1);
                else
                    CorrectSpots{k1} = [];
                    nCorrect(k1) = 0;
                end
                if exist('FalsePos', 'var')
                    FalsePosSpots{k1} = FalsePos;
                    nFalsePos(k1) = size(FalsePos, 1);
                else
                    FalsePosSpots{k1} = [];
                    nFalsePos(k1) = 0;
                end
                if exist('FalseNeg', 'var')
                    FalseNegSpots{k1} = FalseNeg;
                    nFalseNeg(k1) = size(FalseNeg, 1);
                else
                    FalseNegSpots{k1} = [];
                    nFalseNeg(k1) = 0;
                end
                
                Precision(k1) = nCorrect(k1) / (nCorrect(k1) + nFalsePos(k1) + a);
                Recall(k1) = nCorrect(k1) / (nCorrect(k1) + nFalseNeg(k1) + a);
                F1score(k1) = 2 * (Precision(k1) * Recall(k1)) / (Precision(k1) + Recall(k1) + a);
                AbsDifference(k1) = nFalsePos(k1) - nFalseNeg(k1);
                RelDifference(k1) = (nFalsePos(k1) - nFalseNeg(k1)) / nCorrect(k1);
                
                obj.SpotBrainStats.FiveNumberSummaryPrecision = sng_FiveNumberSum(Precision);
                obj.SpotBrainStats.MeanPrecision = mean(Precision);
                
                obj.SpotBrainStats.FiveNumberSummaryRecall = sng_FiveNumberSum(Recall);
                obj.SpotBrainStats.MeanRecall = mean(Recall);
                
                obj.SpotBrainStats.FiveNumberSummaryF1score = sng_FiveNumberSum(F1score);
                obj.SpotBrainStats.MeanF1score = mean(F1score);
                
                obj.SpotBrainStats.FiveNumberSummaryAbsDifference = sng_FiveNumberSum(AbsDifference);
                obj.SpotBrainStats.MeanAbsDifference = mean(AbsDifference);
                
                obj.SpotBrainStats.FiveNumberSummaryRelDifference = sng_FiveNumberSum(RelDifference);
                obj.SpotBrainStats.MeanRelDifference = mean(RelDifference);
                
            end
            
            %store variables in obj.SpotInfo
            [obj.SpotBrainInfo(1:nfishes).LinkDistance] = LinkDistance{:};
            [obj.SpotBrainInfo(1:nfishes).CorrectSpots] = CorrectSpots{:};
            temp = num2cell(nCorrect); [obj.SpotBrainInfo(1:nfishes).nCorrect] = temp{:};
            [obj.SpotBrainInfo(1:nfishes).FalsePosSpots] = FalsePosSpots{:};
            temp = num2cell(nFalsePos); [obj.SpotBrainInfo(1:nfishes).nFalsePos] = temp{:};
            [obj.SpotBrainInfo(1:nfishes).FalseNegSpots] = FalseNegSpots{:};
            temp = num2cell(nFalseNeg); [obj.SpotBrainInfo(1:nfishes).nFalseNeg] = temp{:};
            temp = num2cell(Precision); [obj.SpotBrainInfo(1:nfishes).Precision] = temp{:};
            temp = num2cell(Recall); [obj.SpotBrainInfo(1:nfishes).Recall] = temp{:};
            temp = num2cell(F1score); [obj.SpotBrainInfo(1:nfishes).F1score] = temp{:};
            temp = num2cell(AbsDifference); [obj.SpotBrainInfo(1:nfishes).AbsDifference] = temp{:};
            temp = num2cell(RelDifference); [obj.SpotBrainInfo(1:nfishes).RelDifference] = temp{:};
        end