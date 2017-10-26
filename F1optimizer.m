 

function Maxf1
 

            MinSpotSizeT = MinSpotSize * (l3/10);
            MinProbabilityT = MinProbability * (l4/10);
                
            %%% assign new threshold logical array to obj.SpotParameters
            for k5 = stacknumbers
                temp = num2cell([obj.SpotParameters{k5}.Area] >= MinSpotSizeT);
                [obj.SpotParameters{k5}.LargerThan] = temp{:};

                temp = num2cell([obj.SpotParameters{k5}.Area] <= MaxSpotSize);
                [obj.SpotParameters{k5}.SmallerThan] = temp{:};

                temp = num2cell([obj.SpotParameters{k5}.ColorProbability] >= MinProbabilityT);
                [obj.SpotParameters{k5}.MinProbability] = temp{:};
            end

            obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinSpotSize',MinSpotSizeT,'');      % size selection 
            obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MaxSpotSize',MaxSpotSize,'');      % size selection 
            obj.zfinput = sng_zfinput(obj.zfinput,0,'SpotDetection','SpotSelection','MinProbability',MinProbabilityT,'');% color selection
            
            obj = obj.SpotVal;
        
                SpotOpt.MeanF1score = mean([obj.SpotInfo.F1score]);

            
            %load([obj.SavePath,'/','SpotOptList5','.mat'],'SpotOptList5');  
            SpotOptList5 = AddToSpotOptList(SpotOptList5,obj); %stores all relevant info into the next field of SpotOptList
            %save([obj.SavePath,'/','SpotOptList5','.mat'],'SpotOptList5');    

end

function f = objectivefcn1(x)
f = 0;
for k = -10:10
    f = f + exp(-(x(1)-x(2))^2 - 2*x(1)^2)*cos(x(2))*sin(2*x(2));
end

Start at x0 = [0.25,-0.25] and search for a minimum of objectivefcn.

x0 = [0.25,-0.25];
x = fminsearch(@objectivefcn1,x0)