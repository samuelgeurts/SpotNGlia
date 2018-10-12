classdef SNGInputParameters 
%class containting SpotNGlia input parameters
    properties
        ImageInfo
        Preprocession
        ExtendedDeptofField
        Registration
        BrainSegmentation
        SpotDetection
        
        StoredInputParameters
    end
    
    
    methods
        
        function objI = SNGInputParameters
            
            if isempty(objI.ImageInfo)
                objI = objI.defaultInput;
            end
        end
            
        function objI = defaultInput(objI)
            
            %default since 2017
            objI.ImageInfo.levels = 1;
            objI.ImageInfo.iterations = 10;
            objI.ImageInfo.scale = 0.0625;
            objI.ImageInfo.threshold = 0.9600;
            objI.ImageInfo.threshold2 = 2;
            objI.ImageInfo.sorting = [];
            
            objI.Preprocession.onoff = 0;
            objI.Preprocession.sigmalp = 1;
            objI.Preprocession.sigmahp = 4;
            objI.Preprocession.scaleC = 0.2500;
            objI.Preprocession.levelsC = 2;
            objI.Preprocession.iterationsC = 10;
            objI.Preprocession.scaleS = 0.25;
            objI.Preprocession.levelsS = 3;
            objI.Preprocession.iterationsS = 20;
            
            objI.ExtendedDeptofField.variancedisksize = 7;
            
            objI.Registration.Method = 'TriangleSmooth';
            objI.Registration.Smooth = 1;
            objI.Registration.ChannelMethod = 'cuboid';
            objI.Registration.AngleSteps1 = 100;
            objI.Registration.Scale1 = 0.0625;
            objI.Registration.AngleRange2 = 0.0100;
            objI.Registration.AngleSteps2 = 50;
            objI.Registration.FinalCenteredYCrop = 1000;
            objI.Registration.ScaleRange = [0.6, 1.3];
            objI.Registration.ScaleSteps = 200;
            objI.Registration.Scale2 = 0.0625;
            objI.Registration.RemoveTail = 600;
            objI.Registration.Scale3 = 0.25;
            objI.Registration.AffineMethod = 'translation';
            objI.Registration.Levels = 3;
            objI.Registration.Iterations = 40;
            
            objI.BrainSegmentation.Method = 'TriangleSmooth';
            objI.BrainSegmentation.Method2 = 4;
            
            objI.SpotDetection.ColorToGrayVector = [0;1;0];
            objI.SpotDetection.ScaleBase = 0.5;
            objI.SpotDetection.ScaleLevels = 7;
            objI.SpotDetection.Kthreshold = 0;
            objI.SpotDetection.MPlevels = [5,6,7];
            objI.SpotDetection.MPthreshold = 200;
            objI.SpotDetection.MinSpotSize = 8.4;
            objI.SpotDetection.MaxSpotSize = 458;
            objI.SpotDetection.MinProbability = 0.066;
            objI.SpotDetection.SpotDistLimit = 10;
            objI.SpotDetection.ComputeOnSlice = 0;
             
        end
            
        function assign(objI,fieldString)
            %This function assigns fields of the given fieldstring to the 
            %workspace. Tthis function maybe can be removed later on as direct
            %calling the propery is also possible.
            fieldNames = fields(objI.(fieldString));
            for iField = 1:numel(fieldNames)
                assignin('caller',fieldNames{iField},objI.(fieldString).(fieldNames{iField}));
            end
            
        end
        
    end
end

            

                






