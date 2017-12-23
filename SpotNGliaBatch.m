classdef SpotNGliaBatch
    
    properties
        FilePaths = [];
        ObjNames = [];
        Objects = [];
        
        MainFolder = '/Users/samuelgeurts/Desktop/PB optimalisatie'
        MainFolderOriginal = '/Volumes/Seagate Expansion Drive/PublicationBatch'
        SaveName = []
        
        BatchStats = [];
        SpotsAnnotated = [];
        SpotsComputed = [];
        
    end
    
    methods
        
        
        function objb = SpotNGliaBatch
            
            %{
                obj = obj.NewPath(1)
                filename = obj.SavePath
                infoname = obj.InfoName
                sngname = obj.SaveName
            %}
            
            %% initialize locations and file names
            
            objb.FilePaths{1} = '2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/Csf1a B1';
            objb.ObjNames{1} = 'SNG_Csf1a B1';
            
            objb.FilePaths{2} = '2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/Csf1b B1';
            objb.ObjNames{2} = 'SNG_Csf1b B1';
            
            objb.FilePaths{3} = '2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/Il 34 B1';
            objb.ObjNames{3} = 'SNG_Il 34 B1';
            
            objb.FilePaths{4} = '2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/WT B1';
            objb.ObjNames{4} = 'SNG_WT B1';
            
            objb.FilePaths{5} = '2016-12-01-Slco2b1-3dpf';
            objb.ObjNames{5} = 'SNG_2016-12-01-Slco2b1-3dpf';
            
            objb.FilePaths{6} = '2017-01-19-f13a1a.1-loc56';
            objb.ObjNames{6} = 'SNG_2017-01-19-f13a1a';
            
            objb.FilePaths{7} = '2017-01-26-fabp11a-Havcr2-3dpf';
            objb.ObjNames{7} = 'SNG_2017-01-26-fabp11a-Havcr2-3dpf';
            
            objb.FilePaths{8} = '2017-02-03-p2rx3a-hcst-3dpf';
            objb.ObjNames{8} = 'SNG_2017-02-03-p2rx3a-hcst-3dpf';
            
            objb.FilePaths{9} = '2016-04-21 ch25h/Ch25h B2';
            objb.ObjNames{9} = 'SNG_Ch25h B2';
            
            objb.FilePaths{10} = '2016-04-21 ch25h/WT B2';
            objb.ObjNames{10} = 'SNG_WT B2';
            
            %objb.FilePaths{11} = '20171026_DKO_NR_3dpf'
            %objb.ObjNames{11} = 'SNG_20171026_DKO_NR_3dpf'
            
        end
        
        function objb = LoadObjects(objb, Original)
            %stores SpotNGlia objects under the struct Object with increasing numbers
            %if Original = true, the files from the harddisk are taken (objb.MainFolderOriginal)
            for k1 = 1:numel(objb.FilePaths)
                if exist('Original)', 'var') && Original
                    object = load([objb.MainFolderOriginal, '/', objb.FilePaths{k1}, '/', objb.ObjNames{k1}, '.mat']);
                else
                    object = load([objb.MainFolder, '/', objb.FilePaths{k1}, '/', objb.ObjNames{k1}, '.mat']);
                end
                %objb.Objects.(['obj', num2str(k1)]) = object.obj;
                objb.Objects{k1} = object.obj;
            end
        end
        
        function UpdateSavePath(objb, Original)
            for k1 = 1:numel(objb.FilePaths)
                %obj = objb.Objects.(['obj', num2str(k1)]);
                obj = objb.Objects{k1};
                if exist('Original)', 'var') && Original
                    obj.SavePath = [objb.MainFolderOriginal, '/', objb.FilePaths{k1}];
                else
                    obj.SavePath = [objb.MainFolder, '/', objb.FilePaths{k1}];
                end
                obj.saveit;
            end
        end
        
        function Excel2StackInfo(objb)
            %% get human counted dat from excel file ans save to obj under obj.StackInfo
            % ugly function for a very specifik file.
            
            excelfile = '/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Multiple Batch Comparison/Multiple Batch Comparison.xlsx';
            
            RangeS = {'B4:B14', 'E4:E16', 'H4:H15', 'K4:K17'};
            RangeD = {2:12, 2:14, 2:13, 2:15};
            
            for k1 = 1:10
                obj = objb.Objects{k1};
                
                
                if k1 <= 4
                    num = num2cell(xlsread(excelfile, 2, RangeS{k1}));
                    [obj.Annotations(RangeD{k1}).Counts] = num{:};
                end
                
                if k1 == 5
                    num1 = num2cell(xlsread(excelfile, 2, 'N3:N13'));
                    num2 = num2cell(xlsread(excelfile, 2, 'Q3:Q16'));
                    [obj.Annotations(1:11).Counts] = num1{:};
                    [obj.Annotations(12:25).Counts] = num2{:};
                end
                
                if k1 == 6
                    num1 = num2cell(xlsread(excelfile, 2, 'T3:T11'));
                    num2 = num2cell(xlsread(excelfile, 2, 'W3:W15'));
                    num3 = num2cell(xlsread(excelfile, 2, 'Z3:Z10'));
                    num4 = num2cell(xlsread(excelfile, 2, 'AC3:AC14'));
                    num5 = num2cell(xlsread(excelfile, 2, 'AF3:AF14'));
                    
                    [obj.Annotations(1:9).Counts] = num1{:};
                    [obj.Annotations(10:21).Counts] = num2{:};
                    [obj.Annotations(23:30).Counts] = num3{:};
                    [obj.Annotations(31:42).Counts] = num4{:};
                    [obj.Annotations(43:53).Counts] = num5{:};
                    [obj.Annotations(54).Counts] = NaN;
                end
                
                if k1 == 7
                    
                    num1 = num2cell(xlsread(excelfile, 2, 'AI3:AI18'));
                    num2 = num2cell(xlsread(excelfile, 2, 'AL3:AL13'));
                    num3 = num2cell(xlsread(excelfile, 2, 'AO3:AO18'));
                    
                    [obj.Annotations(1:16).Counts] = num1{:};
                    [obj.Annotations(17:27).Counts] = num2{:};
                    [obj.Annotations(28:43).Counts] = num3{:};
                end
                
                
                if k1 == 8
                    
                    num1 = num2cell(xlsread(excelfile, 2, 'AR3:AR17'));
                    num2 = num2cell(xlsread(excelfile, 2, 'AU3:AU13'));
                    num3 = num2cell(xlsread(excelfile, 2, 'AX3:AX17'));
                    num4 = num2cell(xlsread(excelfile, 2, 'BA3:BA9'));
                    %num5 = num2cell(xlsread(excelfile,2,'BD3:BD10'));
                    %num6 = num2cell(xlsread(excelfile,2,'BF3:BF5'));
                    
                    [obj.Annotations(1:15).Counts] = num1{:};
                    [obj.Annotations(16:26).Counts] = num2{:};
                    [obj.Annotations(27:41).Counts] = num3{:};
                    [obj.Annotations(42:48).Counts] = num4{:};
                    
                end
                
                if k1 == 9
                    num1 = num2cell(xlsread(excelfile, 2, 'BJ3:BJ13'));
                    [obj.Annotations(1:11).Counts] = num1{:};
                end
                
                if k1 == 10
                    num1 = num2cell(xlsread(excelfile, 2, 'BM3:BM15'));
                    [obj.Annotations(1:13).Counts] = num1{:};
                end
                
                if k1 == 11
                    disp('this batch had unknown annotated data');
                end
                
                for k2 = 1:numel(obj.Annotations)
                    if isempty(obj.Annotations(k2).Counts)
                        obj.Annotations(k2).Counts = NaN;
                    end
                end
                
                obj.saveit
            end
        end
        
        function Computation(objb)
            fprintf('%d ', k1)
            for k1 = 1:10
                objb.Objects{k1}.ZFParameters = sng_zfinput(objb.Objects{k1}.ZFParameters, 0, 'SpotDetection', 'RgbToGray', 'ColorToGrayVector', [0; 1; 0], ''); %color selection
                objb.Objects{k1}.ZFParameters = sng_zfinput(objb.Objects{k1}.ZFParameters, 0, 'SpotDetection', 'Wavelet', 'ScaleBase', 0.333, ''); %color selection
                objb.Objects{k1}.ZFParameters = sng_zfinput(objb.Objects{k1}.ZFParameters, 0, 'SpotDetection', 'MultiProduct', 'MPlevels', 6:11, ''); %color selection
                objb.Objects{k1}.ZFParameters = sng_zfinput(objb.Objects{k1}.ZFParameters, 0, 'SpotDetection', 'MultiProduct', 'MPthreshold', 1500, ''); %color selection
                objb.Objects{k1}.ZFParameters = sng_zfinput(objb.Objects{k1}.ZFParameters, 0, 'SpotDetection', 'SpotSelection', 'MinSpotSize', 27.7, ''); %color selection
                objb.Objects{k1}.ZFParameters = sng_zfinput(objb.Objects{k1}.ZFParameters, 0, 'SpotDetection', 'SpotSelection', 'MaxSpotSize', 417, ''); %color selection
                objb.Objects{k1}.ZFParameters = sng_zfinput(objb.Objects{k1}.ZFParameters, 0, 'SpotDetection', 'SpotSelection', 'MinProbability', 0.0719, ''); %color selection
                
                %zf = load([ob{k1}.obj.SourcePath, '/', 'zfinput.mat']);
                %obj.ZFParameters = zf.zfinput;
            end
            
            cellobject = objb.Objects;
            parfor k1 = 1:10
                k1
                cellobject{k1} = cellobject{k1}.SpotDetection;
            end
            objb.Objects = cellobject;
        end
        %{
function NewPath(objb, mode)
              %% change savepath and fishpath to new path only for
              for k1 = 1:numel(objb.FilePaths)
                  obj = objb.Objects.(['obj', num2str(k1)])
                  tempvar = obj.SavePath
                  obj = obj.NewPath(mode);
                  if ~isequal(obj.SavePath, tempvar)
                      obj.saveit;
                  end
              end
              %%
          end
        %}
        
        function objb = computeTscores(objb)
            for k1 = 1:10
                fprintf('%d ', k1);
                
                obj = objb.Objects{k1};
                obj = obj.TtestVal;
                obj.saveit
                
                objb.BatchStats(k1).ttest = obj.SpotBrainStats.ttest;
                objb.BatchStats(k1).ttestval = obj.SpotBrainStats.ttestval;
                objb.BatchStats(k1).MeanAbsDif = obj.SpotBrainStats.MeanAbsDifference;
                objb.BatchStats(k1).StdAbsDif = obj.SpotBrainStats.StdAbsDifference;
                objb.BatchStats(k1).RMSD = obj.SpotBrainStats.RMSD;
                
                %load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected');
                %ccs = zeros(1, numel(SpotsDetected));
                %ccs = obj.StackInfo.Counts
                %for l1 = 1:numel(SpotsDetected)
                %    ccs(l1) = numel(SpotsDetected{l1});
                %end
                %temp = {ob{k1}.obj.Annotations.Counts};
                objb.SpotsAnnotated{k1} = [obj.Annotations.Counts];
                objb.SpotsComputed{k1} = [obj.StackInfo.Counts];
            end
        end
        
        function objb = HistPar(objb)
            for k1 = 1:10
                %is implemented in obj.Registration
                %fprintf('%d ',k1)
                obj = objb.Objects{k1};
                %obj = obj.HistPar;
                %obj.saveit;
                objb.BatchStats(k1).MeanFishColor = obj.BatchInfo.MeanFishColor;
                objb.BatchStats(k1).StdMeanFishColor = obj.BatchInfo.StdMeanFishColor;
                objb.BatchStats(k1).MeanMaxFishColor = obj.BatchInfo.MeanMaxFishColor;
                objb.BatchStats(k1).StdMaxFishColor = obj.BatchInfo.StdMaxFishColor;
            end
        end
        
        function saveit(objb, name)
            objb.SaveName = name;
            save(strcat(objb.MainFolder, '/', name, '.mat'), 'objb');
        end
                
        function ShowFishHeadHist(objb,batchnumber,fishnumber)
                obj = objb.Objects{batchnumber};                
                ShowFishHeadHist(obj, fishnumber);
        end
        
        function ShowMaxFishHist(objb)
            for k1 = 1:10
                obj = objb.Objects{k1};              
                ShowMaxFishHist(obj);
            end
        end
        
    end
end