classdef SpotNGliaBatch
    
    properties
        FilePaths = cell(1);
        ObjNames = cell{1};
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
            FilePaths{1} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/Csf1a B1'
            ObjNames{1} = 'SNG_Csf1a B1'
            
            FilePaths{2} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/Csf1b B1'
            ObjNames{2} = 'SNG_Csf1b B1'
            
            FilePaths{3} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/Il 34 B1'
            ObjNames{3} = 'SNG_Il 34 B1'
            
            FilePaths{4} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-28 IL34 Csf1a Csf1b Bnc2 3 dpf/WT B1'
            ObjNames{4} = 'SNG_WT B1'
            
            FilePaths{5} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-12-01-Slco2b1-3dpf'
            I{5} = 'INFO_2016-12-01-Slco2b1-3dpf'
            ObjNames{5} = 'SNG_2016-12-01-Slco2b1-3dpf'
            
            FilePaths{6} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2017-01-19-f13a1a.1-loc56'
            ObjNames{6} = 'SNG_2017-01-19-f13a1a'
            
            FilePaths{7} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2017-01-26-fabp11a-Havcr2-3dpf'
            ObjNames{7} = 'SNG_2017-01-26-fabp11a-Havcr2-3dpf'
            
            FilePaths{8} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2017-02-03-p2rx3a-hcst-3dpf'
            ObjNames{8} = 'SNG_2017-02-03-p2rx3a-hcst-3dpf'
            
            FilePaths{9} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-21 ch25h/Ch25h B2'
            ObjNames{9} = 'SNG_Ch25h B2'
            
            FilePaths{10} = '/Volumes/Seagate Expansion Drive/PublicationBatch/2016-04-21 ch25h/WT B2'
            ObjNames{10} = 'SNG_WT B2'
            
            FilePaths{11} = '/Volumes/Seagate Expansion Drive/PublicationBatch/20171026_DKO_NR_3dpf'
            ObjNames{11} = 'SNG_20171026_DKO_NR_3dpf'
            
            
            %% change savepath and fishpath to new path
            for k1 = 1:numel(objb.FilePaths)
                load([FilePaths{k1}, '/', ObjNames{k1}, '.mat']);
                tempvar = obj.SavePath
                obj = obj.NewPath(12);
                if ~isequal(obj.SavePath, tempvar)
                    obj.saveit;
                end
            end
            %%
        end
        
        
        function Excel2StackInfo(obj)
            %% get human counted dat from excel file ans save to obj under obj.StackInfo
            
            excelfile = '/Volumes/Macintosh HD/OneDrive/BIGR Zebrafish/Multiple Batch Comparison/Multiple Batch Comparison.xlsx'
            
            RangeS = {'B4:B14', 'E4:E16', 'H4:H15', 'K4:K17'}
            RangeD = {2:12, 2:14, 2:13, 2:15}
            
            for k1 = 1:11
                load([P{k1}, '/', S{k1}, '.mat']);
                
                if k1 <= 4
                    num = num2cell(xlsread(excelfile, 2, RangeS{k1}))
                    [obj.StackInfo(RangeD{k1}).HumanCountSpots] = num{:}
                end
                
                if k1 == 5
                    num1 = num2cell(xlsread(excelfile, 2, 'N3:N13'))
                    num2 = num2cell(xlsread(excelfile, 2, 'Q3:Q16'))
                    [obj.StackInfo(1:11).HumanCountSpots] = num1{:}
                    [obj.StackInfo(12:25).HumanCountSpots] = num2{:}
                end
                
                if k1 == 6
                    num1 = num2cell(xlsread(excelfile, 2, 'T3:T11'))
                    num2 = num2cell(xlsread(excelfile, 2, 'W3:W15'))
                    num3 = num2cell(xlsread(excelfile, 2, 'Z3:Z10'))
                    num4 = num2cell(xlsread(excelfile, 2, 'AC3:AC14'))
                    num5 = num2cell(xlsread(excelfile, 2, 'AF3:AF14'))
                    
                    [obj.StackInfo(1:9).HumanCountSpots] = num1{:}
                    [obj.StackInfo(10:21).HumanCountSpots] = num2{:}
                    [obj.StackInfo(23:30).HumanCountSpots] = num3{:}
                    [obj.StackInfo(31:42).HumanCountSpots] = num4{:}
                    [obj.StackInfo(43:53).HumanCountSpots] = num5{:}
                end
                
                if k1 == 7
                    
                    num1 = num2cell(xlsread(excelfile, 2, 'AI3:AI18'))
                    num2 = num2cell(xlsread(excelfile, 2, 'AL3:AL13'))
                    num3 = num2cell(xlsread(excelfile, 2, 'AO3:AO18'))
                    
                    [obj.StackInfo(1:16).HumanCountSpots] = num1{:}
                    [obj.StackInfo(17:27).HumanCountSpots] = num2{:}
                    [obj.StackInfo(28:43).HumanCountSpots] = num3{:}
                end
                
                
                if k1 == 8
                    
                    num1 = num2cell(xlsread(excelfile, 2, 'AR3:AR17'))
                    num2 = num2cell(xlsread(excelfile, 2, 'AU3:AU13'))
                    num3 = num2cell(xlsread(excelfile, 2, 'AX3:AX17'))
                    num4 = num2cell(xlsread(excelfile, 2, 'BA3:BA9'))
                    %num5 = num2cell(xlsread(excelfile,2,'BD3:BD10'))
                    %num6 = num2cell(xlsread(excelfile,2,'BF3:BF5'))
                    
                    [obj.StackInfo(1:15).HumanCountSpots] = num1{:}
                    [obj.StackInfo(16:26).HumanCountSpots] = num2{:}
                    [obj.StackInfo(27:41).HumanCountSpots] = num3{:}
                    [obj.StackInfo(42:48).HumanCountSpots] = num4{:}
                    
                end
                
                if k1 == 9
                    num1 = num2cell(xlsread(excelfile, 2, 'BJ3:BJ13'))
                    [obj.StackInfo(1:11).HumanCountSpots] = num1{:}
                end
                
                if k1 == 10
                    num1 = num2cell(xlsread(excelfile, 2, 'BM3:BM15'))
                    [obj.StackInfo(1:13).HumanCountSpots] = num1{:}
                end
                
                if k1 == 11
                    disp('this batch had unknown annotated data')
                end
                
                for k2 = 1:numel(obj.StackInfo)
                    if isempty(obj.StackInfo(k2).HumanCountSpots)
                        obj.StackInfo(k2).HumanCountSpots = NaN
                    end
                end
                
                obj.saveit
            end
        end
        function restant 
            %{
 
         for k1 = 2:11
             k1
             load([P{k1}, '/', S{k1}, '.mat']);
             obj.ShowFish([], true)
         end
            %}
            
            [folder, name, ext] = fileparts(which('obj'));
            S = dbstack('-completenames');
            S(1).file
            
            
            %
            for k2 = 1:8
                
                load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected')
                %load([P{1},'/',I{1},'.mat'], 'SpotsDetected')
                
                nspots = []
                for k1 = 1:numel(SpotsDetected)
                    nspots(k1, 1) = numel(SpotsDetected{k1});
                end
                
                Sheet = [{obj.StackInfo.stackname}', {obj.StackInfo.stacksize}', num2cell(nspots)]
                title = {obj.InfoName, 'images', 'nspots'}
                
                ds = cell2dataset([title; Sheet]);
                export(ds, 'file', [P{k2}, '/', I{k2}, '.csv'], 'delimiter', ',')
            end
            %
            %
            %writes to text file
            fileID = fopen([P{1}, '/', I{1}, '.txt'], 'w');
            %fprintf(fileID,'%s %u\n',obj.StackInfo.stackname, nspots);
            fprintf(fileID, '%s %s %s \n', 'filename', 'images', 'spots');
            for k = 1:12
                fprintf(fileID, '%s %u %u \n', obj.StackInfo(k).stackname, obj.StackInfo(k).stacksize, nspots(k));
            end
            fclose(fileID)
            %
            
            %
            does not work
            fid = fopen([P{1}, '/', I{1}, '.csv'], 'w');
            fprintf(fid, '%s,', Sh{1, 1:end-1});
            fprintf(fid, '%s\n', Sh{1, end});
            fclose(fid);
            dlmwrite('test.csv', Sh(2:end, :), '-append');
            %
        end
    end
end