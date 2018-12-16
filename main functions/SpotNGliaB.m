classdef SpotNGliaB
    
    
    properties
        FishPath = cell(1)
        ObjPath = cell(1)
        ObjName = cell(1)
        Objectn = uint8(0)
        
        SavePath = []
        SaveName = []
    end
    
    properties(Transient = true)
        obj = SpotNGlia([])
    end
    
    methods
        
        function objb = SpotNGliaB
            %constructor
            objb = NewObj(objb);
        end
        
        function objb = NewObj(objb)
            %add paths to the objb
            pathname = [];
            newpath = true;
            while newpath
                n = numel(objb.FishPath);
                if isempty(objb.FishPath{1})
                    n = 0;
                end
                [pathname] = uigetdir(pathname, 'select object location');
                if isequal(pathname, 0)
                    disp('User pressed cancel');
                    newpath = false;
                else
                    objb.FishPath{n+1} = pathname;
                    Objectn = Objectn + 1;
                end
            end
        end
        
        function objb = getSNG(objb)
            %looks for SNG objects at the given path
            for k1 = 1:numel(objb.FishPath)
                dirlist = dir([objb.FishPath{k1}, filesep, 'SNG*.mat']);
                objb.ObjPath{k1} = dirlist.folder;
                objb.ObjName{k1} = dirlist.name;
            end
        end
        
        function objb = all(objb, SNGFunctionName, batchnumbers, input1)
            %executes SponNGlia function on all objects contained by SpotNGliaB object
            %objb = objb.all(@SpotVal)
                        
            if ~exist('batchnumbers', 'var') || isempty(batchnumbers)
                batchnumbers = 1:objb.Objectn;
            end
            
            for k1 = batchnumbers
                disp(num2str(k1))
                if (nargout == 0) && ~exist('input1', 'var')
                    SNGFunctionName(objb.obj(k1));
                elseif (nargout ~= 0) && ~exist('input1', 'var')
                    objb.obj(k1) = SNGFunctionName(objb.obj(k1));
                elseif (nargout == 0) && exist('input1', 'var')
                    SNGFunctionName(objb.obj(k1), input1);
                elseif (nargout ~= 0) && exist('input1', 'var')
                    objb.obj(k1) = SNGFunctionName(objb.obj(k1), input1);
                end
            end
            
            
        end
        
        function objb = allpar(objb, SNGFunctionName, batchnumbers, input1)
            %executes SponNGlia function on all objects contained by SpotNGliaB object
            %objb = objb.all(@SpotVal)
            %uses paralel computing
            
            if ~exist('batchnumbers', 'var') || isempty(batchnumbers)
                batchnumbers = 1:objb.Objectn;
            end
            
            %PARFOR OPTION           
            batchn = numel(batchnumbers);
            for k1 = 1:batchn
                tempobj(k1) = objb.obj(batchnumbers(k1));
            end
            
            if (nargout == 0) && ~exist('input1', 'var')
                parfor k1 = 1:batchn
                    disp(num2str(k1))
                    SNGFunctionName(tempobj(k1));
                end
            elseif (nargout ~= 0) && ~exist('input1', 'var')
                parfor k1 = 1:batchn
                    disp(num2str(k1))
                    tempobj(k1) = SNGFunctionName(tempobj(k1));
                end
            elseif (nargout == 0) && exist('input1', 'var')
                parfor k1 = 1:batchn
                    disp(num2str(k1))
                    SNGFunctionName(tempobj(k1), input1);
                end
            elseif (nargout ~= 0) && exist('input1', 'var')
                parfor k1 = 1:batchn
                    disp(num2str(k1))
                    tempobj(k1) = SNGFunctionName(tempobj(k1), input1);
                end
            end
                        
            for k1 = 1:batchn
                objb.obj(batchnumbers(k1)) = tempobj(k1);
            end
        end
        
        function objb = LoadObjects(objb)
            %stores SpotNGlia objects under the struct Object with increasing numbers
            %if Original = true, the files from the harddisk are taken (objb.MainFolderOriginal)
            for k1 = 1:numel(objb.FishPath)
                temp = load([objb.ObjPath{k1}, filesep, objb.ObjName{k1}]);
                objb.obj(k1) = temp.obj;
            end
        end
        
        function objb = saveit(objb)
            %opens file location and warn if overwrite
            [SaveName, objb.SavePath] = uiputfile([objb.SavePath, objb.SaveName, '.mat'], 'Save file name');
            [~, objb.SaveName] = fileparts(SaveName);
        end
        
        function quicksave(objb)
            %save and auto overwrite file
            save([objb.SavePath, objb.SaveName, '.mat'], 'objb')
        end
        %{
function custom1(obj)
             for k1 = 1:numel(objb.FishPath)
                 obj = objb.Object{k1}.obj;
             end
             for k1 = 1:numel(objb.FishPath)
 
                 obj = obj.NewPath(12);
                 obj = obj.SpotDetection2;
                 obj = obj.LoadAnnotations;
                 obj = obj.SpotVal;
 
                 objb.Object{k1}.obj = obj
             end
 
         end
        %}
        
        
    end
    
end
