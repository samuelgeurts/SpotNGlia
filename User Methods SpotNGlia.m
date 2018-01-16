
%% most important functions

%create SpotNGlia object and select file location
obj = SpotNGlia
%compute stacks
obj = obj.SliceCombination

%check the following variables before proceeding
openvar('obj.ImageInfo')
openvar('obj.StackInfo')

%run all algorithms in order and saves in between
obj = obj.CompleteProgram

%Show fish and apply braincorrections and spotcorrections
obj = obj.CheckFish

%Save object location specified in obj.SavePath
obj.saveit

%% individual algorithm functions which has to be processed in order

obj = obj.PreProcession
obj = obj.ExtendedDeptOfField
obj = obj.Registration
obj = obj.BrainSegmentation
obj = obj.SpotDetection

%% supporting properties that can be set after the SpotNglia object is created

%dont show screen that alerts for stackinfo and imageinfo check [true or false]
obj.ImageInfoChecked_TF = true
%preset obj.Sorting ['Date' or 'Name']
obj.Sorting = 'Date'
%change delimiter to [';' or ',']
obj.Delimiter = ';' 

%% run multiple batches

for k = 1:5
    obj{k} = SpotNGlia
    obj{k}.ImageInfoChecked_TF = true
    obj{k}.Sorting = 'Date'
    obj{k} = obj{k}.SliceCombination
end

for k = 1:5
   obj{k} = obj{k}.CompleteProgram    
end


