%select file location
obj = SpotNGlia

%compute stacks
obj = obj.SliceCombination

%check the following variables before proceeding

openvar('obj.ImageInfo')
openvar('obj.StackInfo')

%run all algorithms and save final images
%when preprocession is completed, half of the time is passed

obj = obj.CompleteProgram

%Show fish and correct brainfishbrains
obj.CheckBrain
