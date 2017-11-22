%select file location
obj = SpotNGlia

%compute stacks
obj = obj.SliceCombination

%check the following file before proceeding

openvar('obj.ImageInfo')
openvar('obj.StackInfo')

%run all algorithms ans save final images
%when preprocession is completed, half of the time is passed

obj = obj.CompleteProgram




%display fishes number 1 to 5
obj.ShowFish(1:5)

%display all fishes
obj.ShowFish

