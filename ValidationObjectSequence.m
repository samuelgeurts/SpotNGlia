
PathsComplete('sp')

SpotPath = uigetdir(SpotPath)

obj = BrainSpotValidation(SpotPath)

obj = BrainVal(obj)
obj = SpotVal(obj)
obj = SpotBrainVal(obj)

obj.show(1)
%Validation.showgraph
showbrain(obj,1:50,'Com')
showbrain(obj,1:50,'Ann')

saveit(obj,'Validation20171025b')


obj2.show



%% SpotNGlia sequence


obj = SpotNGlia(0)

obj.Sorting = 'Name'
obj.ImageInfoChecked_TF = true;

obj = obj.SliceCombination([1,2,5])
obj = obj.PreProcession(1:2)