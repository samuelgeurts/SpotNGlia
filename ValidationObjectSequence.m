
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

