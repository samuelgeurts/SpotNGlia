
PathsComplete('sp')

SpotPath = uigetdir(SpotPath)

Validation = BrainSpotValidation(SpotPath)

Validation = BrainVal(Validation)
Validation = SpotVal(Validation)
Validation = SpotBrainVal(Validation)

%Validation.show
%Validation.showgraph
showbrain(Validation,45,'Com')
showbrain(Validation,45,'Ann')

saveit(Validation,'Validation20171025')

%%

tic
obj = BrainSpotValidation()
toc
tic
%obj.saveit
toc



obj = BrainVal(obj)
obj = SpotVal(obj)
obj = SpotBrainVal(obj)
%obj.show

showbrain(obj,10,'Com')
showbrain(obj,10,'Ann')

%%

obj2 = BrainSpotValidation
%obj = BrainVal(obj)
%obj.show(1)
obj2 = SpotVal(obj)
obj2.show

