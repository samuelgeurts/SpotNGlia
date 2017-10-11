
PathsComplete('sp')

%SpotPath = uigetdir(SpotPath)

Validation = BrainSpotValidation(SpotPath)

Validation = BrainVal(Validation)
Validation = SpotVal(Validation)
Validation = SpotBrainVal(Validation)

Validation.show
Validation.showgraph
showbrain(Validation,1:50)


validation2 = Validation;
Validation2.saveit

%%

tic
obj = BrainSpotValidation
toc
tic
obj.saveit
toc



obj = BrainVal(obj)
obj = SpotVal(obj)
obj = SpotBrainVal(obj)
obj.show

%%

obj2 = BrainSpotValidation
%obj = BrainVal(obj)
%obj.show(1)
obj2 = SpotVal(obj)
obj2.show

