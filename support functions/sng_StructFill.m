function structurename = sng_StructFill(structurename,fieldcell,nrow)    
%fills a row of an existing structure variable with custum number of variables
%if nrow is not given, the next row is taken

%Example 1
%{
structurename = struct('stage',[],'substage',[],'name',[]','value',[])
structurename = sng_StructFill(structurename,{'stageA','substageB','tryout', 0}) 
%}
%Example 2
%{
structurename = struct('stage',[],'substage',[],'name',[]','value',[])
structurename = sng_StructFill(structurename,{'stageA','substageB','tryout', 0},10) 
%}

fieldnms = fieldnames(structurename);


%check if first row exists or is empty
if ~isfield(structurename,fieldnms{1}) || isempty([structurename.(fieldnms{1})])
    nrow = 1;    
elseif ~exist('nrow','var') || isempty(nrow)
    nrow = numel(structurename)+1;
end



for k10 = 1:numel(fieldcell)
    [structurename(nrow).(fieldnms{k10})] = fieldcell{k10};
end




