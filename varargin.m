function PathsComplete = varargin


%{
tp = template
pp = preprocession
ed = extended dept of field
rg = registration
br = brain
sp = spot
roib = brain roi
rois = spot roi

Example 1:
     PathsComplete
Example 2
     PathsComplete(rg bp pp)


%}
  
Basepath = '/Users/samuelgeurts/Dropbox/20170327_3dpf_test';

if strcmp(varargin,'bp')
    assignin('base', 'Basepath', Basepath)
end
    
if strcmp(varargin,'tp')
    TemplatePath = '/Users/samuelgeurts/Dropbox/Template brain/3 dpf';
    assignin('base', 'TemplatePath', TemplatePath)
end
    
%{    
if ~exist(Basepath,'dir')
    disp('base folder')
    Basepath = uigetdir('/Users/samuelgeurts/Dropbox/20170327_3dpf_test','base folder');
    disp('template folder')
    TemplatePath = uigetdir('/Users/samuelgeurts/Dropbox/20170327_3dpf_test','template folder');
end
%}

if strcmp(varargin,'pp')
    PreprocessionPath = [Basepath,'/1_preprocessed'];
    assignin('base', 'PreprocessionPath', PreprocessionPath)
end

if strcmp(varargin,'ed')
    ExtendedDeptOfFieldPath = [Basepath,'/2_combined'];
    assignin('base', 'ExtendedDeptOfFieldPath', ExtendedDeptOfFieldPath)
end

if strcmp(varargin,'rg')
    RegistratedPath = [Basepath,'/3_registrated'];
    assignin('base', 'RegistratedPath', RegistratedPath)
end
    
if strcmp(varargin,'br')
    BrainPath = [Basepath,'/4_brainsegmented'];
    assignin('base', 'BrainPath', BrainPath)
end
    
if strcmp(varargin,'sp')
    SpotPath = [Basepath,'/5_spotdetected'];
    assignin('base', 'SpotPath', SpotPath)
end
    
if strcmp(varargin,'rois')
    RoiMicrogliaPath = [Basepath,'/Roi microglia'];
    assignin('base', 'RoiMicrogliaPath', RoiMicrogliaPath)
end
    
if strcmp(varargin,'roib')
    RoiBrainPath = [Basepath,'/Roi brain'];
    assignin('base', 'RoiBrainPath', RoiBrainPath)
end
    
%{
mkdir(PreprocessionPath);
mkdir(ExtendedDeptOfFieldPath);
mkdir(RegistratedPath);
mkdir(BrainPath);
mkdir(SpotPath);
%}


