function varargout = PathsComplete(varargin)


%{
% tp3 = new template path for 3 dpf
% tp = template
% or = original
% pp = preprocession
% ed = extended dept of field
% rg = registration
% br = brain
% sp = spot
% roib = brain roi
% rois = spot roi


varargin = {'rg','bp','pp'}

Example 1:
     PathsComplete
Example 2
     PathsComplete('rg','bp','pp')

%}

%if at erasmus
if strcmp(getenv('OS'),'Windows_NT') && strcmp(getenv('username'),'260018')
    Basepath = 'C:\Users\260018\Dropbox\20170327_3dpf_test';
    TemplatePath = 'C:\Users/260018\Dropbox\Template brain\3 dpf';
    TemplatePath3dpf = 'C:\Users/260018\Dropbox\Template 3 dpf';        
elseif strcmp(getenv('OS'),'Windows_NT') && strcmp(getenv('username'),'SNGeu')
    Basepath = 'C:\Users\SNGeu\Dropbox\20170327_3dpf_test';
    TemplatePath = 'C:\Users/SNGeu\Dropbox\Template brain\3 dpf';
    TemplatePath3dpf = 'C:\Users/SNGeu\Dropbox\Template 3 dpf';          
%if at mac laptop
elseif  strcmp(getenv('USER'),'samuelgeurts')
    Basepath = '/Users/samuelgeurts/Dropbox/20170327_3dpf_test';
    TemplatePath = '/Users/samuelgeurts/Dropbox/Template brain/3 dpf';
    TemplatePath3dpf = '/Users/samuelgeurts/Dropbox/Template 3 dpf';    
    
elseif  strcmp(getenv('USER'),'Anouk')
    Basepath = '/Users/Anouk/Dropbox/20170327_3dpf_test';
    TemplatePath = '/Users/Anouk/Dropbox/Template brain/3 dpf';
    TemplatePath3dpf = '/Users/Anouk/Dropbox/Template 3 dpf';  
else
    error('unknown platform and user, update PathsComplete.m');
end

%windows10x64
%if strcmp(getenv('OS'),'Windows_NT') || strcmp(getenv('username'),'SNGeu')
%    Basepath = '/Users/samuelgeurts/Dropbox/20170327_3dpf_test';
%end

if nargin == 0
    
    assignin('base', 'Basepath', Basepath)    
    assignin('base', 'TemplatePath', TemplatePath)
    assignin('base', 'TemplatePath3dpf', TemplatePath3dpf)    
    
    OriginalPath = [Basepath,'/original'];
    PreprocessionPath = [Basepath,'/1_preprocessed'];
    ExtendedDeptOfFieldPath = [Basepath,'/2_combined'];
    RegistratedPath = [Basepath,'/3_registrated'];
    BrainPath = [Basepath,'/4_brainsegmented'];
    SpotPath = [Basepath,'/5_spotdetected'];
    RoiMicrogliaPath = [Basepath,'/Roi microglia'];
    RoiBrainPath = [Basepath,'/Roi brain'];

    assignin('base', 'OriginalPath', OriginalPath)
    assignin('base', 'PreprocessionPath', PreprocessionPath)    
    assignin('base', 'ExtendedDeptOfFieldPath', ExtendedDeptOfFieldPath)
    assignin('base', 'RegistratedPath', RegistratedPath)    
    assignin('base', 'BrainPath', BrainPath)
    assignin('base', 'SpotPath', SpotPath)    
    assignin('base', 'RoiMicrogliaPath', RoiMicrogliaPath)
    assignin('base', 'RoiBrainPath', RoiBrainPath)    

else
    
    temp = find(strcmp(varargin,'bp'));
    if temp        
        if nargout == 0
            assignin('base', 'Basepath', Basepath) 
        end
        varargout{temp} = Basepath;
    end
    
    %if any(strcmp(varargin,'bp')) 
    %    assignin('base', 'Basepath', Basepath) 
    %end
    
    temp = find(strcmp(varargin,'tp'));
    if temp
        if nargout == 0
            assignin('base', 'TemplatePath', TemplatePath)
        end
        varargout{temp} = TemplatePath;        
    end
    
    temp = find(strcmp(varargin,'tp3'));
    if temp
        if nargout == 0
            assignin('base', 'TemplatePath3dpf', TemplatePath3dpf)
        end
        varargout{temp} = TemplatePath3dpf;      
    end
    
    %{
    if ~exist(Basepath,'dir')
        disp('base folder')
        Basepath = uigetdir('/Users/samuelgeurts/Dropbox/20170327_3dpf_test','base folder');
        disp('template folder')
        TemplatePath = uigetdir('/Users/samuelgeurts/Dropbox/20170327_3dpf_test','template folder');
    end
    %}
    
    temp = find(strcmp(varargin,'or'));
    if temp
        OriginalPath = [Basepath,'/original'];
        if nargout == 0
            assignin('base', 'OriginalPath', OriginalPath)
        end
        varargout{temp} = OriginalPath;        
    end    
    
    temp = find(strcmp(varargin,'pp')); 
    if temp
        PreprocessionPath = [Basepath,'/1_preprocessed'];
        if nargout == 0
            assignin('base', 'PreprocessionPath', PreprocessionPath)
        end
        varargout{temp} = PreprocessionPath;
    end
    
    temp = find(strcmp(varargin,'ed'));
    if temp
        ExtendedDeptOfFieldPath = [Basepath,'/2_combined'];
        if nargout == 0
            assignin('base', 'ExtendedDeptOfFieldPath', ExtendedDeptOfFieldPath)
        end
        varargout{temp} = ExtendedDeptOfFieldPath;
    end
    
    temp = find(strcmp(varargin,'rg'));
    if temp
        RegistratedPath = [Basepath,'/3_registrated'];
        if nargout == 0
            assignin('base', 'RegistratedPath', RegistratedPath)
        end
        varargout{temp} = RegistratedPath;
    end
    
    temp = find(strcmp(varargin,'br'));
    if temp
        BrainPath = [Basepath,'/4_brainsegmented'];
        if nargout == 0
            assignin('base', 'BrainPath', BrainPath)
        end
        varargout{temp} = BrainPath;
    end
    
    temp = find(strcmp(varargin,'sp'));
    if temp
        SpotPath = [Basepath,'/5_spotdetected'];
        if nargout == 0
            assignin('base', 'SpotPath', SpotPath)
        end
        varargout{temp} = SpotPath;
    end
    
    temp = find(strcmp(varargin,'rois'));
    if temp
        RoiMicrogliaPath = [Basepath,'/Roi microglia'];
        if nargout == 0
            assignin('base', 'RoiMicrogliaPath', RoiMicrogliaPath)
        end
        varargout{temp} = RoiMicrogliaPath;
    end
    
    temp = find(strcmp(varargin,'roib'));
    if temp
        RoiBrainPath = [Basepath,'/Roi brain'];
        if nargout == 0
            assignin('base', 'RoiBrainPath', RoiBrainPath)
        end
        varargout{temp} = RoiBrainPath;
    end

end





%{
mkdir(PreprocessionPath);
mkdir(ExtendedDeptOfFieldPath);
mkdir(RegistratedPath);
mkdir(BrainPath);
mkdir(SpotPath);
%}


