function zfinput = sng_zfinput(zfinput,nrow,varargin)    
%stores parameters with name and other valuable values in a structure variable
%Ex:    zfinput = struct('stage',[],'substage',[],'name',[]','value',[],'sensitivity',[])
%       input = zfinput(input,1,'preprocession','bandpass filtering','sigma',1,'moderate')
%   nrow:   0       - for selecting nect row or replace existing variable 
%           1:inf   - select row

%
%give error if number of input arguments is not 7
narginchk(7,7);    

%check if zfinput first row is empty in case of existing fieldnames or not
if ~isfield(zfinput,'stage') || isempty([zfinput.stage])
    nrow = 1;    
%check if variable already exist, if so it is replaced, otherwise
%added at the next row
elseif strcmp(nrow,'next') || nrow == 0
    na = (find(strcmp({zfinput.stage},varargin{1})));
    nb = (find(strcmp({zfinput.substage},varargin{2})));
    nc = (find(strcmp({zfinput.name},varargin{3})));
    %result is overlapping numbers
    nt = intersect(na,nb);
    nrow = intersect(nt,nc);
    if isempty(nrow)
        nrow = numel(zfinput)+1;
    end
end

zfinput(nrow).stage = varargin{1};
zfinput(nrow).substage = varargin{2};
zfinput(nrow).name = varargin{3};
zfinput(nrow).value = varargin{4};
zfinput(nrow).sensitivity = varargin{5};

end
