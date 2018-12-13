function sng_zfinputAssign(zfinput,fieldstring)
%assign zfinput.value to zfinput.name in workspace

%fieldstring = 'BrainSegmentation'
%sng_zfinputAssign(zfinput,fieldstring)

temp = zfinput(strcmp({zfinput.stage},fieldstring)); 

for k = 1:numel(temp)

    assignin('caller',temp(k).name,temp(k).value)
    
end

