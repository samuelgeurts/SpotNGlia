%cut large mat file 'stackinfo' in parts



PathsComplete('sp')

load([SpotPath,'/','stackinfo.mat'])
%%
n = numel(stackinfo)

%[stackinfobase.stackname(1:50)] = {stackinfo(:).stackname}
%[stackinfobase.stacksize(1:50)] = {stackinfo.stacksize}
%[stackinfobase.imagenames(1:50)] = {stackinfo.imagenames}

%%

temp = {stackinfo.stackname};[stackbase(1:n).stackname] = temp{:}
temp = {stackinfo.stacksize};[stackbase(1:n).stacksize] = temp{:}
temp = {stackinfo.imagenames};[stackbase(1:n).imagenames] = temp{:}

stackinfo_pp = [stackinfo(:).Preprocession]
stackinfo_ed = [stackinfo(:).ExtendedDeptOfField]

temp = {stackinfo.Registration};[stackinfo_rg(1:n).reg] = temp{:}
%or
for l = 1: numel(stackinfo(1).Registration)
    for k = 1:n
    stackinfo_rg(k).(stackinfo(k).Registration(l).name) = stackinfo(k).Registration(l).value
    end
end

stackinfo_br = [stackinfo(:).BrainSegmentation]

temp = {stackinfo.SpotSelection};[stackinfo_sp(1:n).sp] = temp{:}
%or
for l = 1: numel(stackinfo(1).SpotSelection)
    for k = 1:n
    stackinfo_sp(k).(stackinfo(k).SpotSelection(l).name) = stackinfo(k).SpotSelection(l).value
    end
end
