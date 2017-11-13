function Mask = BrainMask(obj, fishnumbers)
%generates brainmask from annotations and store

if ~exist('fishnumbers', 'var')
    fishnumbers = 1:numel(obj.StackInfo);
elseif max(fishnumbers) > numel(obj.StackInfo)
    error('at least one fish does not exist in RegistrationInfo')
end

nfishes = numel(fishnumbers);

if nargout == 0
    MaskDir = uigetdir;
else
    Mask = cell(nfishes,1);
end


for k1 = 1:nfishes
    fn = fishnumbers(k1);
    
    AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
    s = size(AlignedFish);
    bc = obj.Annotations(fn).MidBrain;
    BrainMask = poly2mask(bc(:, 1), bc(:, 2), s(1), s(2));
    
    if nargout == 0
         imwrite(BrainMask, [MaskDir, '/', obj.StackInfo(fn).stackname, '.tif'])
    else
        Mask{k1} = BrainMask;
    end
    
end

