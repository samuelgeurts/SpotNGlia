function [Icombined,edoutput] = ExtendedDeptofFieldLink2(CorrectedSlice,zfinput)

%version ExtendedDeptofFieldLink2
%       make usable for zfinput parameter system

if ~exist('zfinput','var');
    zfinput = struct
    zfinput = sng_zfinput(zfinput,0,'ExtendedDeptOfField','edof','variancedisksize',7,'?'); %sigma for bandpassfilter (lowpass)    
end

sng_zfinputAssign(zfinput,'ExtendedDeptOfField')

[IndexMatrix, variance_sq, Icombined] = sng_StackDOF2(CorrectedSlice,variancedisksize);

edoutput.IndexMatrix = uint8(IndexMatrix);
edoutput.variancedisksize = uint8(variance_sq);


%{
sz = size(CorrectedSlice{1});
for k = 1:numel(CorrectedSlice)
figure;imagesc(uint8(CorrectedSlice{k}))
axis off tight equal
set(gca,'position',[0 0 1 1],'units','normalized')
text(sz(2)*0.95,sz(1)*0.05,num2str(k),'FontSize',20)
end

for k = 1:numel(CorrectedSlice)
figure;imagesc(uint8(variance_sq(:,:,k)));
axis off tight equal
set(gca,'position',[0 0 1 1],'units','normalized')
text(sz(2)*0.95,sz(1)*0.05,num2str(k),'FontSize',20,'Color',[1 1 1])
end

SelectionImage = CorrectedSlice;
for k = 1:numel(CorrectedSlice)
SelectionImage{k}(repmat(slicen,1,1,3) ~= k) = 0;
figure;imagesc(uint8(SelectionImage{k}));
axis off tight equal
set(gca,'position',[0 0 1 1],'units','normalized')
    if k == 1;
        text(sz(2)*0.95,sz(1)*0.05,num2str(k),'FontSize',20,'Color',[0 0 0])
    else
        text(sz(2)*0.95,sz(1)*0.05,num2str(k),'FontSize',20,'Color',[1 1 1])
    end
end

figure;imagesc(Icombined)
set(gca,'position',[0 0 1 1],'units','normalized')
axis off tight equal


%}