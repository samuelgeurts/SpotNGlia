function [Icombined,Parameters] = ExtendedDeptofFieldLink(CorrectedSlice)

if ~exist('zfinput','var');
    zfinput = struct
    zfinput = sng_zfinput(zfinput,0,'Registration','ExtendedDeptOfField','variancedisksize',7,'?'); %sigma for bandpassfilter (lowpass)    
end

rg = zfinput(find(strcmp({zfinput.stage},'Registration'))) 
%fsv = [fs.value];
for k3 = 1:numel(rg)
    eval([rg(k3).name,'_rg=',num2str(rg(k3).value)])
end

variancedisksize_rg = find(strcmp(rg.name,'variancedisksize'))




[Icombined, slicen, variance_sq] = sng_StackDOF(CorrectedSlice,variancedisksize);

Parameters.variancedisksize = variancedisksize;


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