function [Jaccard] = sng_Rgb2BwContrast(array1,array2,vector,imageTF)
%computes histogram of two colorarrays and computes contrast ratio

% vector = bestvec
bin = 255;

cv = combvec([0 255],[0 255],[0 255]);
lims = zeros(size(cv,2),1);
for k = 1:size(cv,2)
    lims(k) = dot(cv(:,k),vector);
end
lims2 = sng_lim(lims);

st = double(array1)*vector; %transformed spot
bt = double(array2)*vector; %transformed background

%st = st* 255/mx; %scaled to 255
%bt = bt* 255/mx;

st = (st - lims2(1)) * 255/(lims2(2)-lims2(1));
bt = (bt - lims2(1)) * 255/(lims2(2)-lims2(1));

SThist = histcounts(st,linspace(0,255,bin+1)); %histogram
BThist = histcounts(bt,linspace(0,255,bin+1));

STnorm = SThist/size(array1,1); %normalized histogram
BTnorm = BThist/size(array2,1);

if imageTF
    [minC2,~,~] = sng_threshold2(STnorm,BTnorm);                %find threshold
    sng_chistplot(STnorm,BTnorm,minC2);set(gcf,'numbertitle','off','name','RGB') 
end

Overlap = min(STnorm,BTnorm);
Combination = max(STnorm,BTnorm);
Jaccard = 1-sum(Overlap)/sum(Combination);

end