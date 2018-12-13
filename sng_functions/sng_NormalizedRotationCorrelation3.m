function [theta,varargout] = sng_NormalizedRotationCorrelation3(I30,angles1,centerofmass,varargin)
%
%Version sng_NormalizedRotationCorrelation3 20170616
%   change in/output for better use in rgoutput
%
%shifts image to center and calculated normalized correlation between the
%image and a flipped and rotated version of itself
%varargin{1} the rotated images specified with number
%
%
%Example [theta,maxangledeg,peaksep,ncorr] = sng_NormalizedRotationCorrelation3(I30,angles,[xcom2,ycom2])
%Example [theta,maxangledeg,peaksep,ncorr,I40,I50] = sng_NormalizedRotationCorrelation3(I30,angles1,[xcom2,ycom2])
%Example [theta,maxangledeg,peaksep,ncorr,I40,I50,I61] = sng_NormalizedRotationCorrelation3(Img2,angles1,[xcom2,ycom2],1:10:100)
%
%

%{
anglesteps = 100
centerofmass = [xcom2,ycom2]
I30 = Img2;
nargin = 3
%}

if nargin == 4 
    I61 = cell(numel(varargin{1}),1);
end

shift = (fliplr(size(I30(:,:,1))/2)-centerofmass);

%translates and crops 1 times the shift ->maybe more accurate
I40 = imtranslate(I30,shift,'OutputView','same');
%no translation and crop 2 times the shift -> smaller images
%I40 = imcrop(I30, [abs(shift),-abs(shift)] - [shift, shift] + [0 0 size(I30(:,:,1)));

%make circle mask
cm = sng_MaxCircleMask(I40);


%{
figure;imagesc(c_mask)
figure;imagesc(I40)
figure;imagesc(I40.*uint8(cm))
%}

I50 = fliplr(I40);

%{
figure;imagesc(I40norm);colorbar
figure;imagesc(inproductI40)
figure;imagesc(I50.*uint8(cm))

%}

ncorr = zeros(numel(angles1),1);
jj = 1;
for j=1:numel(angles1) 
    %I60 = imrotate(I50,angles(j),'bilinear','crop');
    
    I60 = sng_imrotate(I50,angles1(j));
    ncorr(j) = sng_NCCoi(I40,I60,cm);
    
    if nargin == 4 && ismember(j,varargin{1})
        I61{jj}=I60;
        jj = jj+1;
    end    
end

[~,maxcorr] = (max(ncorr));
theta = -(angles1(maxcorr)/2+(pi/2));

%% calculates the relative difference of the first and second maxima
% to give an idication how sure we are of the maximum peak
ncorrmax = imregionalmax(ncorr).*ncorr;
corrpeaks = sort(nonzeros(ncorrmax),'descend');

if numel(corrpeaks) >= 2
    peaksep = 1-abs(corrpeaks(2)/corrpeaks(1));
else
    peaksep = 1-abs(min(ncorr)/corrpeaks(1));
end

maxangle = angles1(maxcorr);

if nargout >= 2;varargout{1} = maxangle;end %maximum angle in degrees
if nargout >= 3;varargout{2} = peaksep;end %relative difference two highes peaks
if nargout >= 4;varargout{3} = ncorr;end
if nargout >= 5;varargout{4} = I40;end
if nargout >= 6;varargout{5} = I50;end
if nargout >= 7;varargout{6} = I61;end

end

%{
figure;plot(ncorr)
for k = 1:numel(I61)
    figure;imagesc(uint8(I61{k}))
end



%}

