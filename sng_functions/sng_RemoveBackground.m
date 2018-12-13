function [Img,varargout] = sng_RemoveBackground(Img,method,varargin)
% this function works for a dark background
%only works for really dark or really light background

%method can be 'ElementThreshold','Afterpeak','Otsu','Triangle','TriangleSmpoth'
%Varargin has different values for different methods

%{
Img =Img(:,:,1);

method = 'TriangleSmooth'
varargin{1} = 1 %smooth

%}

%% generate boundary array
lines = 50; %how much lines from boundary?

top = Img(1:lines,:); %top lines
bottom = Img(end-lines+1:end,:); %bottom lines
left = Img(lines+1:end-lines,1:lines); %left lines minus corner
right = Img(lines+1:end-lines,end-lines+1:end); %right lines mines corner
boundary = double([top(:);bottom(:);left(:);right(:)]); %boundary array

if mean(boundary)>255/2
    boundary = 255-boundary;
    background = 'white';
else
    background = 'black';
end    

%% background removal methods

if strcmp(method,'ElementThreshold')
    sortb = sort(boundary);
    if isempty(varargin{1})
        value = 0.5;
    elseif varargin{1} >= 0 && varargin{1} <= 1
        value = varargin{2};
    else
        value = 0.5;
    end
    threshold = sortb(round(numel(boundary)*value));
end

% removal on base of a chosen value doesnt work as it is different for
% different images. Maybe only a minimum can be set.
% a max can be set as a maximum distance from the max

if strcmp(method,'AfterPeak')
    if isempty(varargin{1})
        value = 5;
    else
        value = varargin{1};
    end
    histb = histcounts(boundary,0:1:255);
    histbsmooth = smooth(histb);
    [~,c] = max(histbsmooth);
    %histdif = histbsmooth(1:end-1)-histbsmooth(2:end);
    %cutoff2 = find(histdif==0, 1, 'first');
    threshold = find(histbsmooth(ceil(c):end) <= value,1,'first')+ceil(c)-1;    
end

if strcmp(method,'Otsu')
    threshold = graythresh(Img)*255;
end

if strcmp(method,'Triangle')
    %value is the peak height
    %when the peak is at zero, better choose a lower peak value
    
    histb = histcounts(boundary,0:1:255);
    [mx,c_max] = max(histb);
    if isempty(varargin{1})
        value = mx;
    else
        value = varargin{1};
    end    
        
    c_end = find(histb ~= 0, 1, 'last');
    triangle = linspace(value,histb(c_end),c_end-c_max+1);
    histbint = histb(c_max:c_end);
    distance = triangle - histbint;
    [~,c_dis] = max(distance);
    threshold = c_dis + c_max-1;
end

if strcmp(method,'TriangleSmooth')
    %value is the smoothing width
    
    if isempty(varargin{1})
        value = 5;
    else
        value = varargin{1};
    end   
    %value = 1
    
    histb = histcounts(boundary,0:1:255);
   
    %smoothing method
    %the histogram is first extended to perform 'better' smoothing. But
    %when all energy is in the first peak, the maximum is shiften to the
    %right due to the average smoothing using zero for values outside the
    %border. But because we want initiale take a value besides the peak
    %value we dont have to solve for that. 
 
    extend(value+1:length(histb(:))+(value)) = histb(:);
    extendsmooth = smooth(extend',value);
    histbsmooth = extendsmooth(1+value:end-value);    
    
    
    [mx,c_max] = max(histbsmooth);

    c_end = find(histbsmooth ~= 0, 1, 'last');
    triangle = linspace(mx,histbsmooth(c_end),c_end-c_max+1);
    histbint = histbsmooth(c_max:c_end);
    distance = triangle - histbint';
    [~,c_dis] = max(distance);
    threshold = c_dis + c_max-1;
    
    %{
    figure;imagesc(Img)
    figure;bar(histb)
    figure;bar(histbsmooth)
    figure;bar(histbint)
    
    %}
    
end

%Values for the rhterhold and triangle
%the point at the threshold
c1 = threshold;
c2 = histbsmooth(threshold);
%the point at the peak
a2 = mx;
a1 = c_max;
%the point het the end of the triangle
b1 = c_end;
b2 = histbsmooth(c_end);
%point halfway the slope
d1 = (b1 - a1)/2 + a1;
d2 = (a2 - b2)/2 + b2;

Threshold = [c1,c2]';
Peak = [a1,a2]';
LastValue = [b1,b2]';
Halfway = [d1,d2]';

Table = table(Threshold,Peak,LastValue,Halfway,'RowNames',{'x','y'});





if strcmp(background,'black')
    Img(Img <= threshold) = 0;
end
if strcmp(background,'white')
    threshold = 255 - threshold;
    Img(Img >= threshold) = 255;
end

if nargout >= 2
    varargout{1} = threshold;
end
   
if nargout >= 3
    varargout{2} = histb;
    varargout{3} = histbsmooth;
end
if nargout >= 4
    varargout{4} = Table;
end



%playing around with creating a figure but decided not to finish it for my report
%Ik heb geprobeerd om de loodrechte lijn te berekenen die de
%gemaximaliseerd afstand aangeeft die de threshold bepaald. Maar dat is
%onzin omdat de loodrechtheid afhankelijk is van de schaal. Loodrechtheid
%is dus geen vereiste.





%{
%attempt to compute a perpendicular line. <= bullshit 
deltax =  c_end - c_max
deltay = histbsmooth(c_end) - mx
s = deltay/deltax %slope
%the point at the threshold

intersectx = - (a2 - c2 - (c1 / s) - (s*a1)) / ( (1/s) + s )
intersectx = (histbsmooth(threshold) - mx + c_max*(deltay/deltax) + threshold*(deltax/deltay))/(deltax/deltay + deltay/deltax)
intersecty = (deltay/deltax) * intersectx + mx
intersecty = (deltax/deltay) * intersectx + histbsmooth(threshold)
%}

%}
