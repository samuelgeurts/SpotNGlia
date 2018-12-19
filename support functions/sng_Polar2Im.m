function Img = sng_Polar2Im(PImg)
%transforms polar image to image using nearest neighbour interpolation
% Img should be 1000 x 1000
% PImg = 500 X pi * 300
%
%   16-02-2017  if loop to select for single or color channel transform

%{
Img=I80;
figure;imagesc(Img)
PImg = sng_Im2Polar3(Img(:,401:1400,:));
figure;imagesc(PImg)
%}

%we start with a matrix 942*500
%which is transformed to a circle in a 1000x1000 square
%we want to know where the pixels in the 1000x1000 image come from

%generate subscript matrix (mesh)
sp = size(PImg)
s = [1000,1000]

[X,Y] = meshgrid(1:s(2),1:s(1));
%centercoords
cen = s(1:2)/2; %center
%make center zero
X2 = X - cen(2);
Y2 = Y - cen(1);

R = sqrt(X2.^2+Y2.^2);
T = atan2(Y2,X2); %branchcut on the left -pi -> pi
%{
min(T(:))
max(T(:))
min(R(:))
max(R(:))
%}
%transform polar coordinates to matrix indices
T2 = (T + pi)/(2*pi) * (sp(1)-1)+1;
R2 = R + 1;

%{
min(T2(:))
max(T2(:))
min(R2(:))
max(R2(:))
%}

R3 = round(R2);
T3 = round(T2);

%{
min(T3(:))
max(T3(:))
min(R3(:))
max(R3(:))
%}

R4=R3;R4(R3 > 500)=500;
T4=T3;T4(R3 > 500)=500;

%{
min(T4(:))
max(T4(:))
min(R4(:))
max(R4(:))
%}

Q = sub2ind(sp,T4(:),R4(:));


if numel(sp) == 3
    PR = PImg(:,:,1);
    PG = PImg(:,:,2);
    PB = PImg(:,:,3);

    R = PR(Q);
    G = PG(Q);
    B = PB(Q);

    Img = reshape([R;G;B],1000,1000,3);
end

if numel(sp) == 2
    Temp = PImg(Q);
    Img = reshape(Temp,1000,1000);
end





%{
figure;imagesc(Img)
%}


end
