function PImg = sng_Im2Polar3(Img,bs)
%transforms image to polar images
% Img should be 1000 x 1000
% PImg = 500 X pi * 300

%{
Img =IsquareMB
%}


%{
figure;imagesc(Img)
Img = Img(:,501:1500,:);
figure;imagesc(Img)
%}

%{
%generate subscript matrix (mesh)
s = size(R);
[X,Y] = meshgrid(1:s(2),1:s(1));
%centercoords
cen = s(1:2)/2; %center
%make center zero
X2 = X - cen(2);
Y2 = Y - cen(1);


%% image to polar
Radius = sqrt(X2.^2+Y2.^2);
Theta = atan2(Y2,X2); %branchcut on the left -pi -> pi

%}

s = size(Img);
cen = s(1:2)/2; %center


%mesh for the interpolation
ra = linspace(0,499,500);
th = linspace(-pi,pi,round(pi*300));
[Ra,Th] = meshgrid(ra,th);

X3 = Ra.*cos(Th)+cen(2);
Y3 = Ra.*sin(Th)+cen(1);

%nearest neigbour, can introduce bilinear
X4 = round(X3);
Y4 = round(Y3);

%bilinear



Q = sub2ind([1000,1000],Y4(:),X4(:));

if ndims(Img) == 2
    PR = Img(Q);
    PImg = reshape(PR,numel(th),numel(ra),1);
elseif ndims(Img) == 3 && s(3) == 3
    R = Img(:,:,1);
    G = Img(:,:,2);
    B = Img(:,:,3);

    PR = R(Q);
    PG = G(Q);
    PB = B(Q);

    PImg = reshape([PR;PG;PB],numel(th),numel(ra),3);
end


%{
figure;imagesc(PImg)
axis equal tight
view(gca,[90 -90]);
%}










