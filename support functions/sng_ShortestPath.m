%{
INorm3 = PolarN;
RC = round([R2, T2]);

%}
%{
 cxy = fliplr(CompleteTemplate.CenterMidBrain)
 coord = [518, 701] %(y,x)
 si = [1000, 1000]
 sf = size(AlignedFish)
 sp = size(Ipolar)
 sn = size(INorm3)
%{
     AlignedFish = imread([obj.SavePath,'/','AlignedFish','/',obj.StackInfo(fn).stackname,'.tif']);
     figure;imagesc(AlignedFish)
     hold on
     scatter(coord(2),coord(1))
%}
 coord2 = coord - cxy %set center to cxy
%{
      X = coord2(2)
      Y = coord2(1)
      [Isquare] = sng_boxaroundcenter(AlignedFish,fliplr(cxy));
      Isquare(Y-5:Y+5,X-5:X+5,1) = 20;
      figure;imshow(uint8(Isquare))
 
      Ipolar = sng_Im2Polar3(Isquare);
      figure;imshow(uint8(Ipolar))
      sp = size(Ipolar)
      si = size(Isquare)
%}
 [T2, R2] = sng_cart2pol2sub(coord2, sp) %tranform to polar than to pixel c
 RC = round([R2, T2])
 RC(2, 1:2) = [200, 500]
%}


function [parameters, I2, J2] = sng_ShortestPath(INorm3, RC, RC2)
%RC = (Row, Col) where path has to go through from left to right
%
%   a single path through [Y,X] is computed
%       Example: [I2, J2] = sng_shortestpath(INorm3, [Y,X])
%   a single path through all points is computed
%       Example: [I2, J2] = sng_shortestpath(INorm3, [Y(:),X(:)])
%   a single path form (Y,X) to (Y2,X2) is computed
%       Example: [I2, J2] = sng_shortestpath(INorm3, [Y,X],[Y2,X2])
%   k paths form (Y(k),X(k)) to (Y2(k),X2(k)) are computed
%       Example: [I2, J2] = sng_shortestpath(INorm3, [Y(k),X(k)],[Y2(k),X2(k)])

%{
RC = round([R2, T2])
INorm3 = PolarN;
%}

%{
          figure;imagesc(INorm3);sng_imfix
          hold on
          scatter(RC(1,2),RC(1,1),'MarkerFaceColor',[1 0 0])
          scatter(RC(2,2),RC(2,1),'MarkerFaceColor',[1 0.5 0])
 
%}

%Ishift = INorm3(:, [Col:end, 1:Col - 1]);
Idouble = [INorm3, INorm3];
sn = size(INorm3);
sd = size(Idouble);

%TODO can be made quicker if not doubled image but something with indices/coordinates
[edges] = sng_SquareGridDiGraph(sd); %compute edges
DGraph = digraph([edges(:, 1)], [edges(:, 2)], Idouble(edges(:, 2))); %compute directer graph


if exist('RC2', 'var')
    a = arrayfun(@(x, y) sub2ind(sd, x, y), RC(:, 1), RC(:, 2));
    b = arrayfun(@(x, y) sub2ind(sd, x, y), RC2(:, 1), RC2(:, 2));
    path = zeros(numel(a), sn(2)+1);
    
    Lengthd = zeros(1,numel(a));
    Area(1,k) = zeros(1,numel(a));
    Hor(1,k) = zeros(1,numel(a));
    Ver(1,k) = zeros(1,numel(a));
    
        for k = 1:numel(a)
        %assumed that the left and right site are connected
        if b(k) <= a(k)
            b(k) = b(k) + prod(sn);
        end
        [path(k, :), Lengthd(k)] = shortestpath(DGraph, a(k), b(k), 'Method', 'acyclic');
        
        [I, J] = ind2sub(sd, path);
        Area(k) = sum(pi*(I.^2)/1000); %<---?
        Hor(k) = I(end) + I(end/2);
        Ver(k) = I(round((3 * sn(2))/4)) + I(round((sn(2))/4));
    end
    
    
else
    %create path stop sequence
    inx = arrayfun(@(x, y) sub2ind(sn, x, y), RC(:, 1), RC(:, 2));
    inx = sort(inx);
    inx = [inx; inx(1) + prod(sn)];
    
    a = inx(1:end-1);
    b = inx(2:end);
    path = cell(1,numel(a));
    for k = 1:numel(a)
        [path{k}, Lengthd] = shortestpath(DGraph, a(k), b(k), 'Method', 'acyclic');
        if isempty(path{k})
            path{k} = b(k);
        else
            path{k} = path{k}(2:end);
        end
    end
    fullpath = cell2mat(path);

    
    
    [I, J] = ind2sub(sd, fullpath);
    
    %shift the line to the orignal image dimensions
    cut = find(J-sn(2) == 1);    
    if isempty(cut)
       %add addition points at end and begin of the image if an uncomputable path crosses the edge 
       Ie = round(I(end-1) + (sn(2)-J(end-1))*((I(end) - I(end-1))/(J(end) - J(end-1))));
       Ib =round(I(end-1) + (sn(2)-J(end-1)+1)*((I(end) - I(end-1))/(J(end) - J(end-1))));
       Je = sn(2);
       Jb = sn(2) + 1;       
       I2 =  [Ib I(end) I(1:end-1) Ie ];
       J2 =  [Jb-sn(2) J(end)-sn(2) J(1:end-1) Je ];
    else        
        I2 = [I(cut:end), I(1:cut-1)];
        J2 = [J(cut:end)-sn(2), J(1:cut-1)];
    end
    %J2 = 1:sn(2);
    
    Area = sum(pi*(I2.^2)/1000);
    Hor = I2(end) + I2(round(numel(I2)/2));
    Ver = I2(round((3 * sn(2))/4)) + I2(round((sn(2))/4));
end


%{
         figure;imagesc(INorm3);sng_imfix
         hold on
         s=scatter(RC(:,2),RC(:,1),'MarkerFaceColor',[1 0 0],'SizeData', 100)
%}
%{
         figure;imagesc(Idouble);sng_imfix
         hold on
         s=scatter(RC(:,2),RC(:,1),'MarkerFaceColor',[1 0 0],'SizeData', 100)
         plot(J,I,'r','LineWidth',2)
         drawnow
%}
%{
         figure;imagesc(INorm3);sng_imfix
         hold on
         plot(J2,I2,'r','LineWidth',2)
%}


parameters = [Lengthd, Area, Hor, Ver];


end

%{
 %transform to pixels comform to theta from -pi to pi
 
 Rpol = I2 - 1
 Tpol = 2 * pi * J2 / sp(1) - pi
 %Tpol2 = linspace(-pi,pi,sp(1))
 
 [X, Y] = pol2cart(Tpol, Rpol)
 
 X2 = X + si(1) / 2
 Y2 = Y + si(2) / 2
 
%{
      figure;imshow(uint8(Isquare))
      hold on
      plot(X2,Y2)
%}
 
 
 X3 = X2 + cxy(2) - si(1) / 2 %set center to cxy
 Y3 = Y2 + cxy(1) - si(1) / 2
 
%{
      figure;imshow(uint8(AlignedFish))
      hold on
      plot(X3,Y3)
%}
%}
