function CheckBrain(obj, fishnumbers)
%application to adjust brain region

if ~exist('fishnumbers', 'var')
    fishnumbers = 1:numel(obj.StackInfo);
elseif max(fishnumbers) > numel(obj.StackInfo)
    error('at least one fish does not exist in RegistrationInfo')
end

nfishes = numel(fishnumbers);

SpotParameters = [];
k1 = 1;
n = 0;
xy = [];

load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo')
load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected')

cxy = fliplr(obj.CompleteTemplate.CenterMidBrain);

%initialize figure with axes and plots
fh = figure;
sh = imshow(ones(obj.CompleteTemplate.Size),...
    'Border',...
    'tight',...
    'InitialMagnification', 50);
hold on
ln = plot(0,0,...
    'Color', [255, 75, 75]/255,...
    'LineWidth', 2);
sc = scatter(0,0, 400,...
    'LineWidth', 2,...
    'MarkerEdgeColor', 1/255*[255, 75, 75]);
str = sprintf(['Fish: %.0f of %.0f','\n',...
    'Computed Spots: %.0f'],0,0,0)

ah = annotation('textbox', [.05, .63, .7, .3],...
    'String', str,...
    'FitBoxToText', 'on', 'FontSize', 15);
ph = plot(0, 0,...
    'bo');



[XbOrg, YbOrg, XsOrg, YsOrg] = IniFish(fishnumbers(k1));


    function [XbOrg, YbOrg, XsOrg, YsOrg] = IniFish(fn)
        
        %figure handles
        %obj.ShowFish(fn);
        
        cmbr = BrainSegmentationInfo(fn).BrainEdge;
        AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(fn).stackname, '.tif']);
        cmbs = reshape([SpotsDetected{fn}.Centroid], 2, numel(SpotsDetected{fn}))';
        
        sh.CData = uint8(AlignedFish);
        ln.XData = cmbr(:, 2);
        ln.YData = cmbr(:, 1);
        sc.XData = cmbs(:, 1);
        sc.YData = cmbs(:, 2);              
        ah.String = sprintf(['Fish: %.0f of %.0f','\n',...
            'Computed Spots: %.0f'],k1,nfishes,size(cmbs, 1))
        ph.XData = [];
        ph.YData = [];
        
        %story original brain and spot values for reset function
        XbOrg = get(ln, 'XData');
        YbOrg = get(ln, 'YData');
        XsOrg = get(sc, 'XData');
        YsOrg = get(sc, 'YData');
        drawnow
    end

load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters');


%buttonfigure parameters
sizb = 70;
barcor = 22; %size of the close-minimalization-figure bar


%{
        figh = findobj('type', 'figure');
        figs = numel(figh);
 
        %error if no figure is available
        if figs == 0
            error('no figure present')
        elseif figs == 1
            fignumber = get(figh, 'Number');
        else
            fignumber = cell2mat(get(figh, 'Number'));
            [~, order] = sort(fignumber);
        end
%}

%%


%create a figure below the place a default figure is placed
%figurenumber is 1 higher than max number of found figures

%posf = figh(end).Position
posf = fh.Position;
posb = [posf(1), posf(2) - sizb - barcor, posf(3), sizb];

uifignumber = fh.Number + 1;
figure(uifignumber)
set(gcf, 'position', posb);
set(gcf, 'Toolbar', 'none');
set(gcf, 'Menubar', 'none');


%[x, y] = ginput(1)

%initialize Correct Brain button
btn = uicontrol('Style', 'pushbutton', 'String', 'Correct Brain', ...
    'Position', [30, 30, 100, 22], ...
    'Callback', {@CorrectBrainButton}); %#ok<NASGU>
%initialize Correct Brain button

btn2 = uicontrol('Style', 'pushbutton', 'String', 'Reset', ...
    'Position', [130, 30, 100, 22], ...
    'Callback', {@ResetButton}); %#ok<NASGU>

btn3 = uicontrol('Style', 'pushbutton', 'String', 'Save and proceed', ...
    'Position', [230, 30, 100, 22], ...
    'Callback', {@SaveButton});

uicontrol(btn3);

%% correct button
    function CorrectBrainButton(~, ~) %evenveel + 2 variabelen als in uicontrol
        %f = figure(figh(figs))
        fn = fishnumbers(k1);
        
        

        PolarN = BrainSegmentationInfo(fn).PolarTransform;
        sp = fliplr(size(PolarN));
        si = size(PolarN);
        
        figure(fh)
        hold on
        
        rectangle('Position', [cxy(2) - 500, cxy(1) - 500, 1000, 1000]);
        
        
        but = 1;
        
        %% pick input and display
        while but == 1
            [x, y, but] = ginput(1); %pick a new input
            
            if but == 1
                n = n + 1;
                xy(n, :) = [x, y];
                ph.XData = xy(:, 1);
                ph.YData = xy(:, 2);
            end
        end
        
        if ~isempty(xy)
            %% transform coordinates to polar
            coord = fliplr(xy); %(y,x)
            coord2 = coord - cxy; %set center to cxy
            %transform to polar coordinates
            [T, R] = cart2pol(coord2(:, 2), coord2(:, 1));
            %transform polar coordinates to pixel coordinates
            T2 = (T + pi) / (2 * pi) * (sp(1) - 1) + 1;
            R2 = R + 1;
            
            %% shortest path
            [~, I2, J2] = sng_ShortestPath(PolarN, round([R2, T2]));
            
            %% transform path to aligned fish
            Rpol = I2 - 1;
            Tpol = 2 * pi * J2 / sp(1) - pi;
            [X, Y] = pol2cart(Tpol, Rpol);
            X2 = X + si(2) / 2;
            Y2 = Y + si(1) / 2;
            X3 = X2 + cxy(2) - si(2) / 2; %set center to cxy
            Y3 = Y2 + cxy(1) - si(1) / 2;
            
            %{
             AlignedFish = imread([obj.SavePath,'/','AlignedFish','/',obj.StackInfo(fn).stackname,'.tif']);
             figure;imagesc(AlignedFish)
             hold on
             scatter(coord(2),coord(1))
 
             X = coord2(2)
             Y = coord2(1)
             [Isquare] = sng_boxaroundcenter(AlignedFish,fliplr(cxy));
             Isquare(Y-5:Y+5,X-5:X+5,1) = 20;
             figure;imshow(uint8(Isquare))
 
             Ipolar = sng_Im2Polar3(Isquare);
             figure;imshow(uint8(Ipolar))
             sp = size(Ipolar)
             si = size(Isquare)
 
               figure;imshow(uint8(Isquare))
               hold on
               plot(X2,Y2)
            %}
            
            %% set the new line object
            set(ln, 'XData', X3)
            set(ln, 'YData', Y3)
            set(ln, 'Marker', 'none')
            
            %% compute new spots insite area
            [rc] = reshape([SpotParameters{fn}.Centroid], 2, numel(SpotParameters{fn}))';
            [in, ~] = inpolygon(rc(:, 1), rc(:, 2), X3, Y3);
            
            %temp = num2cell(in); [Regions1.Insite] = temp{:};
            
            SpotsDetec = SpotParameters{fn}(in' == 1 & ...
                [SpotParameters{fn}.LargerThan] == 1 & ...
                [SpotParameters{fn}.SmallerThan] == 1 & ...
                [SpotParameters{fn}.MinProbability] == 1);
            
            [rc] = reshape([SpotsDetec.Centroid], 2, numel(SpotsDetec))';
            
            set(sc, 'XData', rc(:, 1))
            set(sc, 'YData', rc(:, 2))
        end
        
        figure(uifignumber);
    end


%% reset button
    function ResetButton(~, ~) %evenveel + 2 variabelen als in uicontrol
        
        n = 0;
        xy = [];
        ph.XData = [];
        ph.YData = [];
        
        set(ln, 'XData', XbOrg);
        set(ln, 'YData', YbOrg);
        set(sc, 'XData', XsOrg);
        set(sc, 'YData', YsOrg);
    end

%% save button
    function SaveButton(~, ~) %evenveel + 2 variabelen als in uicontrol

        %%SAVE
        
        n = 0;
        xy = [];
        if k1 < nfishes
            k1 = k1 + 1;
            [XbOrg, YbOrg, XsOrg, YsOrg] = IniFish(fishnumbers(k1));
        end
        
        
    end


%{
         fig = uifigure('Position',posb);
         btn = uibutton(fig,'push',...
             'Position',[30 30 100 22],...
             'ButtonPushedFcn', @(btn,event) CorrectBrainPushed(btn));
%}

%{
function CorrectBrainPushed(btn)
         disp('push')
         [x, y] = ginput(1)
         figure(sliderfignumber);
         end
%}

end