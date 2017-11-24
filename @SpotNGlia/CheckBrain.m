function CheckBrain(obj, ifish)
%application to adjust brain region

h = waitbar(0, 'loading data', 'Name', 'Loading Data');

%fishnumbers = 1:numel(obj.StackInfo);
if ~exist('ifish', 'var') || (ifish > numel(obj.StackInfo))
    ifish = 1;
end

SpotParameters = [];
cxy = fliplr(obj.CompleteTemplate.CenterMidBrain);
nfishes = numel(obj.StackInfo);
checkupOrg(nfishes) = struct('Midbrain', [], 'Spots', []);

waitbar(0.25/7, h,'BrainSegmentationInfo')
load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo')
waitbar(2/7, h,'SpotsDetected')
load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected')
waitbar(3/7, h,'SpotParameters')
load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters');
waitbar(6.75/7, h,'checkup')


%fill adaptive variable containing initial brains ans spots
%checkupOrg stays unchanged, checkup is adaptive

for k = 1:nfishes
    checkupOrg(k).Spots = reshape([SpotsDetected{k}.Centroid], 2, numel(SpotsDetected{k}))'; %#ok<IDISVAR,USENS>
    checkupOrg(k).Midbrain = fliplr(BrainSegmentationInfo(k).BrainEdge);
end

load([obj.SavePath, '/', obj.InfoName, '.mat'], 'checkup')
if ~exist('checkup', 'var')
    checkup = checkupOrg;
    checkup(1).Corrections = [];
    [checkup.Include] = deal(true);    
end
waitbar(7/7, h,'Image')


%initialize figure with axes and plots
fh = figure;
sh = imshow(ones(obj.CompleteTemplate.Size), ...
    'Border', ...
    'tight', ...
    'InitialMagnification', 50);
hold on
ln = plot(0, 0, ...
    'Color', [255, 75, 75]/255, ...
    'LineWidth', 2);
sc = scatter(0, 0, 400, ...
    'LineWidth', 2, ...
    'MarkerEdgeColor', 1/255*[255, 75, 75]);
str = sprintf(['Fish: %.0f of %.0f', '\n', ...
    'Computed Spots: %.0f'], 0, 0, 0);

ah = annotation('textbox', [.05, .63, .7, .3], ...
    'String', str, ...
    'FitBoxToText', 'on', 'FontSize', 15);
ph = plot(0, 0, ...
    'bo');


%[XbOrg, YbOrg, XsOrg, YsOrg] = IniFish(ifish);
IniFish


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


%initialize Buttons
btn = uicontrol('Style', 'pushbutton', 'String', 'Correct Brain', ...
    'Position', [30, 10, 100, 22], ...
    'Callback', {@CorrectBrainButton});%#ok<NASGU>

btn2 = uicontrol('Style', 'pushbutton', 'String', 'Reset', ...
    'Position', [130, 10, 100, 22], ...
    'Callback', {@ResetButton});%#ok<NASGU>

sliderstep = 1 / (nfishes - 1);
sld = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', nfishes, 'Value', ifish, ...
    'Position', [30, 30, 200, 30], ...
    'Callback', @FishSlider, ...
    'SliderStep', [sliderstep, sliderstep], ...
    'Units', 'Normalized');

cbh = uicontrol('Style', 'checkbox', 'String', 'Include Fish', ...
    'Value', checkup(ifish).Include, 'Position', [240, 41, 100, 22], ...
    'Callback', @checkBoxCallback);

btn3 = uicontrol('Style', 'pushbutton', 'String', 'Save', ...
    'Position', [530, 10, 100, 22], ...
    'Callback', {@SaveButton});%#ok<NASGU>
uicontrol(sld);

delete(h)

%% correct button
    function CorrectBrainButton(~, ~) %evenveel + 2 variabelen als in uicontrol
        %f = figure(figh(figs))
        
        
        PolarN = BrainSegmentationInfo(ifish).PolarTransform;
        sp = fliplr(size(PolarN));
        si = size(PolarN);
        
        figure(fh)
        hold on
        rectangle('Position', [cxy(2) - 500, cxy(1) - 500, 1000, 1000]);
        
        but = 1;
        xy = checkup(ifish).Corrections;
        n = size(checkup(ifish).Corrections, 1);
        
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
            %set(ln, 'XData', X3)
            %set(ln, 'YData', Y3)
            %set(ln, 'Marker', 'none')
            
            ln.XData = X3;
            ln.YData = Y3;
            ln.Marker = 'none';
            
            %% compute new spots insite area
            [rc] = reshape([SpotParameters{ifish}.Centroid], 2, numel(SpotParameters{ifish}))';
            [in, ~] = inpolygon(rc(:, 1), rc(:, 2), X3, Y3);
            
            %temp = num2cell(in); [Regions1.Insite] = temp{:};
            
            SpotsDetec = SpotParameters{ifish}(in' == 1 & ...
                [SpotParameters{ifish}.LargerThan] == 1 & ...
                [SpotParameters{ifish}.SmallerThan] == 1 & ...
                [SpotParameters{ifish}.MinProbability] == 1);
            
            [rc] = reshape([SpotsDetec.Centroid], 2, numel(SpotsDetec))';
            
            %% set the new spots
            sc.XData = rc(:, 1);
            sc.YData = rc(:, 2);
            
        end
        updatecheckup
    end


%% reset button
    function ResetButton(~, ~) %evenveel + 2 variabelen als in uicontrol
        
        ph.XData = [];
        ph.YData = [];
        ln.XData = checkupOrg(ifish).Midbrain(:, 1);
        ln.YData = checkupOrg(ifish).Midbrain(:, 2);
        sc.XData = checkupOrg(ifish).Spots(:, 1);
        sc.YData = checkupOrg(ifish).Spots(:, 2);
        cbh.Value = true;

        updatecheckup
        
    end

%% fish slider
    function FishSlider(source, ~) %evenveel + 2 variabelen als in uicontrol
        ifish = round(source.Value);
        IniFish;
    end

%% checkbox include
    function checkBoxCallback(source, ~) %evenveel + 2 variabelen als in uicontrol
        checkup(ifish).Include = source.Value;
        
        if ~source.Value           
            ph.XData = [];
            ph.YData = [];
            ln.XData = [];
            ln.YData = [];
            sc.XData = [];
            sc.YData = [];

            updatecheckup
        else
            ResetButton
        end
    end

%% save button
    function SaveButton(~, ~) %evenveel + 2 variabelen als in uicontrol
        updatecheckup
        save([obj.SavePath, '/', obj.InfoName, '.mat'], 'checkup', '-append');
        obj.buildsheet
    end

%% update checkup variable
    function updatecheckup
        checkup(ifish).Midbrain = [ln.XData', ln.YData'];
        checkup(ifish).Spots = [sc.XData', sc.YData'];
        checkup(ifish).Corrections = [ph.XData', ph.YData'];
        checkup(ifish).Include = cbh.Value;
        uicontrol(sld);
    end

    function IniFish
        AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(ifish).stackname, '.tif']);
        
        ln.XData = checkup(ifish).Midbrain(:, 1);
        ln.YData = checkup(ifish).Midbrain(:, 2);
        sc.XData = checkup(ifish).Spots(:, 1);
        sc.YData = checkup(ifish).Spots(:, 2);
        sh.CData = uint8(AlignedFish);
        ah.String = sprintf(['Fish: %.0f of %.0f', '\n', ...
            'Computed Spots: %.0f'], ifish, nfishes, numel(sc.XData));
        if ~isempty(checkup(ifish).Corrections)
            ph.XData = checkup(ifish).Corrections(:, 1);
            ph.YData = checkup(ifish).Corrections(:, 2);
        else
            ph.XData = [];
            ph.YData = [];
        end
        
        cbh.Value = checkup(ifish).Include; 
    end
end