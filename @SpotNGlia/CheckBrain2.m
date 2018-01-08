function CheckBrain2(obj, ifish)
%application to adjust brain region

h = waitbar(0, 'loading data', 'Name', 'Loading Data');

%fishnumbers = 1:numel(obj.StackInfo);
if ~exist('ifish', 'var') || (ifish > numel(obj.StackInfo))
    ifish = 1;
end
slice_end = cumsum([obj.StackInfo.stacksize]);
slice_start = [1, slice_end(1:end-1) + 1];
islice = slice_start(ifish);
subslice = 1;

SpotParameters = [];
cxy = fliplr(obj.CompleteTemplate.CenterMidBrain);
nfishes = numel(obj.StackInfo);
nslices = sum([obj.StackInfo.stacksize]);
checkupOrg(nfishes) = struct('Midbrain', [], 'Spots', []);


waitbar(0.0684, h, 'BrainSegmentationInfo')
load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo')
waitbar(0.1625, h, 'SpotsDetected')
load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotsDetected')
waitbar(0.2826, h, 'SpotParameters')
load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters');
waitbar(0.7227, h, 'RegistrationInfo')
load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo');
waitbar(0.7877, h, 'checkup')


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
waitbar(1, h, 'Image')


if isfield(obj.Annotations,'Spots')
    SpotAnn = true;
else
    SpotAnn = false;    
end
if isfield(obj.Annotations,'MidBrain')
    MidBrainAnn = true;
else
    MidBrainAnn = false;    
end

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
ln2 = plot(0, 0, ...
    'Color', [75, 75, 255]/255, ...
    'LineWidth', 2,...
    'Visible','off');
sc = scatter(0, 0, 400, ...
    'LineWidth', 2, ...
    'MarkerEdgeColor', 1/255*[255, 75, 75]);
sc2 = scatter(0, 0, 400, ...
    'LineWidth', 2, ...
    'MarkerEdgeColor', 1/255*[75, 75, 255],...
    'SizeData',300,...
    'Visible','off');
str = sprintf(['Fish: %.0f of %.0f', '\n', ...
    'Computed Spots: %.0f'], 0, 0, 0);

ah = annotation('textbox', [.05, .63, .7, .3], ...
    'String', str, ...
    'FitBoxToText', 'on', 'FontSize', 15);
ph = plot(0, 0, ...
    'bo');


%buttonfigure parameters
sizb = 70;
barcor = 22; %size of the close-minimalization-figure bar

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

rad = uicontrol('Style', 'radiobutton', 'String', 'Show Slices', ...
    'Position', [240, 10, 100, 22], ...
    'Value', 0, ...
    'Callback', {@SliceButton});%#ok<NASGU>

rad2 = uicontrol('Style', 'radiobutton', 'String', 'Compare', ...
    'Position', [340, 10, 100, 22], ...
    'Value', 0, ...
    'Callback', {@SliceButton2});%#ok<NASGU>

sliderstep = 1 / (nfishes - 1);
sld = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', nfishes, 'Value', ifish, ...
    'Position', [30, 30, 200, 30], ...
    'Callback', @FishSlider, ...
    'SliderStep', [sliderstep, sliderstep], ...
    'Units', 'Normalized');

sliderstep2 = 1 / (nslices - 1);
sld2 = uicontrol('Style', 'slider', ...
    'Min', 1, 'Max', nslices, 'Value', islice, ...
    'Position', [30, 30, 200, 30], ...
    'Callback', @SliceSlider, ...
    'SliderStep', [sliderstep2, sliderstep2], ...
    'Units', 'Normalized', ...
    'Visible', 'off');

cbh = uicontrol('Style', 'checkbox', 'String', 'Include Fish', ...
    'Value', checkup(ifish).Include, 'Position', [240, 41, 100, 22], ...
    'Callback', @checkBoxCallback);

btn3 = uicontrol('Style', 'pushbutton', 'String', 'Save', ...
    'Position', [530, 10, 100, 22], ...
    'Callback', {@SaveButton});%#ok<NASGU>
uicontrol(sld);

IniFish

delete(h)

    function CorrectBrainButton(~, ~)
        %correct button
        
        PolarN = BrainSegmentationInfo(ifish).PolarTransform;
        sp = fliplr(size(PolarN));
        si = size(PolarN);
        
        figure(fh)
        hold on
        rectangle('Position', [cxy(2) - 500, cxy(1) - 500, 1000, 1000]);
        
        but = 1;
        xy = checkup(ifish).Corrections;
        n = size(checkup(ifish).Corrections, 1);
        
        % pick input and display
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
            % transform coordinates to polar
            coord = fliplr(xy); %(y,x)
            coord2 = coord - repmat(cxy, size(coord, 1), 1); %set center to cxy
            
            %transform to polar coordinates
            [T, R] = cart2pol(coord2(:, 2), coord2(:, 1));
            %transform polar coordinates to pixel coordinates
            T2 = (T + pi) / (2 * pi) * (sp(1) - 1) + 1;
            R2 = R + 1;
            
            % shortest path
            [~, I2, J2] = sng_ShortestPath(PolarN, round([R2, T2]));
            
            % transform path to aligned fish
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
            
            ln.XData = X3;
            ln.YData = Y3;
            ln.Marker = 'none';
            
            % compute new spots insite area
            [rc] = reshape([SpotParameters{ifish}.Centroid], 2, numel(SpotParameters{ifish}))';
            [in, ~] = inpolygon(rc(:, 1), rc(:, 2), X3, Y3);
            
            SpotsDetec = SpotParameters{ifish}(in' == 1 & ...
                [SpotParameters{ifish}.LargerThan] == 1 & ...
                [SpotParameters{ifish}.SmallerThan] == 1 & ...
                [SpotParameters{ifish}.MinProbability] == 1);
            
            [rc] = reshape([SpotsDetec.Centroid], 2, numel(SpotsDetec))';
            
            % set the new spots
            sc.XData = rc(:, 1);
            sc.YData = rc(:, 2);
            ah.String = sprintf(['Fish: %.0f of %.0f', '\n', ...
                'Computed Spots: %.0f'], ifish, nfishes, numel(sc.XData));
        end
        updatecheckup
    end
    function ResetButton(~, ~)
        %reset button
        ph.XData = [];
        ph.YData = [];
        ln.XData = checkupOrg(ifish).Midbrain(:, 1);
        ln.YData = checkupOrg(ifish).Midbrain(:, 2);
        sc.XData = checkupOrg(ifish).Spots(:, 1);
        sc.YData = checkupOrg(ifish).Spots(:, 2);
        cbh.Value = true;
        ah.String = sprintf(['Fish: %.0f of %.0f', '\n', ...
            'Computed Spots: %.0f'], ifish, nfishes, numel(sc.XData));
        updatecheckup
        
    end
    function SliceButton(source, ~)
        %show slice button
        if source.Value
            sld.Visible = 'off';
            sld2.Visible = 'on';
            uicontrol(sld2);
        else
            sld.Visible = 'on';
            sld2.Visible = 'off';
            uicontrol(sld);
        end
        IniFish
        
    end
    function SliceButton2(source, ~)
        %show slice button
        source.Value
        if source.Value
            ln2.Visible = 'on';
            sc2.Visible = 'on';
        else
            ln2.Visible = 'off';
            sc2.Visible = 'off';
        end
        IniFish
        
    end
    function FishSlider(source, ~)
        %fish slider
        ifish = round(source.Value);
        islice = slice_start(ifish);
        subslice = islice - slice_start(ifish) + 1;
        sld2.Value = islice;
        IniFish;
        
    end
    function SliceSlider(source, ~)
        islice = round(source.Value);
        ifish = find(islice >= slice_start, 1, 'last');
        subslice = islice - slice_start(ifish) + 1;
        sld.Value = ifish;
        IniFish;
    end
    function checkBoxCallback(~, ~)
        %checkbox include
        if cbh.Value
            ln.Visible = 'on';
            sc.Visible = 'on';
            ph.Visible = 'on';
        else
            ln.Visible = 'off';
            sc.Visible = 'off';
            ph.Visible = 'off';
        end
        checkup(ifish).Include = cbh.Value;
        if rad.Value
            uicontrol(sld2);
        else
            uicontrol(sld);
        end
    end
    function SaveButton(~, ~)
        %save button
        updatecheckup
        save([obj.SavePath, '/', obj.InfoName, '.mat'], 'checkup', '-append');
        obj.buildsheet
    end
    function updatecheckup
        %update checkup variable
        checkup(ifish).Midbrain = [ln.XData', ln.YData'];
        checkup(ifish).Spots = [sc.XData', sc.YData'];
        checkup(ifish).Corrections = [ph.XData', ph.YData'];
        checkup(ifish).Include = cbh.Value;
        
        if rad.Value
            uicontrol(sld2);
        else
            uicontrol(sld);
        end
        
    end
    function IniFish
        
        if rad.Value
            tform_complete = RegistrationInfo{ifish}(strcmp({RegistrationInfo{ifish}.name}, 'tform_complete')).value;
            t = Tiff([obj.SavePath, '/', 'CorrectedFish', '/', obj.StackInfo(ifish).stackname, '.tif'], 'r');
            setDirectory(t, subslice)
            CorrectedFish = t.readRGBAImage();
            CorrectedFish = im2uint8(CorrectedFish);
            AlignedSlice = imwarp(CorrectedFish, tform_complete, 'FillValues', 255, 'OutputView', obj.CompleteTemplate.ref_temp);
            sh.CData = uint8(AlignedSlice);
            ah.String = sprintf(['Fish: %.0f of %.0f   Slice: %.f of %.f' , '\n', ...
                'Computed Spots: %.0f'], ifish, nfishes,subslice,obj.StackInfo(ifish).stacksize,numel(sc.XData));
        else
            AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(ifish).stackname, '.tif']);
            sh.CData = uint8(AlignedFish);
            ah.String = sprintf(['Fish: %.0f of %.0f', '\n', ...
                'Computed Spots: %.0f'], ifish, nfishes, numel(sc.XData));
        end
        
        %show by hand annotated midbrain and spots if exist
        if SpotAnn
            sc2.XData = obj.Annotations(ifish).Spots(:,1);
            sc2.YData = obj.Annotations(ifish).Spots(:,2);
        end
        if MidBrainAnn
            ln2 = obj.Annotations(ifish).MidBrain(:,1);
            ln2 = obj.Annotations(ifish).MidBrain(:,2);
        end
        
        %update figures
        ln.XData = checkup(ifish).Midbrain(:, 1);
        ln.YData = checkup(ifish).Midbrain(:, 2);
        sc.XData = checkup(ifish).Spots(:, 1);
        sc.YData = checkup(ifish).Spots(:, 2);
        
        if ~isempty(checkup(ifish).Corrections)
            ph.XData = checkup(ifish).Corrections(:, 1);
            ph.YData = checkup(ifish).Corrections(:, 2);
        else
            ph.XData = [];
            ph.YData = [];
        end
        
        cbh.Value = checkup(ifish).Include;       
        checkBoxCallback
    end
end