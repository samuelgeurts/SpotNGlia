function objout = CheckFish(obj, ifish, INFO)
%application to adjust brain region and spots

%used for assignin to the righ workspace variable
%objectname = inputname(1);

global fh sh ln ln2 sc sc2 str ah ph ph2 ph3 rec %figure and axis handles
global uih btn btn2 btn3 btn4 btn5 rad rad2 rad3 sld sld2 cbh cbh2 %uicontrol handles


if ~exist('ifish', 'var') || (ifish > numel(obj.StackInfo))
    ifish = 1;
end

h = waitbar(0, 'Load Template', 'Name', 'Loading Data');

%obj = LoadTemplate(obj);

waitbar(0.03, h, 'Brain, Spot, Registration')

if ~exist('INFO', 'var') || ...
        ~isfield(INFO, 'BrainSegmentationInfo') || ...
        ~isfield(INFO, 'SpotParameters') || ...
        ~isfield(INFO, 'RegistrationInfo')
    INFO = load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo', 'SpotParameters', 'RegistrationInfo');
end

% waitbar(0.03, h, 'BrainSegmentationInfo')
% load([obj.SavePath, '/', obj.InfoName, '.mat'], 'BrainSegmentationInfo')
% waitbar(0.17, h, 'SpotParameters')
% load([obj.SavePath, '/', obj.InfoName, '.mat'], 'SpotParameters');
% waitbar(0.87, h, 'RegistrationInfo')
% load([obj.SavePath, '/', obj.InfoName, '.mat'], 'RegistrationInfo');
%


waitbar(0.93, h, 'Initialize')

slice_end = cumsum([obj.StackInfo.stacksize]);
slice_start = [1, slice_end(1:end-1) + 1];
islice = slice_start(ifish);
subslice = 1;

cxy = fliplr(obj.CompleteTemplate.CenterMidBrain);
nfishes = numel(obj.StackInfo);
nslices = sum([obj.StackInfo.stacksize]);

%check if the field exist and if at least one field is not empty
SpotAnn = isfield(obj.Annotations, 'Spots') && any(~cellfun(@isempty, {obj.Annotations.Spots}));
MidBrainAnn = isfield(obj.Annotations, 'MidBrain') && any(~cellfun(@isempty, {obj.Annotations.MidBrain}));

InitializeCheckup
InitializeImageFigure
InitializeButtonFigure

waitbar(1, h, 'Initialize')

IniFish
uicontrol(sld);
delete(h)
uiwait(uih)

    function InitializeCheckup
        FillComputations(obj, INFO.BrainSegmentationInfo);
        FillCheckup(obj);
        
        if isempty(obj.checkup)
            obj.checkup = obj.Computations;
            obj.checkup(1).Corrections = [];
            obj.checkup(1).SpotAdditions = [];
            obj.checkup(1).SpotRemovals = [];
            [obj.checkup.Include] = deal(true);
            [obj.checkup.Polar] = deal(true);
        end
        
        %if an older checkup file is used SpotNGlia1.4.0 and before, the field 'SpotAdditions' and 'SpotRemovals' are added
        if ~isfield(obj.checkup, 'Counts')
            for k1 = 1:numel(obj.checkup)
                obj.checkup(k1).Counts = size(obj.checkup(k1).Spots, 1);
            end
        end
        if ~isfield(obj.checkup, 'SpotAdditions')
            obj.checkup(1).SpotAdditions = [];
        end
        if ~isfield(obj.checkup, 'SpotRemovals')
            obj.checkup(1).SpotRemovals = [];
        end
        if ~isfield(obj.checkup, 'Polar')
            [obj.checkup.Polar] = deal(true);
        end
    end
    function InitializeImageFigure
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
            'LineWidth', 2, ...
            'Visible', 'off');
        sc = scatter(0, 0, 400, ...
            'LineWidth', 2, ...
            'MarkerEdgeColor', 1/255*[255, 75, 75]);
        sc2 = scatter(0, 0, 400, ...
            'LineWidth', 2, ...
            'MarkerEdgeColor', 1/255*[75, 75, 255], ...
            'SizeData', 300, ...
            'Visible', 'off');
        str = sprintf(['Fish: %.0f of %.0f', '\n', ...
            'Computed Spots: %.0f'], 0, 0, 0);
        ah = annotation('textbox', [.05, .63, .7, .3], ...
            'String', str, ...
            'FitBoxToText', 'on', 'FontSize', 15);
        ph = plot(0, 0, ...
            'bo');
        ph2 = plot(0, 0, ...
            'bx');
        ph3 = plot(0, 0, ...
            'rx');
        rec = rectangle('Position', [cxy(2) - 500, cxy(1) - 500, 1000, 1000], ...
            'Visible', 'off');
    end
    function InitializeButtonFigure
        %createsa figure below the place first figure is placed with buttons
        %this way, the top figure can easyly changed size without change button location
        
        sizb = 70; %vertical size of sub figure below
        barcor = 22; %distance from top figure, 22 is zero distance including clos bar
        
        posf = fh.Position;
        posb = [posf(1), posf(2) - sizb - barcor, posf(3), sizb];
        
        uifignumber = fh.Number + 1;
        uih = figure(uifignumber);
        set(gcf, 'position', posb);
        set(gcf, 'Toolbar', 'none');
        set(gcf, 'Menubar', 'none');
        
        
        btn = uicontrol('Style', 'pushbutton', 'String', 'Correct Brain', ...
            'Position', [30, 10, 100, 22], ...
            'Callback', {@CorrectBrainButton});
        
        btn2 = uicontrol('Style', 'pushbutton', 'String', 'Brain Reset', ...
            'Position', [130, 10, 100, 22], ...
            'Callback', {@ResetButton});
        
        btn3 = uicontrol('Style', 'pushbutton', 'String', 'Correct Spot', ...
            'Position', [340, 10, 100, 22], ...
            'Callback', {@CorrectSpotButton});
        
        btn4 = uicontrol('Style', 'pushbutton', 'String', 'Spot Reset', ...
            'Position', [440, 10, 100, 22], ...
            'Callback', {@ResetButton2});
        
        rad = uicontrol('Style', 'radiobutton', 'String', 'Show Slices', ...
            'Position', [340, 41, 85, 22], ...
            'Value', 0, ...
            'Callback', {@RadioSliceButton});
        
        rad2 = uicontrol('Style', 'radiobutton', 'String', 'Compare', ...
            'Position', [425, 41, 75, 22], ...
            'Value', 0, ...
            'Callback', {@RadioAnnotationButton});
        
        rad3 = uicontrol('Style', 'radiobutton', 'String', 'Info', ...
            'Position', [500, 41, 50, 22], ...
            'Value', 1, ...
            'Callback', {@RadioInfoButton});
        
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
            'Value', obj.checkup(ifish).Include, 'Position', [240, 41, 100, 22], ...
            'Callback', @IncludeCheckBox);
        
        cbh2 = uicontrol('Style', 'radiobutton', 'String', 'Polar', ...
            'Value', obj.checkup(ifish).Polar, 'Position', [240, 10, 100, 22], ...
            'Callback', @PolarCheckBox);
        
        btn5 = uicontrol('Style', 'pushbutton', 'String', 'Save', ...
            'Position', [580, 41, 100, 22], ...
            'Callback', {@SaveButton});
        btn6 = uicontrol('Style', 'pushbutton', 'String', 'Save and Close', ...
            'Position', [580, 10, 100, 22], ...
            'Callback', {@CloseButton});%#ok<NASGU>
        
        
    end
    function CorrectBrainButton(~, ~)
        %correct button
        
        rec.Visible = 'on';
        figure(fh)
        
        
        %figure(fh)
        %hold on
        %rectangle('Position', [cxy(2) - 500, cxy(1) - 500, 1000, 1000]);
        
        but = 1;
        xy = obj.checkup(ifish).Corrections;
        n = size(obj.checkup(ifish).Corrections, 1);
        
        % pick input and display
        while but == 1 || (but == 3)
            [x, y, but] = ginput(1); %pick a new input
            
            if but == 1
                n = n + 1;
                xy(n, :) = [x, y];
                ph.XData = xy(:, 1);
                ph.YData = xy(:, 2);
            end
            if but == 3
                
                [mn, ind] = min(sqrt((ph.XData - x).^2+(ph.YData - y).^2));
                
                if mn <= 20
                    xy(ind, :) = [];
                    ph.XData(ind) = [];
                    ph.YData(ind) = [];
                    n = n - 1;
                end
            end
            if (cbh2.Value == 0) && (size(xy, 1) >= 3)
                
                [X3, Y3] = ComputeBrain(xy);
                rc = obj.SpotsInsiteArea(INFO.SpotParameters{ifish}, [X3, Y3]);
                
                ln.XData = X3;
                ln.YData = Y3;
                ln.Marker = 'none';
                sc.XData = rc(:, 1);
                sc.YData = rc(:, 2);
                
                updatespots
                updateannotation
            end
        end
        
        [X3, Y3] = ComputeBrain(xy);
        rc = obj.SpotsInsiteArea(INFO.SpotParameters{ifish}, [X3, Y3]);
        
        
        ln.XData = X3;
        ln.YData = Y3;
        ln.Marker = 'none';
        sc.XData = rc(:, 1);
        sc.YData = rc(:, 2);
        
        updatespots
        updateannotation
        updatecheckup
        
        rec.Visible = 'off';
        
    end
    function CorrectSpotButton(~, ~)
        %left click adds a spot, right clic removes a spot
        %other button quits function
        
        rec.Visible = 'on';
        figure(fh)
        but = 1;
        
        %pick input and display
        while (but == 1) || (but == 3)
            [x, y, but] = ginput(1); %pick a new input
            
            %left clic to add
            if but == 1
                [mn, ~] = min(sqrt((sc.XData - x).^2+(sc.YData - y).^2));
                if mn >= 10
                    sc.XData = [sc.XData, x];
                    sc.YData = [sc.YData, y];
                    ph2.XData = [ph2.XData, x];
                    ph2.YData = [ph2.YData, y];
                    
                else
                    warning('new spot is very close to a previous annotated spot')
                    beep
                end
            end
            
            %right clic to remove button
            if but == 3
                [mn, mi] = min(sqrt((sc.XData - x).^2+(sc.YData - y).^2));
                if mn <= 20
                    
                    %checks if removed spot is from a previous added spot
                    %if so, spot is removed from added spot list
                    %if not, spot is added to removed spot list
                    [TF, ind] = ismember([sc.XData(mi), sc.YData(mi)], [ph2.XData', ph2.YData'], 'rows');
                    if TF
                        ph2.XData(ind) = [];
                        ph2.YData(ind) = [];
                    else
                        ph3.XData = [ph3.XData, sc.XData(mi)];
                        ph3.YData = [ph3.YData, sc.YData(mi)];
                        
                    end
                    %remove spot from scatterplot
                    sc.XData(mi) = [];
                    sc.YData(mi) = [];
                end
            end
        end
        updateannotation
        updatecheckup
    end
    function ResetButton(~, ~)
        %reset button which resets brain but keeps added and removed spots
        ph.XData = [];
        ph.YData = [];
        ln.XData = obj.Computations(ifish).Midbrain(:, 1);
        ln.YData = obj.Computations(ifish).Midbrain(:, 2);
        sc.XData = obj.Computations(ifish).Spots(:, 1);
        sc.YData = obj.Computations(ifish).Spots(:, 2);
        cbh.Value = true;
        cbh2.Value = true;
        
        updatespots
        updateannotation
        updatecheckup
        
    end
    function ResetButton2(~, ~)
        %reset button which only removes added spots and adds removed spots
        
        %         %remove previous added spots
        %         [TF, ~] = ismember([sc.XData', sc.YData'], [ph2.XData', ph2.YData'], 'rows');
        %         sc.XData(TF) = [];
        %         sc.YData(TF) = [];
        %         %add previous removed spots
        %         sc.XData = [sc.XData, ph3.XData];
        %         sc.YData = [sc.YData, ph3.YData];
        
        
        %recompute spots from corrected brain
        rc = obj.SpotsInsiteArea(INFO.SpotParameters{ifish}, obj.checkup(ifish).Midbrain);
        sc.XData = rc(:, 1);
        sc.YData = rc(:, 2);
        %set added and removed spots to zero
        ph2.XData = [];
        ph2.YData = [];
        ph3.XData = [];
        ph3.YData = [];
        
        updateannotation
        updatecheckup
    end
    function RadioSliceButton(~, ~)
        %show slice button
        if rad.Value
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
    function RadioAnnotationButton(~, ~)
        %show slice button
        if rad2.Value
            ln2.Visible = 'on';
            sc2.Visible = 'on';
        else
            ln2.Visible = 'off';
            sc2.Visible = 'off';
        end
        %IniFish
        
    end
    function RadioInfoButton(~, ~)
        if rad3.Value
            ah.Visible = 'on';
        else
            ah.Visible = 'off';
        end
    end
    function FishSlider(~, ~)
        %fish slider
        ifish = round(sld.Value);
        islice = slice_start(ifish);
        subslice = islice - slice_start(ifish) + 1;
        sld2.Value = islice;
        IniFish;
        
    end
    function SliceSlider(~, ~)
        islice = round(sld2.Value);
        ifish = find(islice >= slice_start, 1, 'last');
        subslice = islice - slice_start(ifish) + 1;
        sld.Value = ifish;
        IniFish;
    end
    function SaveButton(~, ~)
        %save button
        h = waitbar(0, 'update checkup', 'Name', 'Saving');
        updatecheckup
        figure(h)
        waitbar(1/50, h, 'save checkup')
        
        %save([obj.SavePath, '/', obj.InfoName, '.mat'], 'checkup', '-append');
        obj.saveit
        
        waitbar(4/6, h, 'save excel sheet')
        obj = obj.buildsheet;
        waitbar(1, h, 'complete')
        delete(h)
        
        %update the object in workspace, close and reopen CheckFish
        objout = obj;
        close(uih);
        close(fh);
        objout = CheckFish(obj, ifish, INFO);
        
        % a bit ugly way to update the object into the workspace
        %assignin('base',objectname,obj) % Assign the data to the base workspace.
    end
    function CloseButton(~, ~)
        h = waitbar(0, 'update checkup', 'Name', 'Saving');
        updatecheckup
        figure(h)
        waitbar(1/50, h, 'save checkup')
        
        %save([obj.SavePath, '/', obj.InfoName, '.mat'], 'checkup', '-append');
        obj.saveit
        
        waitbar(4/6, h, 'save excel sheet')
        obj.buildsheet;
        waitbar(1, h, 'complete')
        delete(h)
        
        %update the object in workspace, close and reopen CheckFish
        %objout = obj; %not nesseary anymore as obj (SpotNGlia) is a handle class now
        close(uih);
        close(fh);
    end
    function IncludeCheckBox(~, ~)
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
        obj.checkup(ifish).Include = cbh.Value;
        if rad.Value
            uicontrol(sld2);
        else
            uicontrol(sld);
        end
    end
    function PolarCheckBox(~, ~)
        %checkbox Polar
        obj.checkup(ifish).Polar = cbh2.Value;
        
        if size(obj.checkup(ifish).Corrections, 1) <= 2
            warning('at least 3 brain corrections has to me made to select area')
            beep
            %cbh2.Value = 1;
        else
            [X3, Y3] = ComputeBrain([ph.XData', ph.YData']);
            %[X3,Y3] = ComputeBrain(obj.checkup(ifish).Corrections);
            rc = obj.SpotsInsiteArea(INFO.SpotParameters{ifish}, [X3, Y3]);
            
            ln.XData = X3;
            ln.YData = Y3;
            ln.Marker = 'none';
            sc.XData = rc(:, 1);
            sc.YData = rc(:, 2);
            
            updatespots
            updateannotation
            updatecheckup
        end
        
        if rad.Value
            uicontrol(sld2);
        else
            uicontrol(sld);
        end
    end
    function updatespots
        %this function updates the spots according to previous added or removed spots
        
        %remove previous removed spots
        [TF, ~] = ismember([sc.XData', sc.YData'], [ph3.XData', ph3.YData'], 'rows');
        sc.XData(TF) = [];
        sc.YData(TF) = [];
        %add previous added spots
        sc.XData = [sc.XData, ph2.XData];
        sc.YData = [sc.YData, ph2.YData];
    end
    function updatecheckup
        %update checkup variable
        obj.checkup(ifish).Midbrain = [ln.XData', ln.YData'];
        obj.checkup(ifish).Spots = [sc.XData', sc.YData'];
        obj.checkup(ifish).Corrections = [ph.XData', ph.YData'];
        obj.checkup(ifish).Include = cbh.Value;
        obj.checkup(ifish).Polar = cbh2.Value;
        obj.checkup(ifish).SpotAdditions = [ph2.XData', ph2.YData'];
        obj.checkup(ifish).SpotRemovals = [ph3.XData', ph3.YData'];
        obj.checkup(ifish).Counts = numel(sc.XData);
        
        if rad.Value
            uicontrol(sld2);
        else
            uicontrol(sld);
        end
        
    end
    function updateannotation
        if rad.Value
            ah.String = sprintf(['Fish: %.0f of %.0f   Slice: %.f of %.f', '\n', ...
                'Computed Spots: %.0f'], ifish, nfishes, subslice, obj.StackInfo(ifish).stacksize, numel(sc.XData));
        else
            ah.String = sprintf(['Fish: %.0f of %.0f', '\n', ...
                'Computed Spots: %.0f'], ifish, nfishes, numel(sc.XData));
        end
    end
    function IniFish
        
        if rad.Value
            tform_complete = INFO.RegistrationInfo{ifish}(strcmp({INFO.RegistrationInfo{ifish}.name}, 'tform_complete')).value;
            t = Tiff([obj.SavePath, '/', 'CorrectedFish', '/', obj.StackInfo(ifish).stackname, '.tif'], 'r');
            setDirectory(t, subslice)
            CorrectedFish = t.readRGBAImage();
            CorrectedFish = im2uint8(CorrectedFish);
            AlignedSlice = imwarp(CorrectedFish, tform_complete, 'FillValues', 255, 'OutputView', obj.CompleteTemplate.ref_temp);
            sh.CData = uint8(AlignedSlice);
        else
            AlignedFish = imread([obj.SavePath, '/', 'AlignedFish', '/', obj.StackInfo(ifish).stackname, '.tif']);
            sh.CData = uint8(AlignedFish);
            %ah.String = sprintf(['Fish: %.0f of %.0f', '\n', ...
            %    'Computed Spots: %.0f'], ifish, nfishes, numel(sc.XData));
        end
        
        
        %show by hand annotated midbrain and spots if exist
        if SpotAnn
            if ~isnan(obj.Annotations(ifish).Spots(1))
                sc2.XData = obj.Annotations(ifish).Spots(:, 1);
                sc2.YData = obj.Annotations(ifish).Spots(:, 2);
            else
                sc2.XData = [];
                sc2.YData = [];
            end
        end
        if MidBrainAnn
            if ~isnan(obj.Annotations(ifish).MidBrain(1))
                ln2.XData = obj.Annotations(ifish).MidBrain(:, 1);
                ln2.YData = obj.Annotations(ifish).MidBrain(:, 2);
            else
                ln2.XData = [];
                ln2.YData = [];
            end
        end
        
        %update figures
        ln.XData = obj.checkup(ifish).Midbrain(:, 1);
        ln.YData = obj.checkup(ifish).Midbrain(:, 2);
        
        
        %update scatter data: spots,braincorrections,addition,removals
        if ~isempty(obj.checkup(ifish).Spots)
            sc.XData = obj.checkup(ifish).Spots(:, 1);
            sc.YData = obj.checkup(ifish).Spots(:, 2);
        else
            sc.XData = [];
            sc.YData = [];
        end
        if ~isempty(obj.checkup(ifish).Corrections)
            ph.XData = obj.checkup(ifish).Corrections(:, 1);
            ph.YData = obj.checkup(ifish).Corrections(:, 2);
        else
            ph.XData = [];
            ph.YData = [];
        end
        if ~isempty(obj.checkup(ifish).SpotAdditions)
            ph2.XData = obj.checkup(ifish).SpotAdditions(:, 1);
            ph2.YData = obj.checkup(ifish).SpotAdditions(:, 2);
        else
            ph2.XData = [];
            ph2.YData = [];
        end
        if ~isempty(obj.checkup(ifish).SpotRemovals)
            ph3.XData = obj.checkup(ifish).SpotRemovals(:, 1);
            ph3.YData = obj.checkup(ifish).SpotRemovals(:, 2);
        else
            ph3.XData = [];
            ph3.YData = [];
        end
        
        updateannotation
        cbh.Value = obj.checkup(ifish).Include;
        IncludeCheckBox
    end
    function [X3, Y3] = ComputeBrain(xy)
        
        PolarN = INFO.BrainSegmentationInfo(ifish).PolarTransform;
        sp = fliplr(size(PolarN));
        si = size(PolarN);
        
        if ~isempty(xy)
            if cbh2.Value
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
                X3 = X3';
                Y3 = Y3';
                
            else
                %sorting method for coordinates to polinome
                %[xy2,~] = sng_OrderContourCoordinates(xy);
                
                %other sorting method with known center
                xy2 = xy - mean(xy, 1);
                [~, ind] = sort(cart2pol(xy2(:, 1), xy2(:, 2)));
                
                X3 = xy(ind, 1);
                Y3 = xy(ind, 2);
                %fix endpoints
                X3 = [X3; X3(1)];
                Y3 = [Y3; Y3(1)];
                
            end
        end
    end


end