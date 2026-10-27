function [selectedMEPs, selectedIdx] = selectingMEP(allMEP, t, BrainsightErrors)

%{
      data should be a matrix of the MEP wished to be analysed with the
      following format :
               each column is a different MEP (MEP(i,:) - EMG data of ith MEP)
%}

nMEP = size(allMEP, 2);
if nargin < 3 || isempty(BrainsightErrors)
    BrainsightErrors = struct( ...
        'TargetError_mm',   nan(nMEP,1), ...
        'AngularError_deg', nan(nMEP,1), ...
        'TwistError_deg',   nan(nMEP,1));
end

% Only alter the layout when at least one MEP has a Brainsight match
hasBrainsight = any(~isnan(BrainsightErrors.TargetError_mm) ...
    | ~isnan(BrainsightErrors.AngularError_deg) ...
    | ~isnan(BrainsightErrors.TwistError_deg));

if hasBrainsight
    figWidth = 1090; panelWidth = 330; cbWidth = 75;
else
    figWidth = 1000; panelWidth = 240; cbWidth = 120;
end

% Create a figure
f = uifigure('Name', 'MEP Selection', 'Position', [100 100 figWidth 600]);

uilabel(f, ...
    'Text', 'Select the MEPs and click "Export Selected MEPs"', ...
    'Position', [50 565 600 25], ...
    'FontSize', 14, ...
    'FontWeight', 'bold');

uilabel(f, ...
    'Text', 'Please make sure to select only MEP-shaped signals to ensure proper use of the next functions.', ...
    'Position', [50 540 700 20], ...
    'FontSize', 13);

fitYCheckbox = uicheckbox(f, ...
    'Text', 'Fit Y axis on MEP window (-100 to 100 ms)', ...
    'Value', false, ...
    'Position', [50 515 400 20]);

ax = uiaxes('Parent', f, 'Position', [50 130 650 380]);
% position in [%]
hold(ax, 'on');

% Plot all the MEPs
hLines=plot(ax,t',allMEP);
xlabel(ax,'Time (ms)')
ylabel(ax,'Amplitude (V)')
xline(ax,0,'r--','Stimulation','LabelHorizontalAlignment','left');
title(ax, 'MEP Selection');

fitYCheckbox.ValueChangedFcn = @(src,evt)toggleYFitWindow(src,f,ax,hLines);

% Create checkbox panel (empty)
panel = uipanel(f,...
    'Title','MEPs',...
    'Position',[730 50 panelWidth 490],...
    'Scrollable','on');

if hasBrainsight
    legendY = nMEP*26 + 10;
    uilabel(panel, ...
        'Text','P=Position (mm)   T=Tangence (°)   O=Orientation (°)', ...
        'Position',[15 legendY 300 16], ...
        'FontSize',10, 'FontColor',[0.4 0.4 0.4]);
end

% Add all the checkboxes and their state
f.UserData.cb = cell(nMEP,1);
f.UserData.completed = false;
f.UserData.selectedMEPs = [];
f.UserData.selectedIdx = [];
f.UserData.lastSelected = nMEP;
f.UserData.displayMode = "all";
f.UserData.fitYWindow = false;
f.UserData.t = t;
f.UserData.yFitRange = [-100 100];

for i = 1:nMEP
    rowY = nMEP*26-25*i;
    f.UserData.cb{i} = uicheckbox(panel,...
        'Text',sprintf('MEP %d',i),...
        'Value',true,...
        'Position',[15 rowY cbWidth 20],...
        'ValueChangedFcn',...
        @(src,evt)toggleMEP(src,i,f,ax,hLines));

    if hasBrainsight
        err = struct( ...
            'TargetError_mm',   BrainsightErrors.TargetError_mm(i), ...
            'AngularError_deg', BrainsightErrors.AngularError_deg(i), ...
            'TwistError_deg',   BrainsightErrors.TwistError_deg(i));
        uilabel(panel, ...
            'Text', formatBrainsightHTML(err), ...
            'Interpreter','html', ...
            'Position',[95 rowY 215 20]);
    end
end

% Buttons

uibutton(f,...
    'Text','Select All',...
    'Position',[50 50 110 40],...
    'ButtonPushedFcn',...
    @(src,evt)selectAll(f,ax,hLines));

uibutton(f,...
    'Text','Deselect All',...
    'Position',[170 50 110 40],...
    'ButtonPushedFcn',...
    @(src,evt)deselectAll(f,ax,hLines));

uilabel(f,...
    'Text','Display Mode:',...
    'Position',[300 60 80 22],...
    'HorizontalAlignment', 'right');

uidropdown(f,...
    'Items',{'Show All Selected','Show Last Selected Only'},...
    'Value','Show All Selected',...
    'Position',[390 60 150 22],...
    'ValueChangedFcn',...
    @(src,evt)changeDisplayMode(src,f,ax,hLines));

uibutton(f, 'Text', 'Export Selected MEPs', ...
    'Position', [560 50 140 40], ...
    'BackgroundColor', [0.8 0.95 0.8], ...
    'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(src,evt)extractingSelectedMEPs(allMEP,f));

% Wait for the user to complete selection

updateDisplay(f,hLines,ax);
uiwait(f)

% Retrieve results
if isvalid(f)
    selectedMEPs = f.UserData.selectedMEPs;
    selectedIdx  = f.UserData.selectedIdx;
    delete(f);
else
    selectedMEPs = [];
    selectedIdx = [];
end

end

%% Function that will display or not MEP

function toggleMEP(src,idx,f,ax,hLines)
if src.Value
    f.UserData.lastSelected = idx;
else
    if f.UserData.lastSelected == idx
        f.UserData.lastSelected = [];
        for k = length(f.UserData.cb):-1:1
            if f.UserData.cb{k}.Value
                f.UserData.lastSelected = k;
                break;
            end
        end
    end
end
updateDisplay(f,hLines,ax);
end

%% Functions to change the display mode

function changeDisplayMode(src,f,ax,hLines)
if strcmp(src.Value,'Show All Selected')
    f.UserData.displayMode = "all";
else
    f.UserData.displayMode = "last";
end
updateDisplay(f,hLines,ax);
end

%% Function to toggle fitting the Y axis to the MEP window only

function toggleYFitWindow(src,f,ax,hLines)
f.UserData.fitYWindow = logical(src.Value);
updateDisplay(f,hLines,ax);
end

function updateDisplay(f,hLines,ax)
cb = f.UserData.cb;
selected = cellfun(@(x) logical(x.Value), cb, 'UniformOutput', true);

switch f.UserData.displayMode
    case "all"
        for k = 1:length(hLines)
            if selected(k)
                hLines(k).Visible = 'on';
            else
                hLines(k).Visible = 'off';
            end
        end

    case "last"
        for k = 1:length(hLines)
            hLines(k).Visible = 'off';
        end

        idx = f.UserData.lastSelected;
        if ~isempty(idx) && idx <= length(hLines)
            if selected(idx)
                hLines(idx).Visible = 'on';
            end
        end
end

applyYFit(f,ax,hLines);
drawnow
end

%% Function to set the Y axis limits, either auto or fit to the MEP window

function applyYFit(f,ax,hLines)
if isempty(ax) || ~isvalid(ax)
    return
end
if ~f.UserData.fitYWindow
    ax.YLimMode = 'auto';
    return
end

t = f.UserData.t;
winMask = t >= f.UserData.yFitRange(1) & t <= f.UserData.yFitRange(2);
isVisible = strcmp({hLines.Visible}, 'on');

windowData = [];
for k = find(isVisible)
    ydata = hLines(k).YData;
    windowData = [windowData, ydata(winMask)]; %#ok<AGROW>
end

if isempty(windowData)
    ax.YLimMode = 'auto';
    return
end

yMin = min(windowData);
yMax = max(windowData);
if yMax == yMin
    pad = max(abs(yMin), 1) * 0.1;
else
    pad = (yMax - yMin) * 0.05;
end
ax.YLim = [yMin - pad, yMax + pad];
end

%% Functions to select and deselect all MEPs

function selectAll(f,ax,hLines)
for k = 1:length(f.UserData.cb)
    f.UserData.cb{k}.Value = true;
end
f.UserData.lastSelected = length(f.UserData.cb);
updateDisplay(f,hLines,ax);
end

function deselectAll(f,ax,hLines)
for k = 1:length(f.UserData.cb)
    f.UserData.cb{k}.Value = false;
end
f.UserData.lastSelected = [];
updateDisplay(f,hLines,ax);
end

%% Function that returns only the selected MEPs

function extractingSelectedMEPs(allMEP,f)
cb = f.UserData.cb;
selected = cellfun(@(x) logical(x.Value), cb, 'UniformOutput', true);

% if only want the MEP signals:
% % Extract only selected MEPs
% selectedMEPs = MEP(:, selected);   % rows = trials, cols = time points
%
% % Transpose so that each column is one trial
% resultMatrix = selectedMEPs';      % now: (time points × trials)
%
% % Display size in command window
% disp(size(resultMatrix));
% % Export to workspace
% assignin('base', 'SelectedMEPs', resultMatrix);

% Collect all the samples of selected MEPs
f.UserData.selectedIdx = find(selected);
f.UserData.selectedMEPs = allMEP(:,selected);
f.UserData.completed = true;

% Export to workspace
% assignin('base', 'SelectedMEPs', selectedMEPs);
% assignin('base', 'SelectedIdx', selectedIdx);

% Display confirmation in command window
fprintf('Exported %d selected MEPs.\n', sum(selected));
uiresume(f)
end

%% Brainsight error display helpers

function html = formatBrainsightHTML(err)
% Builds the HTML snippet shown next to a MEP checkbox: Position (P, mm),
% Tangence (T, deg) and Orientation (O, deg) errors.
if isnan(err.TargetError_mm) && isnan(err.AngularError_deg) && isnan(err.TwistError_deg)
    html = '<html><font color="#808080">&mdash;</font></html>';
    return
end
cP = errorColorHex(err.TargetError_mm,   3, 5,  false);
cT = errorColorHex(err.AngularError_deg, 10,15, false);
cO = errorColorHex(err.TwistError_deg,   10,15, true);
html = sprintf(['<html>P:<font color="%s">%.1fmm</font>&nbsp;' ...
                 'T:<font color="%s">%.1f&deg;</font>&nbsp;' ...
                 'O:<font color="%s">%.1f&deg;</font></html>'], ...
    cP, err.TargetError_mm, cT, err.AngularError_deg, cO, err.TwistError_deg);
end

function hex = errorColorHex(val, greenMax, redMin, signed)
% TargetError <=3/5mm, AngularError/TwistError <=10/15deg
if isnan(val)
    hex = '#808080';
    return
end
if signed
    v = abs(val);
else
    v = val;
end
if v <= greenMax
    hex = '#458556';
elseif v <= redMin
    hex = '#BF9000';
else
    hex = '#BA4F50';
end
end