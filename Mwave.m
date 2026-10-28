%% M-WAVE P2P TOOL — LOCAL PNS DIAGNOSTICS / GUI LAYOUT REFINEMENT
% 1) Load one CFS or MAT file
% 2) Display raw EMG signal in ms scale (to make it apply filter change the line 133 to true)
% 3) Start with a small ms window for better visualization
% 4) Scroll through signal using horizontal slider or mouse wheel
% 5) Zoom/pan manually if needed
% 6) Press "Select M-wave", then click START and END
% 7) Calculates P2P and AUC for each selected window
% 8) Chronological Undo Last Edit for manual selections and markers
% 9) Saves CSV only on Done; closing X cancels without saving
% 10) PNS marker starts at ADC edge 
clc; clear; close all
addpath("functions");
addpath("matcfs64c");
%% ================= USER SETTINGS =================
defaultDir   = pwd;  % Start file selection in the current working folder
EMG_field = '';  % Select from available channels when the recording is loaded
% Display settings
initialWindow_ms = 1000;   % width of the initial view in ms
scrollStep_ms    = 500;    % mouse wheel movement in ms
% Save folder
saveFolder = pwd;  % Save outputs in the current working folder
%% ================= PICK FILE =================
[fileName, fileDir] = uigetfile( ...
    {'*.cfs;*.mat','CFS or MAT (*.cfs, *.mat)'; '*.cfs','CFS (*.cfs)'; '*.mat','MAT (*.mat)'}, ...
    'Select file containing M-waves', defaultDir);
if isequal(fileName,0)
    error('No file selected.');
end
fullPath = fullfile(fileDir, fileName);
%% ================= LOAD ORIGINAL EMG (NO FILTER) =================
[t, x, Fs, ~, tADC, adcSignal, EMG_field] = loadOneEMG(fullPath, EMG_field);
fprintf('\nLoaded file:\n%s\n', fullPath);
fprintf('EMG channel: %s\n', EMG_field);
fprintf('Sampling frequency: %.1f Hz\n', Fs);
fprintf('EMG filtering: disabled (original signal).\n\n');
%% ================= RUN TOOL =================
[Results, PeriodResults, didFinish] = runFullSignalMwaveTool(t, x, Fs, fileName, EMG_field, initialWindow_ms, scrollStep_ms, tADC, adcSignal);
if ~didFinish
    fprintf('Session cancelled: no results saved.\n');
    return;
end
%% ================= ONE CLEAN EXPORT =================
% One CSV with one row per manually selected M-wave.
[~, baseName] = fileparts(fileName);
if ~exist(saveFolder, 'dir'), mkdir(saveFolder); end
outCSV = fullfile(saveFolder, sprintf('%s_%s_Mwave_results.csv', baseName, EMG_field));
% Keep the clinical/experimental measurements and a simple distinction
% between manual and automatic intervals; omit plotting and diagnostic data.
if isempty(Results)
    Export = table('Size',[0 10], ...
        'VariableTypes',{'double','string','double','double','double','double','double','double','double','double'}, ...
        'VariableNames',{'Mwave','Selection','PNS_ms','Onset_ms','Offset_ms', ...
        'Duration_ms','P2P_uV','AbsoluteAUC_uVms','Max_uV','Min_uV'});
else
    n = height(Results);
    selection = repmat("Manual",n,1);
    pns = nan(n,1);
    Export = table(Results.Mwave,selection,pns,Results.tStart_ms, ...
        Results.tEnd_ms,Results.Duration_ms,Results.P2P_uV, ...
        Results.AbsoluteAUC_uVms,Results.Max_uV,Results.Min_uV, ...
        'VariableNames',{'Mwave','Selection','PNS_ms','Onset_ms','Offset_ms', ...
        'Duration_ms','P2P_uV','AbsoluteAUC_uVms','Max_uV','Min_uV'});
end
writetable(Export,outCSV);
fprintf('Saved one M-wave results file:\n%s\n',outCSV);
disp(Export);
%% =======================================================================
%% ========================= LOCAL FUNCTIONS =============================
%% =======================================================================
function [t, x, Fs, usedFilter, tADC, adcSignal, EMG_field] = loadOneEMG(fullPath, EMG_field)
    [~, ~, ext] = fileparts(fullPath);
    ext = lower(ext);
    switch ext
        case '.cfs'
            tmp = readCFSfile(fullPath);
        case '.mat'
            S = load(fullPath);
            if isfield(S,'tmp')
                tmp = S.tmp;
            else
                tmp = S;
            end
        otherwise
            error('Unsupported file type: %s', ext);
    end
    if ~isstruct(tmp)
        error('Loaded file is not a struct.');
    end
    f = fieldnames(tmp);
    if isfield(tmp,'tmp')
        raw = tmp.tmp;
    elseif numel(f)==1
        raw = tmp.(f{1});
    else
        raw = tmp;
    end
    % ADC1_2 is a PNS timing reference, kept separate from all EMG processing.
   tADC = [];
   adcSignal = [];
   if isfield(raw, 'ADC1_2')
       adcChannel = raw.ADC1_2;
       if isfield(adcChannel,'dat') && isfield(adcChannel,'FreqS') && ...
               isnumeric(adcChannel.FreqS) && adcChannel.FreqS > 0
           adcSignal = double(adcChannel.dat(:));
           tADC = (0:numel(adcSignal)-1)' / double(adcChannel.FreqS);
       else
           warning('ADC1_2 found, but dat or FreqS is invalid. ADC trace disabled.');
       end
   else
       warning('ADC1_2 not found. The EMG tool will still work normally.');
   end
   raw_fields = fieldnames(raw);
    data = struct();
    for i = 1:numel(raw_fields)
        oldName = raw_fields{i};
        if startsWith(oldName,'EMG_')
            newName = oldName;
        else
            newName = regexprep(oldName, '_\d+$', '');
        end
        data.(newName) = raw.(oldName);
    end
    % Offer only EMG channels present in this recording.
    emgs = fieldnames(data);
    emgs = emgs(startsWith(emgs,'EMG_'));
    if isempty(emgs)
        error('No EMG channels found in %s.', fullPath);
    end
    if isempty(EMG_field)
        [~, recordingName, recordingExt] = fileparts(fullPath);
        [selection, confirmed] = listdlg( ...
            'PromptString', {'Select the EMG channel to analyze:', ...
                [recordingName recordingExt]}, ...
            'Name', 'Select EMG channel', ...
            'SelectionMode', 'single', ...
            'ListString', emgs, ...
            'ListSize', [320 200]);
        if ~confirmed || isempty(selection)
            error('EMG:SelectionCancelled', 'Analysis cancelled: no EMG channel selected.');
        end
        EMG_field = emgs{selection};
    end

    if ~isfield(data, EMG_field)
        emgs = fieldnames(data);
        emgs = emgs(startsWith(emgs,'EMG_'));
        error('EMG "%s" not found. Available EMG fields: %s', ...
            EMG_field, strjoin(emgs, ', '));
    end
    sig = data.(EMG_field).dat(:);
    Fs  = data.(EMG_field).FreqS;
    t   = (0:numel(sig)-1)' / Fs;
    usedFilter = false;
    x = sig; % Preserve original acquired EMG samples; no digital filter.

end
function [Results, PeriodResults, didFinish] = runFullSignalMwaveTool(t, x, Fs, fileName, EMG_field, initialWindow_ms, scrollStep_ms, tADC, adcSignal)
    t_ms = t(:) * 1000;
    x_uV = -x(:) * 1e6;
    totalDuration_ms = t_ms(end);
    if initialWindow_ms > totalDuration_ms
        initialWindow_ms = totalDuration_ms;
    end
    fig = figure( ...
        'Color','w', ...
        'Name', sprintf('M-wave P2P Tool: %s', fileName), ...
        'NumberTitle','off', ...
        'WindowScrollWheelFcn', @mouseScroll, ...
        'CloseRequestFcn', @cancelTool);
    ax = axes('Parent', fig, ...
        'Position', [0.07 0.52 0.90 0.40]);
    % Use independent left/right amplitude scales; both traces share milliseconds.
   % ADC1_2 is plotted only and never used in P2P calculations.
   yyaxis(ax, 'left');
   hold(ax,'on');
   hEMG = plot(ax, t_ms, x_uV, 'k');
   if ~isempty(adcSignal)
       yyaxis(ax, 'right');
       hADC = plot(ax, tADC(:)*1000, adcSignal(:), ...
           'Color', [0.85 0.2 0.15], 'LineWidth', 1, 'Visible', 'off');
       ylabel(ax, 'PNS trigger (ADC1_2, raw units)');
       ax.YColor = [0.85 0.2 0.15];
       % Initially hide right-side axis labels along with the trigger trace.
       ax.YTickLabel = [];
       yyaxis(ax, 'left');
   else
       hADC = gobjects(0);
   end
   yyaxis(ax, 'left');
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'off');
    xlabel(ax, 'Time (ms)');
    ylabel(ax, 'EMG (\muV)');
    title(ax, { ...
        sprintf('%s | %s', fileName, EMG_field), ...
        'Mouse wheel/slider: scroll X. X +/-: horizontal zoom. Y +/-: vertical zoom.'}, ...
        'Interpreter','none');
    xlim(ax, [0 initialWindow_ms]);
    yMin = min(x_uV);
    yMax = max(x_uV);
    if yMin == yMax
        yMin = yMin - 1;
        yMax = yMax + 1;
    end
    pad = 0.1 * (yMax - yMin);
    ylim(ax, [yMin-pad, yMax+pad]);
    % Visible on demand; never modifies x_uV or any measurement.
   hADCToggle = uicontrol('Parent',fig, 'Style','checkbox', ...
       'String','Show PNS stimulus (ADC1_2)', ...
       'Units','normalized','Position',[0.63 0.095 0.30 0.055], ...
       'Value',0,'BackgroundColor','w', ...
       'Enable', 'on', ...
       'Callback',@toggleADC);
   if isempty(adcSignal), set(hADCToggle,'Enable','off'); end
   %% ===== Buttons =====
    uicontrol('Style','pushbutton', ...
        'String','Manual M-wave Selection', ...
        'Units','normalized', ...
        'Position',[0.07 0.17 0.18 0.06], ...
        'FontWeight','bold', ...
        'Callback',@selectMwave);
    uicontrol('Style','pushbutton', ...
        'String','Undo Last Edit', ...
        'Units','normalized', ...
        'Position',[0.27 0.17 0.15 0.06], ...
        'Callback',@undoLastEdit);
    uicontrol('Style','pushbutton', ...
        'String','Full view', ...
        'Units','normalized', ...
        'Position',[0.44 0.17 0.15 0.06], ...
        'Callback',@fullView);
    uicontrol('Style','pushbutton', ...
        'String','MS view', ...
        'Units','normalized', ...
        'Position',[0.61 0.17 0.12 0.06], ...
        'Callback',@msView);
    uicontrol('Style','pushbutton', ...
        'String','Done', ...
        'Units','normalized', ...
        'Position',[0.78 0.17 0.15 0.06], ...
        'Callback',@doneTool);
    % Automatic detection and log controls intentionally removed.
    %% Independent axis zoom buttons (keep the other axis fixed)
    btnLabels = {'X +','X -','Y +','Y -','Reset Y'};
    btnCalls = {@xZoomIn,@xZoomOut,@yZoomIn,@yZoomOut,@resetY};
    btnPos = [0.07 0.17 0.27 0.37 0.47];
    btnW = [0.085 0.085 0.085 0.085 0.13];
    for k = 1:numel(btnLabels)
        uicontrol('Parent',fig,'Style','pushbutton', ...
            'String',btnLabels{k},'Units','normalized', ...
            'Position',[btnPos(k) 0.095 btnW(k) 0.055], ...
            'Callback',btnCalls{k});
    end
    infoText = gobjects(0);
    %% ===== Horizontal scroll bar =====
    if totalDuration_ms > initialWindow_ms
        sliderMax = totalDuration_ms - initialWindow_ms;
    else
        sliderMax = 1;
    end
    hSlider = uicontrol('Style','slider', ...
        'Units','normalized', ...
        'Position',[0.07 0.035 0.86 0.035], ...
        'Min',0, ...
        'Max',sliderMax, ...
        'Value',0, ...
        'SliderStep',[min(scrollStep_ms/sliderMax,1), min(initialWindow_ms/sliderMax,1)], ...
        'Callback',@sliderMoved);
    %% ===== Store data =====
    S = struct();
    S.hADC = hADC;
   S.hADCToggle = hADCToggle;
   S.fig = fig;
    S.ax = ax;
    S.hEMG = hEMG;
    S.original_uV = x_uV;
    S.methodName = 'Original';
    S.infoText = infoText;
    S.hSlider = hSlider;
    S.t_ms = t_ms;
    S.x_uV = x_uV;
    S.Fs = Fs;
    S.fileName = fileName;
    S.EMG_field = EMG_field;
    S.initialWindow_ms = initialWindow_ms;
    S.scrollStep_ms = scrollStep_ms;
    S.totalDuration_ms = totalDuration_ms;
    S.originalYLim = ylim(ax);
    S.Mwave = [];
    S.tStart_ms = [];
    S.tEnd_ms = [];
    S.Duration_ms = [];
    S.P2P_uV = [];
    S.Max_uV = [];
    S.Min_uV = [];
    S.tMax_ms = [];
    S.tMin_ms = [];
    S.handles = {};
    S.showDuration = true;
    S.periodHandles = {};
    S.periodDrag = 0;
    S.periodOn = NaN; S.periodOff = NaN;
    S.periodRecords = struct('Pulse_ms',{},'Onset_ms',{},'Offset_ms',{},'Source',{},'Method',{},'MeasurementIndex',{});
    S.activePeriodRecord = 0;
    S.correctedWindows = zeros(0,2); % corrected time ranges [start_ms, end_ms]
    S.liuHistory = {}; % exact EMG samples before each applied correction
    S.pnsHandle = gobjects(0);
    S.pnsDrag = false;
    S.pnsOverrides = zeros(0,2); % [original ADC edge, adjusted time]
    S.undoStack = {}; % full chronological undo of every edit
    S.artifactEnds = zeros(0,2); % [PNS time, detected artifact endpoint]
    guidata(fig, S);
    refreshPNSMarker();
    zoom(fig, 'off');
    pan(fig, 'off');
    fprintf('\n======================================================\n');
    fprintf('M-WAVE TOOL\n');
    fprintf('======================================================\n');
    fprintf('Time axis is in ms.\n');
    fprintf('Use the horizontal scroll bar to move through the signal.\n');
    fprintf('You can also use the mouse wheel to move left/right.\n');
    fprintf('Press "Manual M-wave Selection", then click START and END.\n');
    fprintf('Use "Undo Last Edit" if needed.\n');
    fprintf('======================================================\n\n');
    didFinish = false;
    uiwait(fig);
    if ~isgraphics(fig)
        % X / window close: cancel without exporting or overwriting files.
        Results = table();
        PeriodResults = table();
        return;
    end
    S = guidata(fig);
    didFinish = isfield(S,'didFinish') && S.didFinish;
    delete(fig);
    if ~didFinish
        Results = table();
        PeriodResults = table();
        return;
    end
    Results = table( ...
        repmat(string(fileName), numel(S.Mwave), 1), ...
        repmat(string(EMG_field), numel(S.Mwave), 1), ...
        S.Mwave(:), ...
        S.tStart_ms(:), ...
        S.tEnd_ms(:), ...
        S.Duration_ms(:), ...
        S.P2P_uV(:), ...
        S.Max_uV(:), ...
        S.Min_uV(:), ...
        S.tMax_ms(:), ...
        S.tMin_ms(:), ...
        'VariableNames', { ...
        'File', ...
        'EMG_Channel', ...
        'Mwave', ...
        'tStart_ms', ...
        'tEnd_ms', ...
        'Duration_ms', ...
        'P2P_uV', ...
        'Max_uV', ...
        'Min_uV', ...
        'tMax_ms', ...
        'tMin_ms'});
    % The original manual M-wave selection workflow also exports absolute AUC.
    if ~isempty(Results)
        aucManual=nan(height(Results),1);
        for ii=1:height(Results)
            mask=S.t_ms>=Results.tStart_ms(ii)&S.t_ms<=Results.tEnd_ms(ii);
            preMask=S.t_ms>=Results.tStart_ms(ii)-50 & ...
                    S.t_ms<=Results.tStart_ms(ii)-5;
            baseline=0;
            if nnz(preMask)>=5,baseline=median(S.x_uV(preMask));end
            if nnz(mask)>=2
                aucManual(ii)=trapz(S.t_ms(mask),abs(S.x_uV(mask)-baseline));
            end
        end
        Results.AbsoluteAUC_uVms=aucManual;
    end
    % Automatic periods are not generated in this manual-only version.
    PeriodResults=table('Size',[0 8], ...
        'VariableTypes',{'double','double','double','double','double','double','string','string'}, ...
        'VariableNames',{'PNS_ms','Onset_ms','Offset_ms','Duration_ms', ...
        'AbsoluteAUC_uVms','P2P_uV','BoundarySource','AppliedArtifactMethod'});
    %% ===================================================================
    function toggleADC(src, ~)
       D = guidata(fig);
       if isempty(D.hADC) || ~isgraphics(D.hADC)
           return;
       end
       yyaxis(D.ax,'right');
       if logical(get(src,'Value'))
           set(D.hADC,'Visible','on');
           D.ax.YColor = [0.85 0.2 0.15];
           % Restore right-axis numeric tick labels from current ticks.
           D.ax.YTickMode = 'auto';
           D.ax.YTickLabelMode = 'auto';
           ylabel(D.ax,'PNS trigger (ADC1_2, raw units)');
       else
           set(D.hADC,'Visible','off');
           D.ax.YTickLabel = [];
           ylabel(D.ax,'');
       end
       yyaxis(D.ax,'left');
   end
   %% ========================= CALLBACKS ===============================
    %% ===================================================================
    function xZoomIn(~, ~)
        scaleX(0.5);
    end
    function xZoomOut(~, ~)
        scaleX(2.0);
    end
    function yZoomIn(~, ~)
        scaleY(0.5);
    end
    function yZoomOut(~, ~)
        scaleY(2.0);
    end
    function resetY(~, ~)
        D = guidata(fig);
        yyaxis(D.ax,'left');
       ylim(D.ax, D.originalYLim);
    end
    function scaleX(factor)
        D = guidata(fig);
        xl = xlim(D.ax);
        width = min(D.totalDuration_ms, max(2 * 1000 / D.Fs, diff(xl) * factor));
        middle = mean(xl);
        left = min(max(middle - width / 2, D.t_ms(1)), D.totalDuration_ms - width);
        xlim(D.ax, [left left + width]);
        syncSlider(D);
        refreshPNSMarker();
    end
    function scaleY(factor)
        D = guidata(fig);
        yyaxis(D.ax,'left');
       yl = ylim(D.ax);
        midpoint = mean(yl);
        halfRange = diff(yl) * factor / 2;
        ylim(D.ax, [midpoint - halfRange midpoint + halfRange]);
    end
    function syncSlider(D)
        %% Slider range must follow the current zoomed time window.
        xl = xlim(D.ax);
        maxStart = max(0, D.totalDuration_ms - diff(xl));
        if maxStart <= 0
            set(D.hSlider,'Min',0,'Max',1,'Value',0,'Enable','off');
        else
            smallStep = min(max(D.scrollStep_ms / maxStart,0.001),1);
            pageStep = min(max(diff(xl) / maxStart,smallStep),1);
            set(D.hSlider,'Enable','on','Min',0,'Max',maxStart, ...
                'Value',min(max(xl(1),0),maxStart), ...
                'SliderStep',[smallStep pageStep]);
        end
    end
    function selectMwave(src, ~)
        fig = ancestor(src, 'figure');
        S = guidata(fig);
        zoom(fig, 'off');
        pan(fig, 'off');
        appendLog('Click START then END around the M-wave.');
        yyaxis(S.ax, 'left');
        axes(S.ax);
        fprintf('\nClick START of M-wave...\n');
        [xs, ~] = ginput(1);
        fprintf('Click END of M-wave...\n');
        [xe, ~] = ginput(1);
        a = min(xs, xe);
        b = max(xs, xe);
        a = max(a, S.t_ms(1));
        b = min(b, S.t_ms(end));
        idx = S.t_ms >= a & S.t_ms <= b;
        if sum(idx) < 2
            warning('Selected window too short. Try again.');
            appendLog('Selection too short. Try again.');
            zoom(fig, 'off');
            pan(fig, 'off');
            return;
        end
        rememberFrame(S);
        S=guidata(fig);
        seg = S.x_uV(idx);
        tseg = S.t_ms(idx);
        [maxVal, idxMax] = max(seg);
        [minVal, idxMin] = min(seg);
        p2p = maxVal - minVal;
        if isempty(S.Mwave)
            newNumber = 1;
        else
            newNumber = max(S.Mwave) + 1;
        end
        S.Mwave(end+1,1)       = newNumber;
        S.tStart_ms(end+1,1)   = a;
        S.tEnd_ms(end+1,1)     = b;
        S.Duration_ms(end+1,1) = b - a;
        S.P2P_uV(end+1,1)      = p2p;
        S.Max_uV(end+1,1)      = maxVal;
        S.Min_uV(end+1,1)      = minVal;
        S.tMax_ms(end+1,1)     = tseg(idxMax);
        S.tMin_ms(end+1,1)     = tseg(idxMin);
        yyaxis(S.ax,'left');
       yl = ylim(S.ax);
        hPatch = patch(S.ax, ...
            [a b b a], ...
            [yl(1) yl(1) yl(2) yl(2)], ...
            [0.85 0.85 0.85], ...
            'EdgeColor','none', ...
            'FaceAlpha',0.25);
        % Avoid uistack with yyaxis: it may fail with mixed chart objects.
        % The translucent patch is noninteractive and does not affect P2P.
        set(hPatch, 'HitTest', 'off');
        hStart = xline(S.ax, a, '--b', sprintf('M%d start', newNumber), ...
            'LabelVerticalAlignment','bottom');
        hEnd = xline(S.ax, b, '--r', sprintf('M%d end', newNumber), ...
            'LabelVerticalAlignment','bottom');
        hMax = plot(S.ax, tseg(idxMax), maxVal, 'ro', ...
            'MarkerFaceColor','r');
        hMin = plot(S.ax, tseg(idxMin), minVal, 'bo', ...
            'MarkerFaceColor','b');
        hText = text(S.ax, a, yl(2), ...
            sprintf(' M%d: P2P %.1f \muV | duration %.2f ms', newNumber, p2p, b-a), ...
            'VerticalAlignment','top', ...
            'FontWeight','bold', ...
            'Color','k');
        % Keep handles in a cell array: Patch, ConstantLine, Line, and Text
        % are different MATLAB graphics classes and cannot safely be concatenated.
        S.handles{end+1} = {hPatch, hStart, hEnd, hMax, hMin, hText};
        set(hPatch,'Visible',mwOnOff(S.showDuration));
        set(hStart,'Visible',mwOnOff(S.showDuration));
        set(hEnd,'Visible',mwOnOff(S.showDuration));
        fprintf('\nM-wave %d selected:\n', newNumber);
        fprintf('  Start    = %.2f ms\n', a);
        fprintf('  End      = %.2f ms\n', b);
        fprintf('  Duration = %.2f ms\n', b-a);
        fprintf('  Max      = %.1f µV at %.2f ms\n', maxVal, tseg(idxMax));
        fprintf('  Min      = %.1f µV at %.2f ms\n', minVal, tseg(idxMin));
        fprintf('  P2P      = %.1f µV\n\n', p2p);
        appendLog(sprintf( ...
            'Selected M-wave %d | P2P = %.1f µV | Continue scrolling/zooming or select next M-wave.', ...
            newNumber, p2p));
        guidata(fig, S);
        zoom(fig, 'off');
        pan(fig, 'off');
    end
    function appendLog(~)
        % No analysis log: retained only for existing non-UI callbacks.
    end
    function [pulse,edge,report]=currentPNSTime(D)
        pulse=NaN;edge=NaN;report='';
        if isempty(tADC) || isempty(adcSignal),report='No ADC1_2 channel.';return;end
        xl=xlim(D.ax);
        [pulses,report]=mwDetectVisibleADCPulses(tADC(:)*1000,adcSignal(:),xl);
        if isempty(pulses),return;end
        [~,k]=min(abs(pulses-mean(xl)));edge=pulses(k);pulse=edge+0;
        if ~isempty(D.pnsOverrides)
            [delta,j]=min(abs(D.pnsOverrides(:,1)-edge));
            if delta<=5,pulse=D.pnsOverrides(j,2);end % Adjustments store manually positioned PNS time
        end
    end
    function refreshPNSMarker()
        if ~isgraphics(fig),return;end
        D=guidata(fig);
        if isfield(D,'pnsDrag') && D.pnsDrag,return;end
        if isfield(D,'pnsHandle')
            validMarkers=D.pnsHandle(isgraphics(D.pnsHandle));
            if ~isempty(validMarkers),delete(validMarkers);end
        end
        D.pnsHandle=gobjects(0);
        [pulse,~,~]=currentPNSTime(D);
        if isfinite(pulse)
            yyaxis(D.ax,'left');
            D.pnsHandle=xline(D.ax,pulse,':','PNS','Color',[.9 .08 .08], ...
                'LineWidth',1.8,'HitTest','on','ButtonDownFcn',@beginPNSDrag);
        end
        guidata(fig,D);
    end
    function beginPNSDrag(~,~)
        D=guidata(fig);
        [~,edge,~]=currentPNSTime(D);
        if ~isfinite(edge),return;end
        rememberFrame(D);D=guidata(fig);
        D.pnsDrag=true;D.pnsDragEdge=edge;
        guidata(fig,D);
        set(fig,'WindowButtonMotionFcn',@dragPNS,'WindowButtonUpFcn',@finishPNSDrag);
    end
    function dragPNS(~,~)
        D=guidata(fig);
        if ~D.pnsDrag || ~isscalar(D.pnsHandle) || ~all(isgraphics(D.pnsHandle)),return;end
        cp=get(D.ax,'CurrentPoint');xl=xlim(D.ax);
        value=min(max(cp(1,1),xl(1)),xl(2));
        set(D.pnsHandle,'Value',value);
    end
    function finishPNSDrag(~,~)
        D=guidata(fig);
        if ~D.pnsDrag,return;end
        D.pnsDrag=false;
        if isscalar(D.pnsHandle) && all(isgraphics(D.pnsHandle))
            newTime=get(D.pnsHandle,'Value');
            j=find(abs(D.pnsOverrides(:,1)-D.pnsDragEdge)<0.001,1);
            if isempty(j),D.pnsOverrides(end+1,:)=[D.pnsDragEdge newTime];
            else,D.pnsOverrides(j,2)=newTime;end
            appendLog(sprintf('PNS marker adjusted: ADC edge %.3f ms -> %.3f ms. Subsequent detection uses adjusted time.',D.pnsDragEdge,newTime));
        end
        guidata(fig,D);
        set(fig,'WindowButtonMotionFcn',[],'WindowButtonUpFcn',[]);
    end
    function rememberFrame(D)
        % Every individual edit gets its own independent complete pre-edit state.
        snap=struct();
        names={'x_uV','Mwave','tStart_ms','tEnd_ms','Duration_ms', ...
            'P2P_uV','Max_uV','Min_uV','tMax_ms','tMin_ms', ...
            'periodRecords','correctedWindows','liuHistory','methodName', ...
            'artifactEnds','pnsOverrides'};
        for q=1:numel(names),snap.(names{q})=D.(names{q});end
        D.undoStack{end+1}=snap;
        guidata(fig,D);
    end

    function undoLastEdit(~,~)
        D=guidata(fig);
        if isempty(D.undoStack)
            appendLog('Undo: no remaining edits.');return;
        end
        snap=D.undoStack{end};D.undoStack(end)=[];
        names=fieldnames(snap);
        for q=1:numel(names),D.(names{q})=snap.(names{q});end
        D.pnsDrag=false;D.periodDrag=0;D.activePeriodRecord=0;
        D.periodOn=NaN;D.periodOff=NaN;
        set(fig,'WindowButtonMotionFcn',[],'WindowButtonUpFcn',[]);
        set(D.hEMG,'YData',D.x_uV);
        D=redrawMeasurements(D);
        guidata(fig,D);
        refreshPNSMarker();
        appendLog(sprintf('Undo: restored previous signal and selections (%d edits remain).',numel(D.undoStack)));
    end

    function D=redrawMeasurements(D)
        % Delete all selection annotations and reconstruct from restored arrays.
        for q=1:numel(D.handles)
            h=D.handles{q};
            for r=1:numel(h),if isgraphics(h{r}),delete(h{r});end,end
        end
        for q=1:numel(D.periodHandles)
            if isgraphics(D.periodHandles{q}),delete(D.periodHandles{q});end
        end
        D.handles={};D.periodHandles={};
        for q=1:numel(D.Mwave)
            a=D.tStart_ms(q);b=D.tEnd_ms(q);
            idx=find(D.t_ms>=a&D.t_ms<=b);
            if numel(idx)<2,D.handles{end+1}={};continue;end
            yyaxis(D.ax,'left');yl=ylim(D.ax);
            [mx,imx]=max(D.x_uV(idx));[mn,imn]=min(D.x_uV(idx));
            tt=D.t_ms(idx);
            hp=patch(D.ax,[a b b a],[yl(1) yl(1) yl(2) yl(2)], ...
                [.85 .85 .85],'EdgeColor','none','FaceAlpha',.20,'HitTest','off');
            hs=xline(D.ax,a,'--b',sprintf('M%d start',D.Mwave(q)));
            he=xline(D.ax,b,'--r',sprintf('M%d end',D.Mwave(q)));
            hmx=plot(D.ax,tt(imx),mx,'ro','MarkerFaceColor','r');
            hmn=plot(D.ax,tt(imn),mn,'bo','MarkerFaceColor','b');
            ht=text(D.ax,a,yl(2),sprintf(' M%d: P2P %.1f \muV | duration %.2f ms', ...
                D.Mwave(q),D.P2P_uV(q),b-a),'VerticalAlignment','top','FontWeight','bold');
            D.handles{end+1}={hp,hs,he,hmx,hmn,ht};
        end
    end

    function sliderMoved(src, ~)
        fig = ancestor(src, 'figure');
        S = guidata(fig);
        currentWindow = diff(xlim(S.ax));
        newStart = get(S.hSlider, 'Value');
        newEnd = newStart + currentWindow;
        if newEnd > S.totalDuration_ms
            newEnd = S.totalDuration_ms;
            newStart = max(0, newEnd - currentWindow);
        end
        xlim(S.ax, [newStart newEnd]);
            syncSlider(S);
        guidata(fig, S);
        refreshPNSMarker();
        appendLog(sprintf('Current view: %.1f to %.1f ms', newStart, newEnd));
        guidata(fig, S);
    end
    function mouseScroll(fig, event)
        S = guidata(fig);
        xl = xlim(S.ax);
        currentWindow = diff(xl);
        shift = event.VerticalScrollCount * S.scrollStep_ms;
        newStart = xl(1) + shift;
        newEnd = xl(2) + shift;
        if newStart < 0
            newStart = 0;
            newEnd = currentWindow;
        end
        if newEnd > S.totalDuration_ms
            newEnd = S.totalDuration_ms;
            newStart = max(0, newEnd - currentWindow);
        end
        xlim(S.ax, [newStart newEnd]);
            syncSlider(S);
        guidata(fig, S);
        refreshPNSMarker();
        appendLog(sprintf('Current view: %.1f to %.1f ms', newStart, newEnd));
        guidata(fig, S);
    end
    function fullView(src, ~)
        fig = ancestor(src, 'figure');
        S = guidata(fig);
        xlim(S.ax, [S.t_ms(1), S.t_ms(end)]);
        syncSlider(S);
        yMin = min(S.x_uV);
        yMax = max(S.x_uV);
        if yMin == yMax
            yMin = yMin - 1;
            yMax = yMax + 1;
        end
        pad = 0.1 * (yMax - yMin);
        yyaxis(S.ax,'left');
       ylim(S.ax, [yMin-pad, yMax+pad]);
        guidata(fig, S);
        refreshPNSMarker();
        appendLog('Full signal view.');
        guidata(fig, S);
    end
    function msView(src, ~)
        fig = ancestor(src, 'figure');
        S = guidata(fig);
        xl = xlim(S.ax);
        currentStart = xl(1);
        newEnd = currentStart + S.initialWindow_ms;
        if newEnd > S.totalDuration_ms
            newEnd = S.totalDuration_ms;
            currentStart = max(0, newEnd - S.initialWindow_ms);
        end
        xlim(S.ax, [currentStart newEnd]);
            syncSlider(S);
        guidata(fig, S);
        refreshPNSMarker();
        appendLog(sprintf('MS view: %.1f to %.1f ms', currentStart, newEnd));
        guidata(fig, S);
    end
    function cancelTool(src, ~)
        % Closing the window is cancellation, never equivalent to Done.
        if isgraphics(src)
            delete(src); % Releases uiwait and aborts the export in the caller.
        end
    end
    function doneTool(src, ~)
        fig = ancestor(src, 'figure');
        S = guidata(fig);
        fprintf('\nSelection finished. Total selected M-waves: %d\n', numel(S.Mwave));
        S.didFinish = true;
        guidata(fig, S);
        uiresume(fig);
    end
end

%% ===== Data-driven early artifact candidate detection =====
function [times,diagnostic] = mwDetectADCPulses(t_ms,adc)
    % Identify brief, high-slope ADC1_2 transitions. Cannot infer true
    % stimulus artifact duration or physiological response onset from ADC.
    times=[]; diagnostic='';
    t_ms=t_ms(:); adc=adc(:);
    good=isfinite(t_ms)&isfinite(adc);
    t_ms=t_ms(good); adc=adc(good);
    if numel(adc)<10 || any(diff(t_ms)<=0)
        diagnostic='ADC samples are missing or time is not monotonic.';return;
    end
    da=diff(adc);
    absDa=abs(da);
    med=median(absDa);
    robustSigma=1.4826*median(abs(absDa-med));
    pk=max(absDa);
    if pk<=max(eps(max(abs(adc))),1e-10)
        diagnostic='ADC trace is nearly constant.';return;
    end
    % Both relative and noise-dependent limits guard against noise edges.
    cutoff=max(0.28*pk,med+8*robustSigma);
    idx=find(absDa>=cutoff);
    if isempty(idx)
        diagnostic='No sharp ADC transition exceeds threshold.';return;
    end
    % Edges from one electrical pulse may occur within milliseconds.
    % Prefer first transition of each cluster as candidate pulse onset.
    keep=[true;diff(t_ms(idx))>15];
    idx=idx(keep);
    times=t_ms(idx+1);
    diagnostic=sprintf('%d candidate stimulus edges.',numel(times));
end

function [times,diagnostic] = mwDetectVisibleADCPulses(t_ms,adc,view)
    % Event detection is performed locally, NOT relative to the largest
    % ADC transition elsewhere in the recording. Trigger only; not EMG.
    times=[];
    t_ms=t_ms(:); adc=double(adc(:));
    valid=isfinite(t_ms)&isfinite(adc);
    t_ms=t_ms(valid); adc=adc(valid);
    if numel(t_ms)<10
        diagnostic='PNS: ADC1_2 has insufficient samples.';return;
    end
    % Include nearby samples for event edge detection when view cuts it.
    within=t_ms>=view(1)-5 & t_ms<=view(2)+5;
    tt=t_ms(within); xx=adc(within);
    if numel(tt)<10
        diagnostic='PNS: no ADC samples in the visible time range.';return;
    end
    d=abs(diff(xx));
    positive=d(isfinite(d));
    if isempty(positive)||max(positive)<1e-10
        diagnostic='PNS: local ADC1_2 is flat.';return;
    end
    noiseMedian=median(positive);
    noiseSD=1.4826*median(abs(positive-noiseMedian));
    limit=max(0.25*max(positive), noiseMedian+8*noiseSD);
    candidates=find(d>=limit);
    if isempty(candidates)
        diagnostic=sprintf('PNS: no ADC edge passed local threshold %.4g raw units/sample.',limit);return;
    end
    candidateTimes=tt(candidates+1);
    % One event may generate a falling edge and a rebound. Group within 15 ms.
    keep=[true; diff(candidateTimes)>15];
    candidateTimes=candidateTimes(keep);
    times=candidateTimes(candidateTimes>=view(1) & candidateTimes<=view(2));
    if isempty(times)
        diagnostic=sprintf('PNS: %d local edge cluster(s), none in current %.1f–%.1f ms view (threshold %.4g).', ...
            numel(candidateTimes),view(1),view(2),limit);
    else
        diagnostic=sprintf('PNS: %d visible pulse(s), at %s ms (local ADC slope threshold %.4g).', ...
            numel(times),strjoin(compose('%.2f',times),', '),limit);
    end
end

function visible = mwOnOff(state)
    if state, visible='on'; else, visible='off'; end
end
