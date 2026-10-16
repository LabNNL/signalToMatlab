function [MEPnew] = OnsetOffset_Inspection(MEP)

mepFields = fieldnames(MEP);
mepNames = mepFields(startsWith(mepFields,'MEP'));
nb_mep = numel(mepNames);

fig = uifigure('Name','Onset-Offset Inspection','Position',[100 100 900 600]);
currentIndex = 1;

decal = 40;                                                                 % Used to align sliders with the axes
modifmax=50;                                                                % Maximum manual offset
fitYWindow = false;                                                         % Whether the Y axis is fit to the -100/100ms window only
yFitRange = [-100 100];
envelopeVisible = true;                                                     % Whether the RMS envelope curve is shown
showRawSignal = false;                                                      % Whether the raw MEP signal is overlaid on the envelope

axesPos = [20 160 680 400];
btnX = axesPos(1)+axesPos(3)+30;

% MEP data stored in arrays
x = MEP.Meta.Time_ms;
onsetArr = zeros(nb_mep,1);
offsetArr = zeros(nb_mep,1);
spArr = zeros(nb_mep,1);
limOnArr = zeros(nb_mep,1);
limOffArr = zeros(nb_mep,1);
hasRaw = false(nb_mep,1);
includedArr = true(nb_mep,1);

    for k = 1:nb_mep
        hasRaw(k) = isfield(MEP.(mepNames{k}),'EMG');
        modifmax_on = modifmax;
        modifmax_off = modifmax;

        onset=MEP.Meta.OnOff_ms(k,1);
        offset=MEP.Meta.OnOff_ms(k,2);
        sp=offset+MEP.(mepNames{k}).Silentperiod;
        onset(isnan(onset)) = 0;
        offset(isnan(offset)) = 0;
        sp(isnan(sp)) = 0;
        sp(sp>400) = 400;

        if onset == 0
            modifmax_on = 100;
        end
        if offset == 0
            modifmax_off = 150;
        end

        onsetArr(k) = onset;
        offsetArr(k) = offset;
        spArr(k) = sp;
        limOnArr(k) = onset+modifmax_on;
        limOffArr(k) = offset+modifmax_off;
    end

    ax = uiaxes(fig,'Position',axesPos);
    ax.Box = 'on';

    yyaxis(ax,'left');
    plotEnv = plot(ax,x,MEP.(mepNames{1}).Enveloppe,'Color','k','LineWidth',1);
    ax.XLim = [min(x) 400];
    ax.YAxis(1).Color = 'k';
    ylabel(ax,'RMS Envelope (V)');

    plotRaw = gobjects(1);
    ylineRawZero = gobjects(1);
    if any(hasRaw)
        yyaxis(ax,'right');
        plotRaw = plot(ax,x,NaN(size(x)),'Visible','off','Color',[0.6 0.6 0.6],'LineWidth',0.8);
        ylabel(ax,'Raw MEP signal (V)');
        ax.YAxis(2).Color = [0.6 0.6 0.6];
        ylineRawZero = yline(ax, 0, 'Color', [0.85 0.85 0.85], 'LineWidth', 1, 'LineStyle', '--', 'Visible', 'off');
        yyaxis(ax,'left');
    end

    xlabel(ax,'Time (ms)');
    title(ax,mepNames{1},'Interpreter','none');

    % xlines - stim, onset, offset & silent period
    xline(ax, 0, 'Color', 'r', 'LineWidth', 1.5, 'LineStyle', ':', 'Label', 'Stimulation', 'LabelHorizontalAlignment', 'left');
    xline1 = xline(ax, onsetArr(1), 'Color', '#6B43E5', 'LineWidth', 1.7, 'LineStyle', '--');
    xline2 = xline(ax, offsetArr(1), 'Color', '#E54379', 'LineWidth', 1.7, 'LineStyle', '--');
    xline3 = xline(ax, spArr(1), 'Color', '#AA50DE', 'LineWidth', 1.7, 'LineStyle', '--');

    % ylines - onset & offset thresholds
    yline1 = yline(ax, MEP.(mepNames{1}).Thresholds.on, 'Color', '#6B43E5', 'LineWidth', .7, 'LineStyle', '--');
    yline2 = yline(ax, MEP.(mepNames{1}).Thresholds.off, 'Color', '#E54379', 'LineWidth', .7, 'LineStyle', '--');

    % Slider 1 (Onset)
    slider1 = uislider(fig, ...
        'Position',[axesPos(1)+decal 120 axesPos(3)-decal 3], ...
        'Limits',[0 limOnArr(1)], ...
        'Value',onsetArr(1));
    slider1.ValueChangingFcn = @(s,e) updateXline1(e.Value);
    slider1.ValueChangedFcn = @(s,e) updateXline1(e.Value);

    % Slider 2 (Offset)
    slider2 = uislider(fig, ...
        'Position',[axesPos(1)+decal 80 axesPos(3)-decal 3], ...
        'Limits',[0 limOffArr(1)], ...
        'Value',offsetArr(1));
    slider2.ValueChangingFcn = @(s,e) updateXline2(e.Value);
    slider2.ValueChangedFcn = @(s,e) updateXline2(e.Value);

    % Slider 3 (Silent Period)
    slider3 = uislider(fig, ...
        'Position',[axesPos(1)+decal 40 axesPos(3)-decal 3], ...
        'Limits',[0 400], ...
        'Value',spArr(1));
    slider3.ValueChangingFcn = @(s,e) updateXline3(e.Value);
    slider3.ValueChangedFcn = @(s,e) updateXline3(e.Value);

    uilabel(fig,'Text','Onset','Position',[axesPos(1)-10 112 50 20],'FontWeight','bold','FontColor','#6B43E5');
    uilabel(fig,'Text','Offset','Position',[axesPos(1)-10 72 50 20],'FontWeight','bold','FontColor','#E54379');
    uilabel(fig,'Text','Sil.Per.','Position',[axesPos(1)-10 32 50 20],'FontWeight','bold','FontColor','#AA50DE');

    % Checkbox - include/exclude this MEP from the exported structure
    checkboxIncluded = uicheckbox(fig, ...
        'Text','Included if checked', ...
        'Position',[btnX 405 150 22], ...
        'Value',true, ...
        'FontWeight','bold', ...
        'ValueChangedFcn',@(src,evt) updateIncluded(src.Value));

    checkboxRowY = axesPos(2)+axesPos(4)+5;

    % Checkbox - show/hide the RMS envelope
    uicheckbox(fig, ...
        'Text','Show MEP envelope', ...
        'Position',[axesPos(1) checkboxRowY 190 22], ...
        'Value',true, ...
        'ValueChangedFcn',@(src,evt) toggleEnvelope(src.Value));

    % Checkbox - overlay the raw MEP signal
    uicheckbox(fig, ...
        'Text','Show raw MEP signal', ...
        'Position',[axesPos(1)+200 checkboxRowY 200 22], ...
        'Value',false, ...
        'ValueChangedFcn',@(src,evt) toggleShowRaw(src.Value));

    % Checkbox - fit the Y axis to the MEP window only
    uicheckbox(fig, ...
        'Text','Fit Y axis on MEP window (-100 to 100 ms)', ...
        'Position',[axesPos(1)+410 checkboxRowY 400 22], ...
        'Value',false, ...
        'ValueChangedFcn',@(src,evt) toggleFitY(src.Value));

    % Initial display
    showCurve(1);

    % Buttons to switch curves
    uibutton(fig,'Text','⏮️ First','Position',[btnX 350 120 40],'ButtonPushedFcn',@(btn,event) jumpTo(1));
    uibutton(fig,'Text','⬅️ Previous','Position',[btnX 280 120 40],'ButtonPushedFcn',@(btn,event) switchCurve(-1));
    uibutton(fig,'Text','Next ➡️','Position',[btnX 210 120 40],'ButtonPushedFcn',@(btn,event) switchCurve(1));
    uibutton(fig,'Text','Last ⏭️','Position',[btnX 140 120 40],'ButtonPushedFcn',@(btn,event) jumpTo(nb_mep));
    uibutton(fig,'Text','Finish','Position',[btnX 70 120 40], ...
    'FontColor','r','BackgroundColor',fig.Color,'ButtonPushedFcn',@(btn,event) finishCallback());

    lblInfo = uilabel(fig,'Text',sprintf('MEP %d / %d',currentIndex,nb_mep), ...
    'Position',[btnX 440 130 30],'FontSize',14,'FontWeight','bold');

%%  Functions
%%
    function updateXline1(val)
        onsetArr(currentIndex) = val;
        xline1.Value = val;
    end

    function updateXline2(val)
        offsetArr(currentIndex) = val;
        xline2.Value = val;
    end

    function updateXline3(val)
        spArr(currentIndex) = val;
        xline3.Value = val;
    end

    function updateIncluded(val)
        includedArr(currentIndex) = logical(val);
    end

    function saveCurrent()
        onsetArr(currentIndex) = xline1.Value;
        offsetArr(currentIndex) = xline2.Value;
        spArr(currentIndex) = xline3.Value;
        includedArr(currentIndex) = logical(checkboxIncluded.Value);
    end

    function setSlider(s, lim, val)
        s.Value = 0;
        s.Limits = [0 lim];
        s.Value = val;
    end

    function showCurve(idx)
        % Load the data of MEP idx into the shared axes, lines and sliders
        mepData = MEP.(mepNames{idx});

        title(ax,mepNames{idx},'Interpreter','none');
        plotEnv.YData = mepData.Enveloppe;
        if isgraphics(plotRaw)
            if hasRaw(idx)
                plotRaw.YData = mepData.EMG;
            else
                plotRaw.YData = NaN(size(x));
            end
        end
        ax.XLim = [min(x) 400];

        yline1.Value = mepData.Thresholds.on;
        yline2.Value = mepData.Thresholds.off;
        xline1.Value = onsetArr(idx);
        xline2.Value = offsetArr(idx);
        xline3.Value = spArr(idx);

        setSlider(slider1, limOnArr(idx), onsetArr(idx));
        setSlider(slider2, limOffArr(idx), offsetArr(idx));
        setSlider(slider3, 400, spArr(idx));
        checkboxIncluded.Value = includedArr(idx);

        updateEnvelopeVisibility();
        updateRawVisibility(idx);
        applyYFit(idx);
    end

    function toggleFitY(value)
        fitYWindow = logical(value);
        applyYFit(currentIndex);
    end

    function toggleShowRaw(value)
        showRawSignal = logical(value);
        updateRawVisibility(currentIndex);
    end

    function toggleEnvelope(value)
        envelopeVisible = logical(value);
        updateEnvelopeVisibility();
    end

    function updateEnvelopeVisibility()
        onOff = onOffText(envelopeVisible);
        plotEnv.Visible = onOff;
        yline1.Visible = onOff;
        yline2.Visible = onOff;
        ax.YAxis(1).Visible = 'on';
        ax.YAxis(1).Label.Visible = onOff;
        if envelopeVisible
            ax.YAxis(1).TickValuesMode = 'auto';
            ax.YAxis(1).TickLabelsMode = 'auto';
        else
            ax.YAxis(1).TickValues = [];
        end
    end

    function updateRawVisibility(idx)
        if ~isgraphics(plotRaw)
            return
        end
        showRaw = showRawSignal && hasRaw(idx);
        plotRaw.Visible = onOffText(showRaw);
        ylineRawZero.Visible = onOffText(showRaw);
        ax.YAxis(2).Visible = 'on';
        ax.YAxis(2).Label.Visible = onOffText(showRaw);
        if showRaw
            ax.YAxis(2).Color = [0.6 0.6 0.6];
            ax.YAxis(2).TickValuesMode = 'auto';
            ax.YAxis(2).TickLabelsMode = 'auto';
        else
            ax.YAxis(2).Color = 'k';
            ax.YAxis(2).TickValues = [];
        end
    end

    function s = onOffText(tf)
        if tf
            s = 'on';
        else
            s = 'off';
        end
    end

    function applyYFit(idx)
        yyaxis(ax,'left');
        if ~fitYWindow
            ax.YLimMode = 'auto';
        else
            fitAxisToWindow(plotEnv);
        end

        if isgraphics(plotRaw)
            yyaxis(ax,'right');
            if ~fitYWindow || ~hasRaw(idx)
                ax.YLimMode = 'auto';
            else
                fitAxisToWindow(plotRaw);
            end
            yyaxis(ax,'left');
        end

        function fitAxisToWindow(plotObj)
            xdata = plotObj.XData;
            ydata = plotObj.YData;
            winMask = xdata >= yFitRange(1) & xdata <= yFitRange(2);
            windowData = ydata(winMask);

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
    end

    function switchCurve(dir)
        saveCurrent();
        currentIndex = currentIndex + dir;
        if currentIndex > nb_mep, currentIndex = 1; end
        if currentIndex < 1, currentIndex = nb_mep; end
        showCurve(currentIndex);
        lblInfo.Text = sprintf('MEP %d / %d',currentIndex,nb_mep);
    end

    function jumpTo(idx)
        saveCurrent();
        currentIndex = idx;
        showCurve(currentIndex);
        lblInfo.Text = sprintf('MEP %d / %d',currentIndex,nb_mep);
    end

    function finishCallback()
        saveCurrent();
        positions = [onsetArr offsetArr];
        positions_idx = zeros(nb_mep,2);
        MEPnew=MEP;
        for i = 1:nb_mep
            positions_idx(i,1) = find(abs(x - onsetArr(i)) == min(abs(x - onsetArr(i))), 1);
            positions_idx(i,2) = find(abs(x - offsetArr(i)) == min(abs(x - offsetArr(i))), 1);
            MEPnew.(mepNames{i}).OnOff_ms=positions(i,:);
            MEPnew.(mepNames{i}).OnOff_idx=positions_idx(i,:);
            MEPnew.(mepNames{i}).Silentperiod=spArr(i)-offsetArr(i);
            MEPnew.(mepNames{i}).Included=includedArr(i);
        end
        MEPnew.Meta.OnOff_ms=positions;
        uiresume(fig);
        delete(fig);
    end

    uiwait(fig);

end
