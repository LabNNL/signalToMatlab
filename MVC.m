%% MVC PRE/POST (entire selected window) 
% What it does:
% 1) Pick PRE file + POST file (CFS or MAT)
% 2) For each file: click START then END for each MVC (default 3 MVCs)
% 3) For each selected MVC window, computes on the ENTIRE window:
%    - P2P (µV)
%    - RMS (µV)  [recommended primary for MVC EMG amplitude]
%    - MAV (mean absolute value, µV)
%    - iEMG (integral of rectified EMG, µV·s)
%    - Duration (s)
% 4) Summarizes across 3 MVCs (mean±SD)
% 5) Computes POST-PRE change (absolute and %), using Mean RMS as primary
% 6) Saves CSV + MAT (in PRE folder)

clc; clear; close all
addpath("functions");
addpath("matcfs64c");

%% ===== User settings =====
defaultDir   = "D:\POST-DOC\PHASE2\CP participants\PAS1_EXO_01(Idris)\post-pre\right leg\TA\MVC";  % <<< change
EMG_field    = "EMG_6";                     % <<< change
doFilter     = true;                        % uses filtrage(20–450 Hz) if available
nTrials      = 3;                           % MVC repetitions

%% ===== Pick PRE and POST files =====
[preFile, preDir] = uigetfile( ...
    {'*.cfs;*.mat','CFS or MAT (*.cfs, *.mat)'; '*.cfs','CFS (*.cfs)'; '*.mat','MAT (*.mat)'}, ...
    'Select PRE file', defaultDir);
if isequal(preFile,0), error('No PRE file selected.'); end
prePath = fullfile(preDir, preFile);

[postFile, postDir] = uigetfile( ...
    {'*.cfs;*.mat','CFS or MAT (*.cfs, *.mat)'; '*.cfs','CFS (*.cfs)'; '*.mat','MAT (*.mat)'}, ...
    'Select POST file', preDir);
if isequal(postFile,0), error('No POST file selected.'); end
postPath = fullfile(postDir, postFile);

%% ===== Process PRE =====
[PreTrials, PreSummary, preMeta] = processOneMVCfile_entire(prePath, EMG_field, doFilter, nTrials, "PRE");

%% Keep PRE and POST preprocessing identical if the user chose raw EMG after a filter failure.
doFilter = preMeta.doFilter;

%% ===== Process POST =====
[PostTrials, PostSummary, postMeta] = processOneMVCfile_entire(postPath, EMG_field, doFilter, nTrials, "POST");

%% ===== Combine + compute PRE/POST change (primary = Mean RMS across trials) =====
preVal  = PreSummary.Mean_RMS_uV;
postVal = PostSummary.Mean_RMS_uV;

deltaAbs = postVal - preVal;
deltaPct = 100 * (postVal - preVal) / preVal;

Change = table(preVal, postVal, deltaAbs, deltaPct, ...
    'VariableNames', {'PRE_MeanRMS_uV','POST_MeanRMS_uV','Delta_uV','Delta_percent'});

disp("=== PRE summary ===");  disp(PreSummary);
disp("=== POST summary ==="); disp(PostSummary);
disp("=== Change (POST - PRE) ==="); disp(Change);

% ===== Save outputs (in PRE folder) =====
saveDir = preDir;
[~, preBase]  = fileparts(preFile);
[~, postBase] = fileparts(postFile);

% outCSV_trials = fullfile(saveDir, sprintf('%s_to_%s_%s_MVC_trials.csv', preBase, postBase, EMG_field));
% writetable([PreTrials; PostTrials], outCSV_trials);
% 
% outCSV_sum = fullfile(saveDir, sprintf('%s_to_%s_%s_MVC_summary.csv', preBase, postBase, EMG_field));
% writetable([PreSummary; PostSummary], outCSV_sum);

outCSV_change = fullfile(saveDir, sprintf('%s_to_%s_%s_MVC_change.csv', preBase, postBase, EMG_field));
writetable(Change, outCSV_change);

% outMAT = fullfile(saveDir, sprintf('%s_to_%s_%s_MVC.mat', preBase, postBase, EMG_field));
% save(outMAT, ...
%     'PreTrials','PostTrials','PreSummary','PostSummary','Change', ...
%     'preMeta','postMeta','EMG_field','doFilter','nTrials', ...
%     '-v7.3');

%fprintf('\nSaved:\n  %s\n  %s\n  %s\n', outCSV_trials, outCSV_sum, outMAT);

%% =======================================================================
%% ========================= LOCAL FUNCTIONS =============================
%% =======================================================================

function [TrialsTable, SummaryTable, meta] = processOneMVCfile_entire(fullPath, EMG_field, doFilter, nTrials, label)
    % Load EMG
    [t, x, Fs, usedFilter] = loadOneEMG(fullPath, EMG_field, doFilter);

    % Select segments + compute metrics (ENTIRE selected window)
    TrialsTable = selectAndComputeTrials_entire(t, x, Fs, nTrials, label);

    % Summary across trials (mean±SD)
    SummaryTable = table(string(label), ...
        mean(TrialsTable.P2P_uV,'omitnan'),  std(TrialsTable.P2P_uV,'omitnan'), ...
        mean(TrialsTable.RMS_uV,'omitnan'),  std(TrialsTable.RMS_uV,'omitnan'), ...
        mean(TrialsTable.MAV_uV,'omitnan'),  std(TrialsTable.MAV_uV,'omitnan'), ...
        mean(TrialsTable.iEMG_uVs,'omitnan'),std(TrialsTable.iEMG_uVs,'omitnan'), ...
        mean(TrialsTable.Duration_s,'omitnan'), std(TrialsTable.Duration_s,'omitnan'), ...
        'VariableNames', {'Session', ...
        'Mean_P2P_uV','SD_P2P_uV', ...
        'Mean_RMS_uV','SD_RMS_uV', ...
        'Mean_MAV_uV','SD_MAV_uV', ...
        'Mean_iEMG_uVs','SD_iEMG_uVs', ...
        'Mean_Dur_s','SD_Dur_s'});

    meta = struct();
    meta.fullPath  = fullPath;
    meta.EMG_field = EMG_field;
    meta.Fs        = Fs;
    meta.doFilter  = usedFilter;
end

function [t, x, Fs, usedFilter] = loadOneEMG(fullPath, EMG_field, doFilter)
    [~, ~, ext] = fileparts(fullPath);
    ext = lower(ext);

    switch ext
        case '.cfs'
            tmp = readCFSfile(fullPath);
        case '.mat'
            S = load(fullPath);
            if isfield(S,'tmp'), tmp = S.tmp; else, tmp = S; end
        otherwise
            error('Unsupported file type: %s', ext);
    end

    if ~isstruct(tmp), error('Loaded tmp is not a struct.'); end

    f = fieldnames(tmp);
    if isfield(tmp,'tmp')
        raw = tmp.tmp;
    elseif numel(f)==1
        raw = tmp.(f{1});
    else
        raw = tmp;
    end

    % Reindex fields
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

    if ~isfield(data, EMG_field)
        emgs = fieldnames(data);
        emgs = emgs(startsWith(emgs,'EMG_'));
        error('EMG "%s" not found. Available: %s', EMG_field, strjoin(emgs, ', '));
    end

    sig = data.(EMG_field).dat(:);
    Fs  = data.(EMG_field).FreqS;
    t   = (0:numel(sig)-1)'/Fs;

    usedFilter = false;
    if doFilter && exist("filtrage","file")==2
        try
            x = filtrage(sig, Fs, 20, 450);
            usedFilter = true;
        catch ME
            % MATLAB error 5006 can occur at butter() because MATLAB's
            % required services are unavailable. Never silently switch the
            % analysis to raw EMG: ask the researcher first.
            msg = sprintf(['20–450 Hz filtering failed:\n\n%s\n\n' ...
                'Use UNFILTERED EMG for BOTH PRE and POST, or cancel?\n' ...
                'The MVC calculations themselves will not change.'], ME.message);
            choice = questdlg(msg, 'MVC filtering unavailable', ...
                'Use raw EMG', 'Cancel', 'Cancel');
            if ~strcmp(choice,'Use raw EMG')
                error('MVC:FilterUnavailable', ...
                    'MVC cancelled: the requested 20–450 Hz filter failed. Original error: %s', ME.message);
            end
            warning('MVC:RawEMG', ...
                'Using unfiltered EMG for PRE and POST because the requested filter failed.');
            x = sig;
        end
    elseif doFilter
        error('MVC:FilterMissing', ...
            'Filtering requested, but filtrage.m was not found on the MATLAB path.');
    else
        x = sig;
    end
end

function TrialsTable = selectAndComputeTrials_entire(t, x, Fs, nTrials, label)
    figure('Color','w','Name',sprintf('%s MVC selection', label));
    plot(t, x); grid on; box off;
    xlabel('Time (s)'); ylabel('EMG (V)');
    title(sprintf('%s: Select %d MVC windows (click START then END)', label, nTrials));
    hold on;

    Trial     = (1:nTrials).';
    Session   = repmat(string(label), nTrials, 1);
    tStart    = nan(nTrials,1);
    tEnd      = nan(nTrials,1);
    Duration  = nan(nTrials,1);
    P2P_uV    = nan(nTrials,1);
    RMS_uV    = nan(nTrials,1);
    MAV_uV    = nan(nTrials,1);
    iEMG_uVs  = nan(nTrials,1);

    for k = 1:nTrials
        fprintf('\n[%s] MVC %d/%d: click START then END...\n', label, k, nTrials);
        [xs, ~] = ginput(1);
        [xe, ~] = ginput(1);

        a = max(min(xs,xe), t(1));
        b = min(max(xs,xe), t(end));

        idx = (t>=a) & (t<=b);
        seg = x(idx);

        if numel(seg) < 2
            warning('[%s] MVC %d too short. Skipping.', label, k);
            continue;
        end

        % Convert to µV for reporting
        seg_uV = seg * 1e6;

        tStart(k)   = a;
        tEnd(k)     = b;
        Duration(k) = (b-a);

        % Metrics on ENTIRE selected window
        P2P_uV(k)   = max(seg_uV) - min(seg_uV);
        RMS_uV(k)   = sqrt(mean(seg_uV.^2));
        MAV_uV(k)   = mean(abs(seg_uV));
        iEMG_uVs(k) = trapz(abs(seg_uV)) / Fs;  % µV·s

        % Visual feedback
        xline(a,'--',sprintf('T%d start',k),'LabelVerticalAlignment','bottom');
        xline(b,'--',sprintf('T%d end',k),'LabelVerticalAlignment','bottom');

        yl = ylim;
        patch([a b b a],[yl(1) yl(1) yl(2) yl(2)],0.9*[1 1 1],'EdgeColor','none','FaceAlpha',0.15);
        ylim(yl);

        fprintf('  RMS=%.1f µV | MAV=%.1f µV | P2P=%.1f µV | iEMG=%.1f µV·s | Dur=%.2f s\n', ...
            RMS_uV(k), MAV_uV(k), P2P_uV(k), iEMG_uVs(k), Duration(k));
    end

    TrialsTable = table(Session, Trial, tStart, tEnd, Duration, P2P_uV, RMS_uV, MAV_uV, iEMG_uVs, ...
        'VariableNames', {'Session','Trial','tStart_s','tEnd_s','Duration_s','P2P_uV','RMS_uV','MAV_uV','iEMG_uVs'});
end