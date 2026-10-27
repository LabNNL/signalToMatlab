function [MEP, summary] = detectMEPOnsetOffset(MEP, varargin) 
% detectMEPOnsetOffset
% Detects onset/offset (and latency, duration, P2P) for each MEP_* in the MEP struct.
%
% USAGE (simple):
%   [MEP, T] = detectMEPOnsetOffset(MEP);
%
% USAGE (with overrides):
%   P = struct('rms_ms',5,'base_ms',[-100 -20],'search_ms',[10 100], ...
%              'k_on',2.5,'k_off',1.5,'hold_on_ms',5,'hold_off_ms',5);
%   [MEP, T] = detectMEPOnsetOffset(MEP, 'Params', P, 'Fs', Fs, 'Time', time);
%
% Inputs (Name-Value):
%   'Params'   : struct of parameters (see defaults below)
%   'Fs'       : sampling frequency if MEP.Meta.Fs is missing
%   'Time'     : time vector (ms) if MEP.Meta.Time_ms is missing (size = nSamp)
%   'EMGPrefix': prefix of EMG fields (default: 'EMG')
%
% Outputs:
%   MEP     : enriched MEP struct (OnOff_ms, OnOff_idx, Baseline, Thresholds, etc.)
%   summary : summary table for each MEP

% ---------- Default parameters ----------
Pdef = struct( ...
    'rms_ms',      5, ...           % RMS window length (ms)
    'base_ms',     [-100 -20], ...  % baseline window (pre-stim, ms)
    'search_ms',   [10 100], ...    % search window for onset/offset (ms)
    'k_on',        2, ...           % onset threshold = mean + k_on * SD
    'k_off',       1.5, ...         % offset threshold (hysteresis)
    'hold_on_ms',  5, ...           % minimum duration above threshold (ms)
    'hold_off_ms', 5 ...            % minimum duration below threshold (ms)
);
EMGPrefix = 'EMG';
Fs = [];
time = [];

% ---------- Parse Name-Value pairs ----------
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'params'
            Pdef = overrideStruct(Pdef, varargin{i+1});
        case 'fs'
            Fs = varargin{i+1};
        case 'time'
            time = varargin{i+1};
        case 'emgprefix'
            EMGPrefix = varargin{i+1};
        otherwise
            error('Unknown parameter: %s', varargin{i});
    end
end
P = Pdef;

% Infer the length of a MEP from the first MEP_* field
names = fieldnames(MEP);
names = names(startsWith(names, 'MEP_'));
if isempty(names)
    warning('No MEP_* field found. Nothing to process.');
    summary = table();
    return;
end
% Find the first EMG_* field to determine nSamp
firstLab = names{1};
flds1 = fieldnames(MEP.(firstLab));
iEmg = find(startsWith(flds1, EMGPrefix), 1, 'first');
if isempty(iEmg)
    error('No %s* field found in %s.', EMGPrefix, firstLab);
end
nSamp = numel(MEP.(firstLab).(flds1{iEmg})); % nb of points

if isempty(time)
    % If no time vector provided, try MEP.Meta.Time_ms
    if isfield(MEP, 'Meta') && isfield(MEP.Meta, 'Time_ms') ...
            && numel(MEP.Meta.Time_ms) == nSamp
        time = MEP.Meta.Time_ms(:);
    else
        % Default: generate evenly spaced time covering [-100, +500] ms
        win_ms = [-100, 500];
        time = linspace(win_ms(1), win_ms(2), nSamp).';
    end
else
    time = time(:);
    if numel(time) ~= nSamp
        error('"Time" length (%d) does not match nSamp (%d).', numel(time), nSamp);
    end
end

% Index corresponding to stimulation time (0 ms)
stimIdx0 = find(time >= 0, 1, 'first');
if isempty(stimIdx0)
    % Fallback: ~100 ms after beginning (assuming window = [-100, +500])
    stimIdx0 = round(0.1 * Fs) + 1;
end

% ---------- Pre-computations related to Fs ----------
N_rms      = max(1, round(P.rms_ms      * Fs / 1000)); % nb of points for RMS window
N_hold_on  = max(1, round(P.hold_on_ms  * Fs / 1000)); % nb of consecutive points > threshold_on
N_hold_off = max(1, round(P.hold_off_ms * Fs / 1000)); % IDEM for threshold_off

bIdx = (time >= P.base_ms(1))  & (time <= P.base_ms(2)); % baseline indices
sIdx = (time >= P.search_ms(1))& (time <= P.search_ms(2)); % search indices
assert(any(bIdx) && any(sIdx), ...
    'Invalid baseline/search windows. Check P.base_ms and P.search_ms.');

% ---------- Loop over each MEP ----------
OnOff_all = nan(numel(names), 2);
Summ = struct('Label',{},'On_ms',{},'Off_ms',{},'Latency_ms',{}, ...
              'Duration_ms',{},'P2P_uV',{}, 'AUC_uVms',{},'SPduration_ms',{}, 'BaseMean',{},'BaseSD',{}, ...
              'ThrOn',{},'ThrOff',{},'ThrSP',{});

for k = 1:numel(names)
    lab = names{k}; % ex: 'MEP_01'
    flds = fieldnames(MEP.(lab));
    emgField = '';
    for f = 1:numel(flds)
        if startsWith(flds{f}, EMGPrefix)
            emgField = flds{f}; break
        end
    end
    if isempty(emgField)
        warning('No %s* field in %s. Skipping.', EMGPrefix, lab);
        continue
    end

    sig = MEP.(lab).(emgField)(:); % raw EMG signal

    % 1) RMS envelope
    env = sqrt(movmean(sig.^2, N_rms));

    % 2) Baseline & thresholds
    mu = mean(env(bIdx));
    sd = std(env(bIdx));
    thr_on  = mu + P.k_on  * sd;
    thr_off = mu + P.k_off * sd;
    thr_sp  = mu; 

    % 3) Onset detection (choix du plus gros burst au-dessus du seuil)
    env_s = env(sIdx);                         % enveloppe dans la fenêtre de recherche
    s_ofs = find(sIdx, 1, 'first') - 1;        % offset pour revenir aux indices globaux

    logic = env_s > thr_on;                    % points au-dessus du seuil onset

    if any(logic)
        % Détection de tous les runs consécutifs au-dessus du seuil
        d = diff([false; logic(:); false]);
        starts_all = find(d == 1);
        ends_all   = find(d == -1) - 1;
        lens       = ends_all - starts_all + 1;

        % On garde seulement les runs qui durent au moins N_hold_on échantillons
        valid = lens >= N_hold_on;
        starts = starts_all(valid);
        ends   = ends_all(valid);

        if isempty(starts)
            idx_on = NaN;
        else
            % Pour chaque run, on calcule le pic de l'enveloppe RMS
            peakVals = arrayfun(@(a,b) max(env_s(a:b)), starts, ends);

            % On choisit le run ayant le pic maximal = "vrai" MEP
            [~, iBest] = max(peakVals);
            on_rel = starts(iBest);

            % Conversion en index global
            idx_on = s_ofs + on_rel;
        end
    else
        idx_on = NaN;
    end

    % 4) Offset detection
    idx_off = NaN;
    if ~isnan(idx_on)
        env_aft = env(idx_on:end);
        off_rel = firstRun(env_aft < thr_off, N_hold_off);
        if ~isempty(off_rel)
            idx_off = idx_on + off_rel - 1;
        end
    end

    % 5) Silent period detection
    idx_sp = NaN;
    if ~isnan(idx_off)
        decal=100;                                                          % Right offset to avoid values above threshold post-MEP
        env_sp = env(idx_off+decal:end);                                    % First index post-offset below the sp threshold
        sp_rel = firstRun(env_sp > thr_sp, 0);                              % First index above threshold
        idx_sp = idx_off + decal + sp_rel -1;
        if isempty(idx_sp)
            idx_sp = NaN;
        end
    end

    % 6) Save results (ms / idx / metrics)
    on_ms  = NaN; off_ms = NaN; sp_ms = NaN;
    if ~isnan(idx_on),  on_ms  = time(idx_on);  end
    if ~isnan(idx_off), off_ms = time(idx_off); end
    if ~isnan(idx_sp), sp_ms = time(idx_sp)-off_ms; end

    MEP.(lab).OnOff_ms       = [on_ms, off_ms];
    MEP.(lab).OnOff_idx      = [idx_on, idx_off];
    MEP.(lab).Baseline.Mean  = mu;
    MEP.(lab).Baseline.SD    = sd;
    MEP.(lab).Thresholds.on  = thr_on;
    MEP.(lab).Thresholds.off = thr_off;
    MEP.(lab).Thresholds.sp  = thr_sp;
    MEP.(lab).Enveloppe      = env;
    MEP.(lab).Silentperiod   = sp_ms;

    OnOff_all(k,:) = [on_ms, off_ms];

    % Summary
    Summ(k).Label       = lab;
    Summ(k).BaseMean    = mu;
    Summ(k).BaseSD      = sd;
    Summ(k).ThrOn       = thr_on;
    Summ(k).ThrOff      = thr_off;
    Summ(k).ThrSP       = thr_sp;
    Summ(k).Signal      = MEP.(lab).EMG;

end

% ---------- Save Meta & Summary ----------
MEP.Meta.Params    = P;
MEP.Meta.Fs        = Fs;
MEP.Meta.Time_ms   = time(:).';
MEP.Meta.OnOff_ms  = OnOff_all;
MEP.Meta.MEP_Order = string(names);

% Onset-Offset inspection and manual modification

MEP=OnsetOffset_Inspection(MEP);

for k=1:numel(names)
    
    lab = names{k};
    on_ms=MEP.(lab).OnOff_ms(1,1);
    off_ms=MEP.(lab).OnOff_ms(1,2);
    idx_on=MEP.(lab).OnOff_idx(1,1);
    idx_off=MEP.(lab).OnOff_idx(1,2);
    SPD=MEP.(lab).Silentperiod;
    sig=MEP.(lab).EMG;

    if ~isnan(idx_on) && ~isnan(idx_off) && idx_off > idx_on
        seg = sig(idx_on:idx_off);
        duration_ms  = off_ms - on_ms;
        p2p_uV       = (max(seg)-min(seg)) * 1e6;   % if 'sig' in volts
        AUC_uVms = (trapz(abs(seg)) / Fs) * 1e9;   % µV·ms
    else
        duration_ms  = NaN;
        p2p_uV       = NaN;
        AUC_uVms     = NaN;
    end

    MEP.(lab).Latency_ms    = on_ms;
    MEP.(lab).Duration_ms   = duration_ms;
    MEP.(lab).Peak2Peak_uV  = p2p_uV;
    MEP.(lab).AUC_uVms      = AUC_uVms;
    MEP.(lab).SPduration_ms = SPD;

    Summ(k).On_ms=MEP.Meta.OnOff_ms(k,1);
    Summ(k).Off_ms=MEP.Meta.OnOff_ms(k,2);
    Summ(k).Latency_ms    = on_ms;
    Summ(k).Duration_ms   = duration_ms;
    Summ(k).P2P_uV        = p2p_uV;
    Summ(k).AUC_uVms      = AUC_uVms;
    Summ(k).SPduration_ms = SPD;
end

% ---------- Exclude MEPs unchecked in OnsetOffset_Inspection ----------
included = true(numel(names),1);
for k = 1:numel(names)
    if isfield(MEP.(names{k}), 'Included')
        included(k) = MEP.(names{k}).Included;
    end
end

excludedNames = names(~included);
for k = 1:numel(excludedNames)
    MEP = rmfield(MEP, excludedNames{k});
end

MEP.Meta.OnOff_ms  = MEP.Meta.OnOff_ms(included,:);
MEP.Meta.MEP_Order = MEP.Meta.MEP_Order(included);
Summ = Summ(included);

if isfield(MEP,'All')
    MEP.All = MEP.All(included,:);
end

keptNames = names(included);
for k = 1:numel(keptNames)
    if isfield(MEP.(keptNames{k}), 'Included')
        MEP.(keptNames{k}) = rmfield(MEP.(keptNames{k}), 'Included');
    end
end

summary = struct2table(Summ);

% ---------- Subfunctions ----------
function out = overrideStruct(base, over)
    % Override struct fields in 'base' with values in 'over'
    out = base;
    if ~isstruct(over), return; end
    fn = fieldnames(over);
    for ii = 1:numel(fn)
        out.(fn{ii}) = over.(fn{ii});
    end
end

function idx = firstRun(logicalVec, runLen)
% Returns the index of the beginning of the first run of length >= runLen
% in logicalVec (logical vector). Empty if nothing found.
    if ~islogical(logicalVec), logicalVec = logical(logicalVec); end
    if ~any(logicalVec)
        idx = [];
        return
    end
    % Detect runs
    d = diff([false; logicalVec(:); false]);
    starts = find(d==1);
    ends   = find(d==-1) - 1;
    lens   = ends - starts + 1;
    i = find(lens >= runLen, 1, 'first');
    if isempty(i)
        idx = [];
    else
        idx = starts(i);
    end
end

end