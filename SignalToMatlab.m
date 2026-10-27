%% About this script
%{
    This script lets you select and extract MEPs from EMG data recorded
    with the CED1401 system and Signal software.

    It expects a single .mat file as input. Therefore, the original
    Signal .cfs file must first be converted using Fabien's function:
                                       'readCFSfile.m'
    Note: this function does not work on macOS.

    If exactly one Brainsight .txt export (which may contain stimulations
    distributed across different .mat files) is found alongside the .mat 
    file, neuronavigation errors will be displayed in the MEP selection UI 
    and exported to the final .mat and .csv files.

    * * * * *

    If you run into any issue, please contact me.
    — Mathilde
%}

%% Clearing the environment
clc
clear
close all

addpath("functions");
addpath("matcfs64c");

%% Export .cfs to .mat
% IF ON MACOS, IGNORE THIS SECTION AND RUN THE NEXT ONE

% % 1) Select the .cfs file
% [cfFile, cfPath] = uigetfile({'*.cfs','Signal files (*.cfs)'}, ...
%                               'Select the .cfs file');
% if isequal(cfFile,0)
%     error('No .cfs file selected. Operation cancelled.');
% end
% cfFull = fullfile(cfPath, cfFile);
% 
% % 2) Read the .cfs using Fabien's function
% tmp = readCFSfile(cfFull);
% 
% % 3) Choose the name and location for saving the .mat
% [~, baseName] = fileparts(cfFile);
% defaultMatName = [regexprep(baseName, '\s+', '_') '.mat']; % replace spaces with "_"
% [matFile, matPath] = uiputfile({'*.mat','MAT-file (*.mat)'}, ...
%                                'Save as...', defaultMatName);
% if isequal(matFile,0)
%     error('Save cancelled.');
% end
% matFull = fullfile(matPath, matFile);
% 
% % 4) Save
% save(matFull, 'tmp', '-v7.3');
% fprintf('MAT file saved: %s\n', matFull);

%% Looking for the file & loading the data
% IF ON WINDOWS, YOU CAN DIRECTLY RUN THE ABOVE SECTION BUT RUN THAT ONE
% TOO

[file, file_dir] = uigetfile('*.mat');
str_file = convertCharsToStrings(file);
str_file_dir = convertCharsToStrings(file_dir);
str_file_path = str_file_dir + str_file;
tmp = load(str_file_path);
fprintf('OK — MAT loaded (%s). Fields reindexed.\n', str_file);

%% Selecting the EMG used and cleaning the field names

fields_tmp = fieldnames(tmp);
raw_indexed_data = tmp.(fields_tmp{1});
raw_fields = fieldnames(raw_indexed_data);
emg_indices = find(startsWith(raw_fields, 'EMG'));
nb_EMGs = length(emg_indices);
available_EMGs = raw_fields(emg_indices);
emg_numbers = extractAfter(available_EMGs, 'EMG_');

if nb_EMGs > 1
    msg = sprintf('%d available EMGs. Enter custom numbers to rename them, or use the sensor numbers (e.g., 2 3 4):', nb_EMGs);
    answer = inputdlg(msg, 'EMG Renaming', 1, {strjoin(emg_numbers, ' ')});
    selected_nums = split(strtrim(answer{1}))';
    selected_EMGs = "EMG_" + selected_nums;

    msg_work = sprintf('Which EMG to you want to analyse (%s)?', strjoin(selected_nums, ', '));
    num_work_answer = inputdlg(msg_work, 'EMG Choice', 1, selected_nums(1));
    EMG_field = "EMG_" + strtrim(num_work_answer{1});
else
    EMG_field = available_EMGs{1};
    selected_EMGs = EMG_field;
end

% Reindexing the structure to use labels easier
metaFieldNames = {'Acquisition_time','Acquisition_date','comment','DataSections','FrameStart_s'};
metaFieldNames = metaFieldNames(ismember(metaFieldNames, raw_fields));

data = struct();
for i = 1:numel(metaFieldNames)
    data.(metaFieldNames{i}) = raw_indexed_data.(metaFieldNames{i});
end

channelFieldNames = setdiff(raw_fields, [metaFieldNames(:); available_EMGs(:)], 'stable');
for i = 1:numel(channelFieldNames)
    oldName = channelFieldNames{i};
    newName = oldName(1:end-2);   % removes the last 2 chars ('_X')
    data.(newName) = raw_indexed_data.(oldName);
end
if nb_EMGs > 1
    for i = 1:length(selected_EMGs)
        data.(selected_EMGs{i}) = raw_indexed_data.(available_EMGs{i});
    end
elseif nb_EMGs == 1
    data.(EMG_field) = raw_indexed_data.(EMG_field);
end

%% Get signals : EMG / Stim

% Get the EMG signal and filter it
% TODO: there may be multiple EMG channels
%       => decide how to select the appropriate channel

if isfield(data, EMG_field)
    EMG = data.(EMG_field).dat;
    freq_EMG = data.(EMG_field).FreqS;
end

% Using Silvère's filtering function
EMG_filtered = filtrage(EMG, freq_EMG, 20, 1000);

% Get the stimulation signal considering
% the signal was acquired on the ADC0
stim = data.ADC0.dat;
freq_stim = data.ADC0.FreqS;

fprintf('OK — EMG filtered (20–1000 Hz, Fs=%.1f Hz). Stim loaded (Fs=%.1f Hz).\n', freq_EMG, freq_stim);

%% Detect stimulation times and match sampling rates
% Resamples per frame (not globally) to avoid accumulating drift caused by the difference in sampling frequencies

if isfield(data.ADC0,'PointsPerFrame')
    stimPointsPerFrame = data.ADC0.PointsPerFrame;
elseif isfield(data,'DataSections')
    nSectionsFallback = data.DataSections;
    stimPointsPerFrame = repmat(length(stim)/nSectionsFallback, nSectionsFallback, 1);
end

emgPointsPerFrame = data.(EMG_field).PointsPerFrame;
emgFrameBoundaries = cumsum(emgPointsPerFrame);
emgFrameStartIdx = [0;emgFrameBoundaries(1:end-1)]+1;
stimFrameBoundaries = cumsum(stimPointsPerFrame);
stimFrameStartIdx = [0;stimFrameBoundaries(1:end-1)]+1;

nSections = numel(emgPointsPerFrame);
Thr = 0.1*max(stim);                                                        % threshold: stim signal above 10% of max voltage

listOfStim = [];
stimFrameIdx = [];
for f = 1:nSections
    emgLocalLen = emgPointsPerFrame(f);
    stimLocal = stim(stimFrameStartIdx(f):stimFrameBoundaries(f));

    time_stim_local = (0:numel(stimLocal)-1) / freq_stim;                   % local time base of this section's stim
    new_time_stim_local = (0:emgLocalLen-1) / freq_EMG;                     % local time base of this section's EMG

    new_stim_local = interp1(time_stim_local, stimLocal, new_time_stim_local, 'linear');

    logic = new_stim_local > Thr;                                           % values above threshold
    localStims = find(diff(logic)==1) + 1;                                  % rising edges, local to this section

    if ~isempty(localStims)
        globalStims = localStims + emgFrameStartIdx(f) - 1;                 % back to global EMG sample index
        listOfStim = [listOfStim, globalStims]; %#ok<AGROW>
        stimFrameIdx = [stimFrameIdx, repmat(f, 1, numel(globalStims))]; %#ok<AGROW>
    end
end

fprintf('OK — Stim detection completed (per-section resampling, %d sections, %d stims detected).\n', ...
    nSections, numel(listOfStim));

%% Match Brainsight neuronavigation errors

% Runs automatically if exactly one *.txt is found next to the *.mat file

BrainsightTable = table();
matchedRows = nan(numel(listOfStim),1);

hasFrameTiming = isfield(data,'Acquisition_time') && isfield(data,'FrameStart_s') ...
    && ~isempty(data.FrameStart_s) && ~isempty(stimFrameIdx);

bsCandidates = dir(fullfile(char(str_file_dir), '*.txt'));
if ~hasFrameTiming
    fprintf('Brainsight matching skipped (no CFS frame timing available for this .mat).\n');
elseif isempty(bsCandidates)
    fprintf('Brainsight matching skipped (no .txt found in %s).\n', str_file_dir);
elseif numel(bsCandidates) > 1
    warning('Brainsight matching skipped — multiple .txt files found in %s: %s', ...
        char(str_file_dir), strjoin({bsCandidates.name}, ', '));
else
    BrainsightTable = readBrainsightSamples(fullfile(bsCandidates(1).folder, bsCandidates(1).name));
    [matchedRows, MatchLog] = matchBrainsightErrors( ...
        data.Acquisition_time, data.FrameStart_s, stimFrameIdx, BrainsightTable);
    fprintf('OK — Brainsight errors matched (%d/%d stims, file: %s).\n', ...
        MatchLog.nMatched, MatchLog.nTotal, bsCandidates(1).name);
end

% Stim errors (one per column of allMEP), for display in selectingMEP.
BrainsightErrors = struct( ...
    'TargetError_mm',   nan(numel(listOfStim),1), ...
    'AngularError_deg', nan(numel(listOfStim),1), ...
    'TwistError_deg',   nan(numel(listOfStim),1));
validMatch = ~isnan(matchedRows);
if any(validMatch)
    BrainsightErrors.TargetError_mm(validMatch)   = BrainsightTable.TargetError(matchedRows(validMatch));
    BrainsightErrors.AngularError_deg(validMatch) = BrainsightTable.AngularError(matchedRows(validMatch));
    BrainsightErrors.TwistError_deg(validMatch)   = BrainsightTable.TwistError(matchedRows(validMatch));
end

%% Double footswitch: optional gait-phase separation

% If 3 ADC channels are present (ADC0 = stim, ADC1/ADC2 = two footswitches)
% this section offers to split the MEPs by gait-cycle phase.

dataFields = fieldnames(data);
adcFields = dataFields(startsWith(dataFields, 'ADC'));
hasDoubleFS = numel(adcFields) == 3 && all(ismember({'ADC0','ADC1','ADC2'}, adcFields));

phaseChoice = "";

if hasDoubleFS
    choice = questdlg( ...
        sprintf('%d ADC channels detected (stim + 2 footswitches). Split the MEPs by gait-cycle phase ?', numel(adcFields)), ...
        'Double Footswitch', 'Yes', 'No', 'No');

    if strcmp(choice, 'Yes')
        freqFS = data.ADC1.FreqS;
        [GaitPhase, FSInfo] = identifyGaitPhase(data.ADC1.dat, data.ADC2.dat, freqFS, listOfStim, freq_EMG);

        fprintf('OK — Gait phase identified for %d/%d stims (%d unresolved; %d alternation conflicts resolved automatically).\n', ...
            sum(GaitPhase ~= "Unknown"), numel(listOfStim), sum(GaitPhase == "Unknown"), FSInfo.nAltConflict);
        for ch = ["FS1","FS2"]
            n = FSInfo.frontVotes.(ch) + FSInfo.backVotes.(ch) + FSInfo.unknownVotes.(ch);
            fprintf('   %s triggered %d stims — front criterion: %d/%d, back criterion: %d/%d  =>  assigned as %s\n', ...
                ch, n, FSInfo.frontVotes.(ch), n, FSInfo.backVotes.(ch), n, FSInfo.position.(ch));
        end

        phaseAnswer = questdlg( ...
            'Which gait-cycle phase do you want to analyse?', ...
            'Gait Phase Selection', 'Heel-off (~40%)', 'Toe-off (~60%)', 'Heel-off (~40%)');

        if strcmp(phaseAnswer, 'Toe-off (~60%)')
            phaseChoice = "60pct";
        else
            phaseChoice = "40pct";
        end

        keepMask = (GaitPhase == phaseChoice);
        fprintf('OK — Keeping %d/%d stims for phase %s (%d excluded: other phase or unresolved).\n', ...
            sum(keepMask), numel(listOfStim), phaseChoice, numel(listOfStim) - sum(keepMask));

        listOfStim   = listOfStim(keepMask);
        stimFrameIdx = stimFrameIdx(keepMask);

        % Keep the Brainsight data aligned with the filtered stims.
        matchedRows = matchedRows(keepMask);
        BrainsightErrors.TargetError_mm   = BrainsightErrors.TargetError_mm(keepMask);
        BrainsightErrors.AngularError_deg = BrainsightErrors.AngularError_deg(keepMask);
        BrainsightErrors.TwistError_deg   = BrainsightErrors.TwistError_deg(keepMask);
        validMatch = ~isnan(matchedRows);
    end
end

%% Build MEP windows
% Define the window where the MEP should appear (stim -100 ms, stim +500 ms)

MEPWindows = [];
for t = 1:length(listOfStim)
    minus = round(listOfStim(t) - 0.1 * freq_EMG) ; % time of stim - 100ms
    plus = round(listOfStim(t) + 0.5 * freq_EMG) ;  % time of stim + 500ms
    if minus < 1    % ensure the lower bound is within data
        minus = 1;
    end
    if plus > length(EMG_filtered)  % ensure the upper bound is within data
        plus = length(EMG_filtered);
    end

    wdw = [minus, plus];            % time indexes of the window around the stim
    MEPWindows = [MEPWindows; wdw]; %#ok<AGROW>
end
fprintf('OK — MEP windows created.\n')

%% Extract all MEP segments
% Longueur cible (en nombre d'échantillons)
winLen = MEPWindows(1,2) - MEPWindows(1,1) + 1;  % ou max(...) si ça varie

nWin   = size(MEPWindows,1);
allMEP = nan(winLen, nWin);   % pré-allocation avec NaN

for w = 1:nWin
    fstart = MEPWindows(w,1);
    fend   = MEPWindows(w,2);

    EMG_window = EMG_filtered(fstart:fend);
    L = numel(EMG_window);

    % Adding NaN values if window is smaller than the other
    allMEP(1:L, w) = EMG_window;
end

% Create a time vector for plotting (in ms, aligned with the window definition)
time =  linspace(-100, 500, size(allMEP,1));

% Select valid MEPs (manual/GUI function)
[selectedMEPs, selectedIdx] = selectingMEP(allMEP, time, BrainsightErrors);

%% MEPs structure
% Create MEP struct (keep valid MEPs and rename to MEP_01, MEP_02, ... original naming is reported too)
% and MEP_SELECTION struct (report which MEPs were rejected, which were kept,
% and the total number initially detected)

[MEP, MEP_SELECTION] = renumberLogMEP(selectedMEPs, selectedIdx, allMEP);

% Time centering aroung 0 ms = stim index in MEP window
stimIdx0 = round(0.1 * freq_EMG); % sample for 100 ms (0-based)

% Store at global level (only valid)
MEP.Meta.Time_ms = time(:).';
MEP.Meta.StimIdx = stimIdx0 + 1;
MEP.Meta.Fs      = freq_EMG;

% CFS frame timing for Brainsight error matching.
if isfield(data,'Acquisition_time')
    MEP.Meta.Acquisition_time = data.Acquisition_time;
end
if isfield(data,'Acquisition_date')
    MEP.Meta.Acquisition_date = data.Acquisition_date;
end
if isfield(data,'FrameStart_s')
    MEP.Meta.FrameStart_s = data.FrameStart_s;
end
MEP.Meta.StimFrameIdx = stimFrameIdx;
if phaseChoice ~= ""
    MEP.Meta.GaitPhase = char(phaseChoice);
end

% Create a struct with all individual MEPs using original naming
originalNamedMEPs = namingMEP(selectedMEPs, selectedIdx);   % creates a struct,
                                                    % if needed later
fprintf('OK — MEP struct created and renumbered (MEP_01..MEP_%02d). Selection log stored.\n', size(selectedMEPs,2));



%% Analyse of the MEPS

% Detect valid MEP windows + peak-to-peak + latency (automatic)

[MEP, T] = detectMEPOnsetOffset(MEP, 'Fs', freq_EMG);
fprintf('OK — Onset/offset, peak-to-peak (p2p), latency, and AUC extracted automatically.\n');

% Attach the Brainsight errors to the MEP field.
if any(validMatch)
    MEP = attachBrainsightErrors(MEP, matchedRows, BrainsightTable);
end

% Attach the EMG channels that were not analysed (MEP.MEP_XX.EMG_add.EMG_xx).
if nb_EMGs > 1
    MEP = attachAdditionalEMGs(MEP, data, selected_EMGs, EMG_field, MEPWindows, freq_EMG);
end

%% Structure export

[~, baseMatName] = fileparts(char(str_file));  % get .mat file name without extension
if phaseChoice ~= ""
    baseMatName = sprintf('%s_%s', baseMatName, phaseChoice);
end
default=fullfile(char(str_file_dir), sprintf('%s_MEPs',baseMatName));
[matFile, matPath] = uiputfile({'*.mat'},'Save MEP structure as :', default);
if isequal(matFile,0)
    warning('Export canceled.');
else
    out = fullfile(matPath, matFile);
    save(out, "MEP");
    fprintf('MEP structure exported: %s\n', out);
end

%% === CSV export for statistical analysis: 1 row per MEP; columns = P2P, Latency, AUC ===

% 1) Build the table to export (keep also the MEP label)
% Convert T.Signal to a string to export it to csv after:
T.SignalString = cellfun(@(x) sprintf('%g,', x), T.Signal, 'UniformOutput', false);
T.SignalString = cellfun(@(s) s(1:end-1), T.SignalString, 'UniformOutput', false);

nMEP = height(T);
TargetError_mm   = nan(nMEP,1);
AngularError_deg = nan(nMEP,1);
TwistError_deg   = nan(nMEP,1);
Distance_mm      = nan(nMEP,1);
for i = 1:nMEP
    lbl = T.Label{i};
    if isfield(MEP,lbl) && isfield(MEP.(lbl),'Brainsight')
        TargetError_mm(i)   = MEP.(lbl).Brainsight.TargetError_mm;
        AngularError_deg(i) = MEP.(lbl).Brainsight.AngularError_deg;
        TwistError_deg(i)   = MEP.(lbl).Brainsight.TwistError_deg;
        Distance_mm(i)      = MEP.(lbl).Brainsight.Distance_mm;
    end
end

ExportTab = table( ...
    T.Label, ...
    T.P2P_uV, ...
    T.Latency_ms, ...
    T.AUC_uVms, ...
    T.SPduration_ms, ...
    T.SignalString, ...
    TargetError_mm, ...
    AngularError_deg, ...
    TwistError_deg, ...
    Distance_mm, ...
    'VariableNames', {'MEP_Label','P2P_uV','Latency_ms','AUC_uVms','SP_ms','Raw_signal', ...
                       'TargetError_mm','AngularError_deg','TwistError_deg','Distance_mm'});

% 2) Propose a default file name (same folder as the .mat)
defaultCSV = fullfile(char(str_file_dir), sprintf('%s_MEP_metrics.csv', baseMatName));

% 3) Save location
[csvFile, csvPath] = uiputfile({'*.csv','CSV file (*.csv)'}, 'Save MEP metrics as...', defaultCSV);
if isequal(csvFile,0)
    warning('CSV export canceled by user.');
else
    outCSV = fullfile(csvPath, csvFile);
    writetable(ExportTab, outCSV);
    fprintf('CSV exported: %s (N=%d MEPs)\n', outCSV, height(ExportTab));
end
