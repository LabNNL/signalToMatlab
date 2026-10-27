function [matchedRows, MatchLog] = matchBrainsightErrors(Acquisition_time, FrameStart_s, StimFrameIdx, BrainsightTable, varargin)

% Matches the acquisition time of each detected Signal stim to a sample from the Brainsight session .txt export.

% Every contiguous window of Brainsight samples (of the same length as the 
% number of detected stims) is scored by how well its inter-sample interval
% pattern matches the Signal stims' interval pattern.
%
% USAGE :
%   [matchedRows, MatchLog] = matchBrainsightErrors(Acquisition_time, FrameStart_s, StimFrameIdx, BrainsightTable)
%   [matchedRows, MatchLog] = matchBrainsightErrors(..., 'Tolerance_s',2, ...)
%
% INPUTS :
%   Acquisition_time : "HH:MM:SS" string from readCFSfile
%   FrameStart_s     : per-frame elapsed time since Acquisition_time
%   StimFrameIdx     : for each DETECTED stim, the CFS frame it originated from
%   BrainsightTable  : from readBrainsightSamples
%
% Name-Value parameters:
%   'AnchorWindow_s'   (1200) max plausible clock offset between the two
%                              computers, used to prune implausible matches
%   'Tolerance_s'      (2)    max allowed residual (s) for a match to be
%                              accepted
%   'ConfidenceMargin' (2)    best candidate's residual must be at least this
%                              many times better than the runner-up (N>1)
%
% OUTPUTS :
%   matchedRows : one row index into BrainsightTable per detected stim (same
%                 order/length as StimFrameIdx), or all-NaN if no confident
%                 match could be found for this acquisition.
%   MatchLog    : struct with nTotal, nMatched, offset_s, residual_s, warnings

p = inputParser;
addParameter(p,'AnchorWindow_s',1200);
addParameter(p,'Tolerance_s',2);
addParameter(p,'ConfidenceMargin',2);
parse(p,varargin{:});
opt = p.Results;

N = numel(StimFrameIdx);
MatchLog = struct('nTotal',N,'nMatched',0,'offset_s',NaN,'residual_s',NaN,'warnings',{{}});

hasTiming = ~isempty(Acquisition_time) && ~isempty(FrameStart_s) && ~isempty(StimFrameIdx);
if ~hasTiming
    msg = 'No CFS frame timing available; Brainsight matching skipped for all stims.';
    warning('matchBrainsightErrors:noTiming','%s',msg);
    MatchLog.warnings{end+1} = msg;
    matchedRows = nan(N,1);
    return
end

t0_s = seconds(timeofday(datetime(Acquisition_time,'InputFormat','HH:mm:ss')));
frameTimes_s = t0_s + FrameStart_s(:);
stimTimes_s  = frameTimes_s(StimFrameIdx(:));                               % one per detected stim, temporal order
bsTimes_s    = seconds(BrainsightTable.Time);

cfsIntervals = diff(stimTimes_s);
nBS = numel(bsTimes_s);

% Search the contiguous Brainsight window that best matches the Signal stim-interval pattern.
best.residual = Inf; best.pos = []; best.offset = NaN;
secondResidual = Inf;
for s = 1:(nBS-N+1)
    window = bsTimes_s(s:s+N-1);
    offset = median(stimTimes_s - window);
    if abs(offset) > opt.AnchorWindow_s
        continue
    end
    if N > 1
        residual = sqrt(mean((cfsIntervals - diff(window)).^2));
    else
        residual = abs(offset);                                             % no pattern to compare; rely on AnchorWindow_s only
    end
    if residual < best.residual
        secondResidual = best.residual;
        best.residual = residual;
        best.pos = s;
        best.offset = offset;
    elseif residual < secondResidual
        secondResidual = residual;
    end
end

confident = ~isempty(best.pos) && best.residual <= opt.Tolerance_s ...
    && (~isfinite(secondResidual) || best.residual <= secondResidual/opt.ConfidenceMargin);

if ~confident
    msg = sprintf(['No confident Brainsight match found (best residual = %.2fs, N=%d stims). ' ...
    'Leaving Brainsight fields as NaN for all stims in this file.'], best.residual, N);
    warning('matchBrainsightErrors:noConfidentMatch','%s',msg);
    MatchLog.warnings{end+1} = msg;
    matchedRows = nan(N,1);
    return
end

matchedRows = (best.pos + (0:N-1)).';
MatchLog.offset_s = best.offset;
MatchLog.residual_s = best.residual;
MatchLog.nMatched = N;

end
