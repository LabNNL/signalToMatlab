function [GaitPhase, FSInfo] = identifyGaitPhase(FS1, FS2, freqFS, listOfStim, freq_EMG)

% Identifie pour chaque stim le footswitch (avant/arrière du pied) ayant 
% declenché la stimulation, puis en déduit la phase du cycle de marche.
%
% Méthode : recherche du dernier front descendant de FS1/FS2 précédant
% chaque stim (avec un filtre médian pour supprimer le bruit et l'artéfact
% électrique sub-ms généré par la TMS). En cas d'ambiguité entre les deux 
% channels, l'alternance imposée par le sequenceur fait foi. La position 
% avant/arrière de chaque canal est ensuite déduite par un vote sur 
% l'ensemble des stims : au moment du front d'un FS, l'autre FS est déjà à
% zéro (FS avant) ou bien encore actif (FS talon).
%
% Outputs :
% GaitPhase : string array (1 par stim) : "40pct", "60pct" ou "Unknown"
% FSInfo : struct de diagnostic (déclencheur par stim, votes avant/arrière,
% position globale assignée a FS1/FS2, nb de conflits d'alternance)

nStim = numel(listOfStim);

medOrder = 9;     % ~1 ms à 10 kHz : supprime l'artéfact TMS sans altérer les vrais fronts
searchWin_s = 0.5;   % recherche du front sur les 500 ms précédant la stim
padSamples = 50;    % marge post-stim pour éviter l'effet de bord de medfilt1 sur la stim

mx = max(max(FS1), max(FS2));
thrHigh = 0.6 * mx;
thrLow = 0.3 * mx;

triggerFS = strings(nStim, 1);
otherState = strings(nStim, 1);

lastTrig = "";
nAltConflict = 0;

for k = 1:nStim
    tSec = (listOfStim(k) - 1) / freq_EMG;
    stimIdxFS = round(tSec * freqFS) + 1;

    winStart = max(1, stimIdxFS - round(searchWin_s * freqFS));
    winEnd   = min(numel(FS1), stimIdxFS + padSamples);

    seg1 = medfilt1(FS1(winStart:winEnd), medOrder);
    seg2 = medfilt1(FS2(winStart:winEnd), medOrder);

    localStimIdx = stimIdxFS - winStart + 1;
    seg1 = seg1(1:localStimIdx);
    seg2 = seg2(1:localStimIdx);

    e1 = lastDescEdge(seg1, thrHigh, thrLow);
    e2 = lastDescEdge(seg2, thrHigh, thrLow);

    cand = "";
    if ~isempty(e1) && isempty(e2)
        cand = "FS1";
    elseif isempty(e1) && ~isempty(e2)
        cand = "FS2";
    elseif ~isempty(e1) && ~isempty(e2)
        if e1 >= e2, cand = "FS1"; else, cand = "FS2"; end
    end

    expected = "";
    if lastTrig == "FS1"
        expected = "FS2";
    elseif lastTrig == "FS2"
        expected = "FS1";
    end

    if cand ~= "" && expected ~= "" && cand ~= expected                     % Désaccord entre le front détecté et l'alternance attendue
        nAltConflict = nAltConflict + 1;
        trig = expected;                                                    % l'alternance fait foi.
    elseif cand ~= ""
        trig = cand;
    else
        trig = expected;
    end

    if trig == ""
        triggerFS(k)  = "Unknown";
        otherState(k) = "unknown";
        lastTrig = "";
        continue
    end

    triggerFS(k) = trig;
    lastTrig = trig;

    if trig == "FS1"
        e = e1; other = seg2;
    else
        e = e2; other = seg1;
    end

    if ~isempty(e)
        otherState(k) = stateAt(other, e, thrHigh, thrLow);
    else
        otherState(k) = "unknown";
    end
end

% Vote avant/arriere : au moment du front d'un FS, l'autre FS est déjà à
% zéro (FS avant) ou bien encore actif (FS arrière).
frontVotes   = struct('FS1', 0, 'FS2', 0);
backVotes    = struct('FS1', 0, 'FS2', 0);
unknownVotes = struct('FS1', 0, 'FS2', 0);
for k = 1:nStim
    if triggerFS(k) == "Unknown"
        continue
    end
    ch = char(triggerFS(k));
    switch otherState(k)
        case "low"
            frontVotes.(ch) = frontVotes.(ch) + 1;
        case "high"
            backVotes.(ch) = backVotes.(ch) + 1;
        otherwise
            unknownVotes.(ch) = unknownVotes.(ch) + 1;
    end
end

score1 = frontVotes.FS1 - backVotes.FS1;
score2 = frontVotes.FS2 - backVotes.FS2;
position = struct();
if score1 >= score2
    position.FS1 = "front"; position.FS2 = "back";
else
    position.FS1 = "back";  position.FS2 = "front";
end

GaitPhase = strings(nStim, 1);
for k = 1:nStim
    if triggerFS(k) == "Unknown"
        GaitPhase(k) = "Unknown";
    elseif position.(char(triggerFS(k))) == "front"
        GaitPhase(k) = "60pct";
    else
        GaitPhase(k) = "40pct";
    end
end

FSInfo = struct( ...
    'triggerFS',    triggerFS, ...
    'otherState',   otherState, ...
    'nAltConflict', nAltConflict, ...
    'frontVotes',   frontVotes, ...
    'backVotes',    backVotes, ...
    'unknownVotes', unknownVotes, ...
    'position',     position);

end

%% Helper : Dernier front descendant

function idx = lastDescEdge(seg, thrHigh, thrLow)
idx = [];
wasHigh = false;
for i = 1:numel(seg)
    if seg(i) > thrHigh
        wasHigh = true;
    elseif seg(i) < thrLow && wasHigh
        idx = i;
        wasHigh = false;
    end
end
end

%% Helper : Etat haut/bas ou indéterminé d'un signal à un échantillon donné

function st = stateAt(seg, idx, thrHigh, thrLow)
v = seg(idx);
if v < thrLow
    st = "low";
elseif v > thrHigh
    st = "high";
else
    st = "unknown";
end
end
