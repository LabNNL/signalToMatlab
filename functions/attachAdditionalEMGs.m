function MEP = attachAdditionalEMGs(MEP, data, selected_EMGs, EMG_field, MEPWindows, freq_EMG)

% Fills MEP.(lab).EMG_add.(EMG_xx) with the EMG channels that were recorded
% but not analysed. Each channel is filtered like the analysed one and cut
% in the same window.

addEMGs = setdiff(string(selected_EMGs), string(EMG_field), 'stable');

names = fieldnames(MEP);
names = names(startsWith(names,'MEP_'));

for e = 1:numel(addEMGs)
    ch = char(addEMGs(e));
    chData = data.(ch);

    if chData.FreqS ~= freq_EMG || numel(chData.dat) ~= numel(data.(EMG_field).dat)
        warning('attachAdditionalEMGs: %s skipped: sampling frequency or length differs from %s.', ch, EMG_field);
        continue
    end

    chFiltered = filtrage(chData.dat, freq_EMG, 20, 1000);

    for k = 1:numel(names)
        lab = names{k};
        win = MEPWindows(MEP.(lab).orig_idx, :);
        seg = chFiltered(win(1):win(2));

        sig = nan(numel(MEP.(lab).EMG), 1);
        sig(1:numel(seg)) = seg;
        MEP.(lab).EMG_add.(ch) = sig;
    end
end

end
