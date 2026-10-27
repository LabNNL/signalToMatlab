function MEP = attachBrainsightErrors(MEP, matchedRows, BrainsightTable)

% Fills MEP.(lab).Brainsight for every MEP, using matchedRows (from 
% matchBrainsightErrors) and BrainsightTable (from readBrainsightSamples).
% NaN/false where the corresponding stim has no row.

names = fieldnames(MEP);
names = names(startsWith(names,'MEP_'));

for k = 1:numel(names)
    lab = names{k};
    t = MEP.(lab).orig_idx;
    row = matchedRows(t);
    if isnan(row)
        MEP.(lab).Brainsight = emptyBrainsight();
    else
        MEP.(lab).Brainsight.TargetError_mm   = BrainsightTable.TargetError(row);
        MEP.(lab).Brainsight.AngularError_deg = BrainsightTable.AngularError(row);
        MEP.(lab).Brainsight.TwistError_deg   = BrainsightTable.TwistError(row);
        MEP.(lab).Brainsight.Distance_mm      = BrainsightTable.Distance(row);
        MEP.(lab).Brainsight.SampleIndex      = BrainsightTable.Index(row);
        MEP.(lab).Brainsight.Time             = BrainsightTable.Time(row);
        MEP.(lab).Brainsight.Matched          = true;
    end
end

end

function bs = emptyBrainsight()
bs.TargetError_mm   = NaN;
bs.AngularError_deg = NaN;
bs.TwistError_deg   = NaN;
bs.Distance_mm      = NaN;
bs.SampleIndex      = NaN;
bs.Time             = NaN;
bs.Matched          = false;
end
