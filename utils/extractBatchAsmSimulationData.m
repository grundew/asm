function [X, varsort, f, fres] = extractBatchAsmSimulationData(datDir, filenamevar)


fn = dir(fullfile(datDir, '*.mat'));
var = parseAsmBatchFilenames(fn, {filenamevar});
[varsort, idsort] = sort(var.(filenamevar));
fnsort = fn(idsort);

for i = 1:length(fnsort)
    matfn = fullfile(datDir, fnsort(i).name);
    mat = load(matfn);
    X(:, i) = mat.pt;
end
f = (0:mat.params.nfft-1)*mat.params.fs/mat.params.nfft;
fres = mat.params.cp/2/mat.params.thickness;
end

