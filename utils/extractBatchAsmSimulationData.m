function [X, varsort, f, fres, p] = extractBatchAsmSimulationData(datDir, filenamevar)
fn = dir(fullfile(datDir, '*.mat'));
var = parseAsmBatchFilenames(fn, {filenamevar});
[varsort, idsort] = sort(var.(filenamevar));
fnsort = fn(idsort);

for i = 1:length(fnsort)
    matfn = fullfile(datDir, fnsort(i).name);
    mat = load(matfn);
    X(:, i) = mat.pt;
    p(i) = mat.params;
end

if isfield(mat.params, 'f')
    f = mat.params.f;
else
    f = (0:mat.params.nfft-1)*mat.params.fs/mat.params.nfft;
end

fres = mat.params.cp/2/mat.params.thickness;
end

