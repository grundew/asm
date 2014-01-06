function [X, varsort, f, fres, p] = extractBatchAsmSimulationData(datDir, filenamevar, prefix)

if nargin < 3
    prefix = '';
end
fn = dir(fullfile(datDir, '*.mat'));
var = parseAsmBatchFilenames(fn, {filenamevar}, prefix);
[varsort, idsort] = sort(var.(filenamevar));
fnsort = fn(idsort);

matfn = fullfile(datDir, fnsort(1).name);
mat = load(matfn);
X = zeros(length(mat.pt), length(fnsort));
X(:, 1) = mat.pt;
p([1, length(fnsort)]) = mat.params;
for i = 2:length(fnsort)
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

