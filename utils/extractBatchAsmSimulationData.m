function [X, ff, p, var] = extractBatchAsmSimulationData(datDir, sortvar)

fn = dir(fullfile(datDir, '*.mat'));
var = parseAsmBatchFilenames(fn);
nfn = length(fn);

if exist('sortvar', 'var')
    assert(any(strcmp(fieldnames(var), sortvar)),...
        'HW:INPUTERROR',...
        'Variable %s doesn''t exist in the file names',...
        sortvar)
    [~, idsort] = sort([var.(sortvar)]);
    var = var(idsort);
    fn = fn(idsort);
end

matfn = fullfile(datDir, fn(1).name);
mat = load(matfn);

X = zeros(length(mat.pt), nfn);
ff = cell(1, nfn);
X(:, 1) = mat.pt;
ff{1} = matfn;
if isfield(mat, 'params')
    p(1) = mat.params;
else
    p = 0;
end

for j = 2:nfn
    matfn = fullfile(datDir, fn(j).name);
    mat = load(matfn);
    X(:, j) = mat.pt;
    if isfield(mat, 'params')
        p(j) = mat.params;
    end
    ff{j} = matfn;
end

end