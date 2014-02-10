function out = extractBatchAsmSimulationData(datDir, sortfilenamevar, splitfilenamevar)

if nargin < 2
    error('HW:INPUTERROR', 'Not enough input paramaters');
end

fn = dir(fullfile(datDir, '*.mat'));
var = parseAsmBatchFilenames(fn);

if ~isfield(var, splitfilenamevar) || ~isfield(var, sortfilenamevar)
    error('HW:INPUTERROR', 'Filename variables does not exist')
end

splitvar = [var.(splitfilenamevar)];
unisplitvar = unique(splitvar);
nsplit = length(unisplitvar);
out(nsplit) = struct('X', [], sortfilenamevar, [], splitfilenamevar, [], 'f', [], 'p', [], 'fn', '');

fnsplit = cell(1, nsplit);
vars = cell(1, nsplit);
splitvars = cell(1, nsplit);
for i = 1:nsplit
    splitvarid = (splitvar==unisplitvar(i));
    fnsplit{i} = fn(splitvarid);
    vars{i} = var(splitvarid);
    splitvars{i} = splitvar(splitvarid);
end

for i = 1:nsplit

    fns = fnsplit{i};
    v = vars{i};
    [vsort, idsort] = sort([v.(sortfilenamevar)]);
    fnsort = fns(idsort);
    matfn = fullfile(datDir, fnsort(1).name);
    mat = load(matfn);

    xx = zeros(length(mat.pt), length(fnsort));
    ff = cell(1, length(fns));
    xx(:, 1) = mat.pt;
    ff{1} = matfn;
    p(1) = mat.params;
    for j = 2:length(fnsort)
        matfn = fullfile(datDir, fnsort(i).name);
        mat = load(matfn);
        xx(:, i) = mat.pt;
        p(j) = mat.params;
        ff{j} = matfn;
    end

    out(i).X = xx;
    out(i).fn = ff;
    out(i).(splitfilenamevar) = splitvars{i};
    out(i).(sortfilenamevar) = vsort;
    out(i).p = p;
    out(i).f = mat.f;
end

end