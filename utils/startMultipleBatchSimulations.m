function startMultipleBatchSimulations(prefix, primvarnames, primvar,...
    secndvarnames, secndvar)

% Check and parse input parameters
assert(length(primvarnames)==length(primvar),...
    'HW:InputError', 'Length of primvarnames and primvar must be the same');
assert(length(secndvarnames)==length(secndvar),...
    'HW:InputError', 'Length of secndvarnames and secndvar must be the same');
assert(ischar(prefix), 'HW:InputError', 'prefix must be a string');

np = length(primvar{1});
ns = length(secndvar{1});
curdir = pwd();

for i = 1:ns
    sv = cell(1, 2*length(secndvar));
    sv(1:2:end) = secndvarnames;
    for kk = 1:length(secndvar)
        sv{2*kk} = secndvar{kk}(i);
    end
    
    outdir = genoutdir(prefix, sv{:});
    mkdir(outdir);
    cd(outdir);
    for j = 1:np
        pv = cell(1, 2*length(primvar));
        pv(1:2:end) = primvarnames;
        for kk = 1:length(secndvar)
            pv{2*kk} = primvar{kk}(j);
        end
        vv = cat(2, pv, sv);
        startAsmSimulation(vv{:});
    end

end

cd(curdir);
end

function d = genoutdir(prefix, varargin)
params = varargin(1:2:end);
values = varargin(2:2:end);

c = cell(1, length(params));
for i = 1:length(params)
    c{i} = sprintf('%s_%d', params{i}, values{i});
end
paramstr = strjoin(c, '_');
d = sprintf('%s_%s', prefix, paramstr);
end