function startMultipleBatchSimulations(prefix, primvarnames, primvar,...
    secndvarnames, secndvar, p)

% Check and parse input parameters
assert(length(primvarnames)==length(primvar),...
    'HW:InputError', 'Length of primvarnames and primvar must be the same');
assert(length(secndvarnames)==length(secndvar),...
    'HW:InputError', 'Length of secndvarnames and secndvar must be the same');
assert(ischar(prefix), 'HW:InputError', 'prefix must be a string');

np = length(primvar{1});
ns = length(secndvar{1});
curdir = pwd();
p.savemat = true;
for i = 1:ns
    
    % Update the secondary variable
    for kk = 1:length(secndvarnames)
        p.(secndvarnames{kk}) = secndvar{kk}(i);
    end
    
    % Make a new directory for the results
    outdir = genoutdir(prefix, secndvarnames, p);
    mkdir(outdir);
    cd(outdir);
    for j = 1:np
        
        % Update primary variable
        for kk = 1:length(primvarnames)
            p.(primvarnames{kk}) = primvar{kk}(j);
        end
        
        startAsmSimulation(p);
        
    end
    
    cd('../');
end

cd(curdir);
end

function d = genoutdir(prefix, params, p)
c = cell(1, length(params));
for i = 1:length(params)
    c{i} = sprintf('%s_%d', params{i}, p.(params{i}));
end
paramstr = strjoin(c, '_');
d = sprintf('%s_%s', prefix, paramstr);
end