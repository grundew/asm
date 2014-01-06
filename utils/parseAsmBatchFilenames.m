function results = parseAsmBatchFilenames(filenames, filenamevars, prefix)

results = struct();
if exist('prefix', 'var') && ~isempty(prefix)
    idstart = length(prefix) + 2;
else
    idstart = 1;
end

for i = 1:length(filenamevars)
    var = filenamevars{i};
    [~, endidx] = arrayfun(@(x) regexp(x.name(idstart:end), filenamevars{i}), filenames);
    assert(length(unique(endidx))==1, 'Somethings fishy with the filenames');
    endidx = endidx + idstart - 1;
    data = arrayfun(@(x) textscan(x.name(endidx:end), '%*s %f %*[^\n]', 'delimiter', '_'), filenames, 'uni', 1);
    results.(var) = cell2mat(data);
end

end