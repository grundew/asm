function results = parseAsmBatchFilenames(filenames, filenamevars)

results = struct();

for i = 1:length(filenamevars)
    var = filenamevars{i};
    [~, endidx] = arrayfun(@(x) regexp(x.name, filenamevars{i}), filenames);
    assert(length(unique(endidx))==1, 'Somethings fishy with the filenames');
    data = arrayfun(@(x) textscan(x.name(endidx+2:end), '%f'), filenames, 'uni', 1);
    results.(var) = cell2mat(data);
end

end