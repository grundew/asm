function results = parseAsmBatchFilenames(filenames)

nfn = length(filenames);

assert(~isempty(nfn)&&nfn>0, 'HW:INPUTERROR', 'No files found');

results(nfn) = struct();
for i = 1:length(filenames)
    [~, fn] = fileparts(filenames(i).name);
    n = length(strfind(fn, '_')) + 1;
    c = textscan(fn, '%s', n, 'Delimiter', '_');
    s = cellfun(@str2num, c{:}, 'uni', 0);
    isstrid = cellfun(@isempty, s, 'uni', 1);
    s(isstrid) = c{1}(isstrid);
    s(end-3:end) = c{1}(end-3:end);
    results(i).prefix = s{1};
    for j = 2:2:length(s)-4
        results(i).(s{j}) = s{j+1};
    end
    timestr = strjoin(s(end-3:end)', '-');
    datev = datevec(timestr, 'dd-mm-yyyy-HHMMSS');
    results(i).time = datestr(datev);
end

end