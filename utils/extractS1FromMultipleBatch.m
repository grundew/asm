function extractS1FromMultipleBatch(inputdir, dirnamevar, filenamevar)

d = dir(inputdir);
d = d([d.isdir]);
idndot = arrayfun(@(x) any(strcmp(x.name, {'.', '..'})), d);
d = d(~idndot);

var = parseAsmBatchFilenames(d, {dirnamevar});
[dirvarsort, idsort] = sort(var.(dirnamevar));
dsort = d(idsort);

fig_s1 = figure;
ax_s1 = axes('Parent', fig_s1);
hold(ax_s1, 'all');
legend(ax_s1, 'show');
for i = 1:length(dsort)
    [X, fnvarsort, f, fres] = extractBatchAsmSimulationData(dsort(i).name, filenamevar);
    fig = figure;
    ax = axes('Parent', fig);
    Xn = abs(X)./repmat(max(abs(X), [], 1), size(X, 1), 1);
    himage = imagesc(f/fres, fnvarsort, Xn', 'Parent', ax);
    ylabel(ax, 'f/f_1');
    xlabel(ax, filenamevar);
    [x, y] = findResonanceInImage(himage, [0.7, 1.2], ylim(ax));
    plot(ax_s1, y, x,...
        'DisplayName', sprintf('%s - %d', dirnamevar, dirvarsort(i)));
end
figure(fig_s1);
end