function extractS1FromMultipleBatch(inputdir, dirnamevar, filenamevar, prefix)

if nargin < 4
    prefix = '';
end
d = dir(inputdir);
d = d([d.isdir]);
idndot = arrayfun(@(x) any(strcmp(x.name, {'.', '..'})), d);
d = d(~idndot);

var = parseAsmBatchFilenames(d, {dirnamevar}, prefix);
[dirvarsort, idsort] = sort(var.(dirnamevar));
dsort = d(idsort);

fig_s1 = figure;
ax_s1 = axes('Parent', fig_s1);
hold(ax_s1, 'all');
legend(ax_s1, 'show');

fig_s1_s = figure;
ax_s1_s = axes('Parent', fig_s1_s);
hold(ax_s1_s, 'all');
legend(ax_s1_s, 'show');
for i = 1:length(dsort)
    [X, fnvarsort, f, fres, p] = extractBatchAsmSimulationData(dsort(i).name, filenamevar);
    fig = figure;
    ax = axes('Parent', fig);
    Xn = abs(X)./repmat(max(abs(X), [], 1), size(X, 1), 1);
    himage = imagesc(f/fres, fnvarsort, Xn', 'Parent', ax);
    xlabel(ax, 'f/f_1');
    ylabel(ax, filenamevar);
    [x, y] = findResonanceInImage(himage, [0.7, 1.2], ylim(ax));
    plot(ax_s1, y, x,...
        'DisplayName', sprintf('%s - %d', dirnamevar, dirvarsort(i)));
    % Fresnel distance
    d_f = y.^2/p(1).cf*fres;
    plot(ax_s1_s, y./d_f, x, ...
        'DisplayName', sprintf('%s - %d', dirnamevar, dirvarsort(i)));
end
figure(fig_s1);
figure(fig_s1_s);
end