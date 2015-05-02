function exportfigure(outputdir, filename, fig, varargin)
% exportfigure(filename, fig)
% 
% Export figure handle fig as a tikz-figure to filename.

switch nargin
    case 2
        fig = gcf;
        floatFormat = '%.10g';
    otherwise
        floatFormat = '%.10g';
end

extraCode = {...
    '\usepackage{siunitx}',...
    '\newlength\figureheight',...
    '\newlength\figurewidth',...
    '\setlength\figureheight{4cm}', ...
    '\setlength\figurewidth{11cm}'};

% Crop data to whats visible
hline = findobj(fig, 'type', 'line');
hax = findobj(fig, 'type', 'axes');
xlim = get(hax, 'XLim');
title(hax, '')
xx = get(hline, 'XData');
yy = get(hline, 'YData');

if iscell(xx)
    for i = 1:length(xx)
        x = xx{i};
        y = yy{i};
        
        idlim = x>=xlim(1) & x<=xlim(2);
        set(hline(i), 'XData', x(idlim));
        set(hline(i), 'YData', y(idlim));
    end
else
    idlim = xx>=xlim(1) & xx<=xlim(2);
    set(hline, 'XData', xx(idlim));
    set(hline, 'YData', yy(idlim));
end
    

fn = fullfile(outputdir, filename);
matlab2tikz('filename', fn,...
    'figurehandle', fig,...
    'height', '\figureheight',...
    'width', '\figurewidth',...
    'extraCode', extraCode,...
    'standalone', true,...
    'showInfo', false,...
    'floatFormat', floatFormat,...
    'checkForUpdates', false,...
    'parseStrings', false,...
    varargin{:});
end
