function [x, y] = findResonanceInImage(himage, dim, xl, yl)

[C, xx, yy] = extractImageData(himage);

if ~exist('dir', 'var')
    dim = 1;
end

% Crop Cdata to limits yl and xl
if ~exist('xl', 'var')
    ax = get(himage, 'Parent');
    xl = xlim(ax);
end

if ~exist('yl', 'var')
    ax = get(himage, 'Parent');
    yl = ylim(ax);
end

idx = xx>=xl(1) & xx<=xl(2);
idy = yy>=yl(1) & yy<=yl(2);
C = C(idy, idx);
x = xx(idx);
y = yy(idy);

[~, idmax] = max(C, [], dim);

if dim==1
    y = y(idmax);
elseif dim==2
    x = x(idmax);
else
    error('HW:InputError', 'dim must be either 1 or 2');
end

end