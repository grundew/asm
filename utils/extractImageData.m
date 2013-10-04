function [C, x, y] = extractImageData(himage)
% [C, x, y] = extractImageData(himage)
%
% Extracts the data in an image object.
% 
% Input:
% himage - Image handle
% 
% Output:
% C - Color data ('CData')
% x - x-axis ('XData')
% y - y-axis ('YData')

assert(~isempty(himage) & ishandle(himage) & strcmp(get(himage, 'Type'), 'image'),...
    'HW:InputError', 'Input or gco needs to be an image handle');

C = get(himage, 'CData');
x = get(himage, 'XData');
y = get(himage, 'YData');
end