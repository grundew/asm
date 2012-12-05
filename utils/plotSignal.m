function [axEcho, axTail] = plotSignal(t, x, ttail, xtail, f, X, Xtail, resonances)
%% Plot whole signal and spectrum
figure
axEcho(1) = subplot(211);
plot(t, real(x))
title('Time signal')
xlabel('Time (s)')
ylabel('Amplitude')
axEcho(2) = subplot(212);
plot(f, db(X), '.')
title('Spectrum')
xlabel('Frequency (Hz)')
ylabel('PSD')

%% Plot the tail and tailspectrum with thickness indications
figure
axTail(1) = subplot(211);
plot(ttail, real(xtail))
title('Tail signal')
xlabel('Time (s)')
ylabel('Amplitude')
axTail(2) = subplot(212);
plot(f, db(Xtail), '.')
yl = get(axTail(2), 'YLim');
resx = repmat(resonances(:)', 2, 1);
resy = repmat(yl(:), 1, length(resonances));
hold(axTail(2), 'all')
plot(axTail(2), resx, resy, 'r--')
title('Tail Spectrum')
xlabel('Frequency (Hz)')
ylabel('PSD')
end