a = 0;
b = pi;
N = 2.^(2:23);
facit = 2;

z = zeros(length(N), 1);
zm = zeros(length(N), 1);
zsimp = zeros(length(N), 1);
tz = zeros(length(N), 1);
tzm = zeros(length(N), 1);
tsimp = zeros(length(N), 1);
zsimp38 = zeros(length(N), 1);
tsimp38 = zeros(length(N), 1);
zsimp38i = zeros(length(N), 1);
tsimp38i = zeros(length(N), 1);

for i = 1:length(N);
    M = 3*N(i);
    H = (b - a)/M;
    x = linspace(a, b, M+1);
    y = sin(x);
    
    % Trapezoidal (built-in)
    tic;
    z(i) = trapz(x, y);
    tz(i) = toc;
    
    % Trapezoidal
    tic;
    I = (0.5*(y(1) + y(end)) + sum(y(2:end-1)));
    for j = 1:10
        xa = j*H/10 + H*(0:M-1);
        zm(i) = zm(i) + sum(sin(xa));
    end
    zm(i) = H*I;
    tzm(i) = toc;
    
    % Simpson formula
    tic
    zsimp(i) = 2*H*(1/6*(y(1) + y(end)) + ...
        1/3*sum(y(3:2:end-2)) + ...
        2/3*sum(y(2:2:end-1)));
    tsimp(i) = toc;
    
    % Simpson's 3/8
    tic
    zsimp38(i) = 3/8*H*(sum(y(1:3:end-2)) +...
        3*(sum(y(2:3:end-2)) + sum(y(3:3:end-1))) + sum(y(4:3:end)));
    tsimp38(i) = toc;
    
    % Improved Simpson's 3/8
    tic
    zsimp38i(i) = 3/8*H*(sum(y(1:3:end-2)) +...
        3*(sum(y(2:3:end-2)) + sum(y(3:3:end-1))) + sum(y(4:3:end))) + ...
        3/160*H*(-sum(y(3:3:end-5)) + 23*sum(y(4:3:end-6)) +...
        58*(sum(y(5:3:end-4)) + sum(y(6:3:end-4))) + 23*sum(y(7:3:end-3)) +...
        - sum(y(8:3:end-2)));
    tsimp38i(i) = toc;
    
end

loglog(N, abs(z-facit), N,...
    abs(zm-facit), '.',...
    N, abs(zsimp-facit), 'x',...
    N, abs(zsimp38-facit),...
    N, abs(zsimp38i-facit))
figure
loglog(N, tz, N, tzm, N, tsimp, N, tsimp38, N, tsimp38i)
ylabel('Time')
xlabel('N')