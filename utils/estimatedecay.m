function tau = estimatedecay(x, t, idxStart, idxEnd)
% Estimate decay
%   Detailed explanation goes here

ta = t(idxStart:idxEnd);
xa = abs(x(idxStart:idxEnd, :));

tau = zeros(1, size(xa, 2));
for i = 1:size(xa, 2)
    p = polyfit(ta, log(xa(:, i)), 1);
    tau(i) = p(1);
end

end