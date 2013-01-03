% Buckingham Pi theorem applied to the transmission coefficient
B = [-3, -3, 1, 1, 0, 0;1, 1, 0, 0, 0, 0; 0, 0, -1, -1, 0, 0];

n = size(B, 2);
r = rank(B);

% Number of diemensionless quantities
k = n - r;


% Pressure
Bp = [1, 1, 0, 0, 0, 0, 0, 0, 0; -3, -3, 1, 0, 1, 1, 1, 1, 0; 0, 0, -1, 0, -1, 0, 0, 0, -1];
nb = size(Bp, 2);
rb = rank(Bp);

% Number of diemensionless quantities
kb = nb - rb