% Define Rastrigin function 
function val = rastrigin(x)
    A = 10;
    n = numel(x);
    val = A * n + sum(x.^2 - A * cos(2 * pi * x));
end
