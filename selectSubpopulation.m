function [subpopulation, L ] = selectSubpopulation(pop, fitness)
    % Parameters
    L_min = 15; % Minimum subpopulation size
    L_max = size(pop, 1) / 2; % Maximum subpopulation size, half of the total population
    D_max = 50; % Maximum expected dimensionality for normalization
    V_max = var(fitness); % Normalization factor for variability
    alpha = 0.5; % Weight for dimensionality influence
    beta = 0.5; % Weight for variability influence
   
    % Problem dimensionality
    D = size(pop, 2);
    
    % Population performance variability
    V = std(fitness);

    % Calculate the dynamic subpopulation size
    fDV = alpha * (D / D_max) + beta * (V / V_max);
    L = max(L_min, min(L_max, ceil(size(pop, 1) * fDV)));
    
    % Select subpopulation
    subpopulation = pop(1:L, :);
end
