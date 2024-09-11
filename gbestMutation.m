function mutatedSolutions = gbestMutation(gbest, dimensions, numRequired)
    % Placeholder for gbest mutation logic
    % This should generate 'numRequired' new solutions based on 'gbest'
    % For demonstration, using simple random perturbation
    sigma = 0.1; % Mutation strength, adjust as necessary
    mutatedSolutions = repmat(gbest, numRequired, 1) + sigma * randn(numRequired, dimensions);
end
