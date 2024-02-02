% Helper functions

function ID = generate_mixed_traffic_flow(n, m)
    % Validate input
    if m > n || n < 0 || m < 0
        error('Invalid input: ensure that 0 <= m <= n');
    end

    % Step 1: Create the vector
    % Initialize a vector with m ones and n-m zeros
    ID = [ones(m, 1); zeros(n - m, 1)];

    % Step 2: Shuffle the vector
    % Generate a random permutation of the vector
    shuffledIndices = randperm(n);
    ID = ID(shuffledIndices);
end