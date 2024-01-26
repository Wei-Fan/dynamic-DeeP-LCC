function HDVs = generate_HDVs(n)
    % Initialize the struct
    HDVs = struct();
    HDVs.type   = 1;
    HDVs.alpha  = 0.4 + (0.8 - 0.4) * rand(n,1);
    HDVs.beta   = 0.7 + (1.1 - 0.7) * rand(n,1);
    HDVs.s_st   = 5;
    HDVs.s_go   = 30 + (40 - 30) * rand(n,1);
    HDVs.v_max  = 30;
    HDVs.s_star = 17.5 + (22.5 - 17.5) * rand(n,1);
end
