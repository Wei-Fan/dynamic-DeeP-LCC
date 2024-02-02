% =========================================================================
%               Dynamic Traffic Flow Simulation
%               100 Vehicles with 20 CAVs
%               Vehicles are added or removed from the main traffic flow
%
% *Attention*: need pre-collected trajectory data to run this simulation. 
% If there are no pre-collected data, please run data_collection.m first. 
% =========================================================================

clc; close all; clear;
addpath('_fcn');
addpath('_data');
warning off;

trial_name = 'one_flow';

% whether output data
output_bool = 1;

% length of data samples
T           = 600;

% Whether update equilibrium velocity
% If not, the CAV has a prescribed equilibrium velocity
update_equilibrium_bool = 0;


% -------------------------------------------------------------------------
%   Parameter setup
% -------------------------------------------------------------------------

h_wait = waitbar(0,'please wait');

% Number of data sets for simulation
data_number         = 1;    % up to 100
% Perturbation amplitude
per_type            = 3;    % 1. sinuoid perturbation 2. small brake perturbation 3. large brake perturbation
                            % 4. larger brake perturbation
                            % 5. Perturbation on a vehicle in the middle of the platoon
sine_amp            = 4; % amplitidue of sinuoid perturbation
brake_amp           = 5; % brake amplitude of brake perturbation

constraint_bool     = 1; % Whether there exist constraints

% Type for HDV car-following model
hdv_type            = 1;    % 1. OVM   2. IDM
% Uncertainty for HDV behavior
acel_noise          = 0.1;  % A white noise signal on HDV's original acceleration

% Parameters in Simulation
total_time          = 200;              % Total Simulation Time
Tstep               = 0.05;             % Time Step
total_time_step     = total_time / Tstep;

% Index for one experiment
computation_time    = zeros(total_time_step, 1);
iteration_num       = zeros(total_time_step, 1);

% Average index for all the experiments
Collected_computation_time  = zeros(data_number, 1);
Collected_iteration_num     = zeros(data_number, 1);

% DeeP-LCC Formulation
Tini                = 20;        % length of past data
N                   = 50;        % length of predicted horizon

% Weight coefficients
weight_choice       = 3;
% case for weight choice in centralized DeeP-LCC

switch weight_choice
    case 1
        weight_v     = 1;        % weight coefficient for velocity error
        weight_s     = 0.5;      % weight coefficient for spacing error
        weight_u     = 0.1;      % weight coefficient for control input
        lambda_g     = 10;       % penalty on ||g||_2^2 in objective
        lambda_y     = 1e4;      % penalty on ||sigma_y||_2^2 in objective
    case 2
        weight_v     = 4;        % weight coefficient for velocity error
        weight_s     = 2;        % weight coefficient for spacing error
        weight_u     = 0.4;      % weight coefficient for control input
        lambda_g     = 10;       % penalty on ||g||_2^2 in objective
        lambda_y     = 1e4;      % penalty on ||sigma_y||_2^2 in objective
 case 3
        weight_v     = 1;        % weight coefficient for velocity error
        weight_s     = 0.5;      % weight coefficient for spacing error
        weight_u     = 0.1;      % weight coefficient for control input
        lambda_g     = 50;       % penalty on ||g||_2^2 in objective
        lambda_y     = 1e4;      % penalty on ||sigma_y||_2^2 in objective
end

% ------------------------------------------
% Parameters in Mixed Traffic
% ------------------------------------------
load(['_data/ID_', trial_name, '_n_100_tmp.mat']);
n_total     = length(ID);            % number of vehicles
m_total     = length(find(ID==1));
Omega_c = find(ID==1);
m_separate  = randi([2,m_total]);
n_separate  = Omega_c(m_separate);

ID1         = ID(1:n_separate - 1);
ID2         = ID(n_separate : end);

Omega1_c    = find(ID1==1);          % position of CAVs
n1          = length(ID1);           % number of vehicles
m1          = length(Omega1_c);      % number of CAVs

Omega2_c    = find(ID2==1);          % position of CAVs
n2          = length(ID2);           % number of vehicles
m2          = length(Omega2_c);      % number of CAVs


% Constraints
acel_max        = 2;
dcel_max        = -5;
spacing_max     = 40;
spacing_min     = 5;
u_limit         = [dcel_max, acel_max];
s_limit         = [spacing_min, spacing_max] - s_star;

% Random setup for OVM
load(['_data/hdv_', trial_name, '_n_100_tmp.mat']);
% Initialize two empty structs
hdv_parameter1 = struct();
hdv_parameter2 = struct();

% Loop through the fields of the input struct
fields = fieldnames(hdv_parameter);
for i = 1:numel(fields)
    field = fields{i};
    data = hdv_parameter.(field);

    % Check if the field is a vector
    if isvector(data)
        n = length(data);
        if n > n2
            % Separate the vector into two parts, n-n2 and n2
            hdv_parameter1.(field) = data(1:n1);
            hdv_parameter2.(field) = data((n1 + 1):end);
        else
            % If the vector has 10 or fewer elements, include it in both structs
            hdv_parameter1.(field) = data;
            hdv_parameter2.(field) = data;
        end
    else
        % If the field is not a vector, include it in both structs
        hdv_parameter1.(field) = data;
        hdv_parameter2.(field) = data;
    end
end

% What is measurable
measure_type = 3;
% 1. Only the velocity errors of all the vehicles are measurable;
% 2. All the states, including velocity error and spacing error are measurable;
% 3. Velocity error and spacing error of the CAVs are measurable,
%    and the velocity error of the HDVs are measurable.

% ------------------
%  size in DeeP-LCC
% ------------------

n1_ctr = 2 * n1;      % number of state variables
m1_ctr = m1;          % number of input variables
p1_ctr = n1 + m1;     % number of output variables

n2_ctr = 2 * n2;      % number of state variables
m2_ctr = m2;          % number of input variables
p2_ctr = n2 + m2;     % number of output variables

for i_data = 1:data_number
    % Load trajectory data
    load(['./_data/trajectory_data_collection/', trial_name, '_data', num2str(i_data), '_T_',num2str(T),'_',num2str(Tini),'_',num2str(N),'_noiseLevel_',num2str(acel_noise),'_tmp.mat']);
    % ---------------------------------------
    %   Scenario initialization
    %----------------------------------------

    % There is one head vehicle at the very beginning
    S1          = zeros(total_time_step, n1 + 1, 3);
    S1(1,1,1)    = 0;
    for i = 2 : n1 + 1
        S1(1,i,1) = S1(1,i-1,1) - hdv_parameter1.s_star(i-1);
    end
    S1(1,:,2)    = v_star * ones(n1 + 1,1);

    S2          = zeros(total_time_step, n2 + 1, 3);
    S2(1,1,1)    = 0;
    for i = 2 : n2 + 1
        S2(1,i,1) = S1(1,i-1,1) - hdv_parameter2.s_star(i-1);
    end
    S2(1,:,2)    = v_star * ones(n2 + 1,1);

    % ------------------
    %  DeeP-LCC Formulation
    % ------------------
    Q1_v         = weight_v*eye(n1);              % penalty for velocity error
    Q1_s         = weight_s*eye(p1_ctr-n1);       % penalty for spacing error
    Q1           = blkdiag(Q1_v,Q1_s);            % penalty for trajectory error
    R1           = weight_u*eye(m1_ctr);          % penalty for control input

    u1           = zeros(m1_ctr,total_time_step); % control input
    x1           = zeros(n1_ctr,total_time_step); % state variables
    y1           = zeros(p1_ctr,total_time_step); % output variables
    pr_status1   = zeros(total_time_step,1);      % problem status
    e1           = zeros(1,total_time_step);      % external input

    Q2_v         = weight_v*eye(n2);              % penalty for velocity error
    Q2_s         = weight_s*eye(p2_ctr-n2);       % penalty for spacing error
    Q2           = blkdiag(Q2_v,Q2_s);            % penalty for trajectory error
    R2           = weight_u*eye(m2_ctr);          % penalty for control input

    u2           = zeros(m2_ctr,total_time_step); % control input
    x2           = zeros(n2_ctr,total_time_step); % state variables
    y2           = zeros(p2_ctr,total_time_step); % output variables
    pr_status2   = zeros(total_time_step,1);      % problem status
    e2           = zeros(1,total_time_step);      % external input

    % ------------------
    %  Separate Hankel Matrix
    % ------------------
    ud1 = ud(1:m1, :);
    ud2 = ud(m1+1:end, :);
    ed1 = ed;
    ed2 = yd(n1, :);
    yd1 = [yd(1:n1, :); yd(n1+n2+1:n1+n2+m1, :)];
    yd2 = [yd(n1+1:n1+n2, :); yd(n1+n2+m1+1:end, :)];

    U1   = hankel_matrix(ud1, Tini + N);
    U1p  = U1(1:Tini * m1_ctr,:);
    U1f  = U1((Tini * m1_ctr + 1):end,:);

    E1   = hankel_matrix(ed1, Tini + N);
    E1p  = E1(1:Tini,:);
    E1f  = E1((Tini + 1):end,:);

    Y1   = hankel_matrix(yd1, Tini + N);
    Y1p  = Y1(1:Tini * p1_ctr,:);
    Y1f  = Y1((Tini * p1_ctr + 1):end,:);

    U2   = hankel_matrix(ud2, Tini + N);
    U2p  = U2(1:Tini * m2_ctr,:);
    U2f  = U2((Tini * m2_ctr + 1):end,:);

    E2   = hankel_matrix(ed2, Tini + N);
    E2p  = E2(1:Tini,:);
    E2f  = E2((Tini + 1):end,:);

    Y2   = hankel_matrix(yd2, Tini + N);
    Y2p  = Y2(1:Tini * p2_ctr,:);
    Y2f  = Y2((Tini * p2_ctr + 1):end,:);

    % ------------------
    %  Reference trajectory
    % ------------------
    r1       = zeros(p1_ctr, total_time_step + N);            % stabilization
    r2       = zeros(p2_ctr, total_time_step + N);            % stabilization

    % ---------------------------------------------------------------------
    %   Simulation starts here
    %----------------------------------------------------------------------
    tic

    % ------------------
    %  Initial trajectory
    % ------------------
    u1ini = zeros(m1_ctr,Tini);
    e1ini = zeros(1,Tini);
    y1ini = zeros(p1_ctr,Tini);

    u2ini = zeros(m2_ctr,Tini);
    e2ini = zeros(1,Tini);
    y2ini = zeros(p2_ctr,Tini);

    for k = 1:Tini-1
        % Update acceleration for flow 1
        acel1 = HDV_dynamics(S1(k,:,:),hdv_parameter1) ...
              - acel_noise + 2*acel_noise*rand(n1,1);

        S1(k,1,3)           = 0;               % the head vehicle
        S1(k,2:end,3)       = acel1;           % all the vehicles using HDV model
        S1(k,Omega1_c+1,3)  = u1ini(:,k);      % the CAV

        S1(k+1,:,2) = S1(k,:,2) + Tstep*S1(k,:,3);
        S1(k+1,1,2) = e1ini(k) + v_star;       % the velocity of the head vehicle
        S1(k+1,:,1) = S1(k,:,1) + Tstep*S1(k,:,2);

        y1ini(:,k) = measure_mixed_traffic(S1(k,2:end,2),S1(k,:,1),ID1,v_star,s_star,measure_type);

        % Update acceleration for flow 2
        acel2 = HDV_dynamics(S2(k,:,:),hdv_parameter2) ...
              - acel_noise + 2*acel_noise*rand(n2,1);

        S2(k,1,3)           = 0;               % the head vehicle
        S2(k,2:end,3)       = acel2;           % all the vehicles using HDV model
        S2(k,Omega2_c+1,3)  = u2ini(:,k);      % the CAV

        S2(k+1,:,2) = S2(k,:,2) + Tstep*S2(k,:,3);
        S2(k+1,1,2) = e2ini(k) + v_star;       % the velocity of the head vehicle
        S2(k+1,:,1) = S2(k,:,1) + Tstep*S2(k,:,2);

        y2ini(:,k) = measure_mixed_traffic(S2(k,2:end,2),S2(k,:,1),ID2,v_star,s_star,measure_type);
    end

    k_end = k+1;
    y1ini(:,k_end) = measure_mixed_traffic(S1(k_end,2:end,2),S1(k_end,:,1),ID1,v_star,s_star,measure_type);
    y2ini(:,k_end) = measure_mixed_traffic(S2(k_end,2:end,2),S2(k_end,:,1),ID2,v_star,s_star,measure_type);

    u1(:,1:Tini) = u1ini;
    e1(:,1:Tini) = e1ini;
    y1(:,1:Tini) = y1ini;

    u2(:,1:Tini) = u2ini;
    e2(:,1:Tini) = e2ini;
    y2(:,1:Tini) = y2ini;

    % ------------------
    %  simulation starts here
    % ------------------
    for k = Tini:total_time_step-1
        % Update acceleration
        tic
        acel1 = HDV_dynamics(S1(k,:,:),hdv_parameter1) ...
              - acel_noise + 2*acel_noise*rand(n1,1);
        acel2 = HDV_dynamics(S2(k,:,:),hdv_parameter2) ...
              - acel_noise + 2*acel_noise*rand(n2,1);

        S1(k,2:end,3) = acel1;     % all the vehicles using HDV model
        S2(k,2:end,3) = acel2;     % all the vehicles using HDV model

        [u1_opt,y1_opt,pr1] = DeeP_LCC(U1p,Y1p,U1f,Y1f,E1p,E1f,u1ini,y1ini,e1ini,Q1,R1,r1(:,k:k+N-1),...
                           lambda_g,lambda_y,u_limit,s_limit);
        [u2_opt,y2_opt,pr2] = DeeP_LCC(U2p,Y2p,U2f,Y2f,E2p,E2f,u2ini,y2ini,e2ini,Q2,R2,r2(:,k:k+N-1),...
                           lambda_g,lambda_y,u_limit,s_limit);
        toc
        computation_time(k) = toc;

        % One-step formulation
        u1(:,k) = u1_opt(1:m1_ctr,1);
        u2(:,k) = u2_opt(1:m1_ctr,1);
        % Update accleration for the CAV
        S1(k,Omega1_c+1,3)   = u1(:,k);
        S2(k,Omega2_c+1,3)   = u2(:,k);
        % Judge whether SD system commands to brake
        brake_vehicle_ID1 = find(acel1==dcel_max);                 % the vehicles that need to brake
        brake_cav_ID1     = intersect(brake_vehicle_ID1,Omega1_c); % the CAVs that need to brake
        if ~isempty(brake_cav_ID1)
            S1(k,brake_cav_ID1+1,3) = dcel_max;
        end

        brake_vehicle_ID2 = find(acel2==dcel_max);                 % the vehicles that need to brake
        brake_cav_ID2     = intersect(brake_vehicle_ID2,Omega2_c); % the CAVs that need to brake
        if ~isempty(brake_cav_ID2)
            S2(k,brake_cav_ID2+1,3) = dcel_max;
        end

        % Update state
        S1(k+1,:,2) = S1(k,:,2) + Tstep*S1(k,:,3);
        S2(k+1,:,2) = S2(k,:,2) + Tstep*S2(k,:,3);

        % -------------
        % Perturbation for the head vehicle
        % -------------
        switch per_type
            case 1
                S1(k+1,1,2) = v_star + sine_amp*sin(2*pi/(10/Tstep)*(k-Tini));
                S1(k+1,:,1) = S1(k,:,1) + Tstep*S1(k,:,2);

                S2(k+1,1,2) = v_star + sine_amp*sin(2*pi/(10/Tstep)*(k-Tini));
                S2(k+1,:,1) = S2(k,:,1) + Tstep*S2(k,:,2);
            case 2
                if (k-Tini)*Tstep <= brake_amp/5
                    S1(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= brake_amp/5+5
                    S1(k+1,1,3) = 1;
                else
                    S1(k+1,1,3) = 0;
                end
                S1(k+1,:,2) = S1(k,:,2) + Tstep*S1(k,:,3);
                S1(k+1,:,1) = S1(k,:,1) + Tstep*S1(k,:,2);

                if (k-Tini)*Tstep <= brake_amp/5
                    S2(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= brake_amp/5+5
                    S2(k+1,1,3) = 1;
                else
                    S2(k+1,1,3) = 0;
                end
                S2(k+1,:,2) = S2(k,:,2) + Tstep*S2(k,:,3);
                S2(k+1,:,1) = S2(k,:,1) + Tstep*S2(k,:,2);
            case 3
                if (k-Tini)*Tstep <= brake_amp/5
                    S1(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= brake_amp/5+3
                    S1(k+1,1,3) = 0;
                elseif (k-Tini)*Tstep <= brake_amp/5+3+5
                    S1(k+1,1,3) = 1;
                else
                    S1(k+1,1,3) = 0;
                end
                S1(k+1,:,2) = S1(k,:,2) + Tstep*S1(k,:,3);
                S1(k+1,:,1) = S1(k,:,1) + Tstep*S1(k,:,2);

                if (k-Tini)*Tstep <= brake_amp/5
                    S2(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= brake_amp/5+3
                    S2(k+1,1,3) = 0;
                elseif (k-Tini)*Tstep <= brake_amp/5+3+5
                    S2(k+1,1,3) = 1;
                else
                    S2(k+1,1,3) = 0;
                end
                S2(k+1,:,2) = S2(k,:,2) + Tstep*S2(k,:,3);
                S2(k+1,:,1) = S2(k,:,1) + Tstep*S2(k,:,2);
            case 4
                if (k-Tini)*Tstep <= 2
                    S1(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= 2+5
                    S1(k+1,1,3) = 2;
                end
                S1(k+1,:,2) = S1(k,:,2) + Tstep*S1(k,:,3);
                S1(k+1,:,1) = S1(k,:,1) + Tstep*S1(k,:,2);

                if (k-Tini)*Tstep <= 2
                    S2(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= 2+5
                    S2(k+1,1,3) = 2;
                end
                S2(k+1,:,2) = S2(k,:,2) + Tstep*S2(k,:,3);
                S2(k+1,:,1) = S2(k,:,1) + Tstep*S2(k,:,2);
            case 5
                if (k-Tini)*Tstep <= brake_amp/2
                    S1(k,4,3) = -2;
                end
                S1(k+1,:,2) = S1(k,:,2) + Tstep*S1(k,:,3);
                S1(k+1,:,1) = S1(k,:,1) + Tstep*S1(k,:,2);

                if (k-Tini)*Tstep <= brake_amp/2
                    S2(k,4,3) = -2;
                end
                S2(k+1,:,2) = S2(k,:,2) + Tstep*S2(k,:,3);
                S2(k+1,:,1) = S2(k,:,1) + Tstep*S2(k,:,2);
        end

        % Record output
        y1(:,k) = measure_mixed_traffic(S1(k,2:end,2),S1(k,:,1),ID1,v_star,s_star,measure_type);
        e1(k)   = S1(k,1,2) - v_star;
        y2(:,k) = measure_mixed_traffic(S2(k,2:end,2),S2(k,:,1),ID2,v_star,s_star,measure_type);
        e2(k)   = S2(k,1,2) - v_star;

        % update past data in control process
        u1ini = u1(:,k-Tini+1:k);
        y1ini = y1(:,k-Tini+1:k);
        e1ini = S1(k-Tini+1:k,1,2) - v_star;
        u2ini = u2(:,k-Tini+1:k);
        y2ini = y2(:,k-Tini+1:k);
        e2ini = S2(k-Tini+1:k,1,2) - v_star;

        fprintf('Simulation number: %d  |  process... %2.2f%% \n',i_data,k/total_time_step*100);

        str=['Num:', num2str(i_data),' | Processing: ',num2str(k/total_time_step*100),'%'];
        waitbar(k/total_time_step,h_wait,str);
    end

    k_end = k+1;
    y1(:,k_end) = measure_mixed_traffic(S1(k_end,2:end,2),S1(k_end,:,1),ID1,v_star,s_star,measure_type);
    y2(:,k_end) = measure_mixed_traffic(S2(k_end,2:end,2),S2(k_end,:,1),ID2,v_star,s_star,measure_type);

    tsim = toc;

    fprintf('Simulation ends at %6.4f seconds \n', tsim);


    % -------------------------------------------------------------------------
    %   Plot Results
    %--------------------------------------------------------------------------

end
close(h_wait);
