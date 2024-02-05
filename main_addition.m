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

trial_name = 'two_flow';

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
load(['_data/ID_', trial_name, '_n1_100_n2_100_tmp.mat']);
n1       = length(ID1);            % number of vehicles
m1       = length(find(ID1==1));
Omega1_c = find(ID1==1);
n2       = length(ID2);            % number of vehicles
m2       = length(find(ID2==1));
Omega2_c = find(ID2==1);

ID       = [ID1; ID2];
n        = n1 + n2;
m        = m1 + m2;
Omega_c  = find(ID==1);

% Constraints
acel_max        = 2;
dcel_max        = -5;
spacing_max     = 40;
spacing_min     = 5;
u_limit         = [dcel_max, acel_max];
s_limit         = [spacing_min, spacing_max] - s_star;

% Random setup for OVM
load(['_data/hdv_', trial_name, '_n1_100_n2_100_tmp.mat']);
% hdv_parameter1;
% hdv_parameter2;
hdv_parameter = struct();

% Loop through the fields of the input struct
fields1 = fieldnames(hdv_parameter1);
for i = 1:numel(fields1)
    field = fields1{i};
    data1 = hdv_parameter1.(field);
    data2 = hdv_parameter2.(field);

    % Check if the field is a vector
    if length(data1) > 1
        hdv_parameter.(field) = [data1; data2];
    else
        % If the field is not a vector, include it in both structs
        hdv_parameter.(field) = data1;
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
n_ctr = 2 * n;      % number of state variables
m_ctr = m;            % number of input variables
p_ctr = n + m;  % number of output variables

for i_data = 1:data_number
    % Load trajectory data
    load(['./_data/trajectory_data_collection/', trial_name, '_data', num2str(i_data), '_T_',num2str(T),'_',num2str(Tini),'_',num2str(N),'_noiseLevel_',num2str(acel_noise),'_tmp.mat']);
    % ---------------------------------------
    %   Scenario initialization
    %----------------------------------------

    % There is one head vehicle at the very beginning
    S          = zeros(total_time_step, n + 1, 3);
    S(1,1,1)    = 0;
    for i = 2 : n + 1
        S(1,i,1) = S(1,i-1,1) - hdv_parameter.s_star(i-1);
    end
    S(1,:,2)    = v_star * ones(n + 1,1);

    % ------------------
    %  DeeP-LCC Formulation
    % ------------------
    Q_v         = weight_v*eye(n);              % penalty for velocity error
    Q_s         = weight_s*eye(p_ctr-n);        % penalty for spacing error
    Q           = blkdiag(Q_v,Q_s);             % penalty for trajectory error
    R           = weight_u*eye(m_ctr);          % penalty for control input

    u           = zeros(m_ctr,total_time_step); % control input
    x           = zeros(n_ctr,total_time_step); % state variables
    y           = zeros(p_ctr,total_time_step); % output variables
    pr_status   = zeros(total_time_step,1);     % problem status
    e           = zeros(1,total_time_step);     % external input

    % ------------------
    %  Separate Hankel Matrix
    % -----------------
    ud  = [ud1; ud2];
    ed  = ed1;
    yd  = [yd1; yd2];
    
    cd  = ed2 - yd1(end, :);
    
    U   = hankel_matrix(ud, Tini + N);
    Up  = U(1:Tini * m_ctr,:);
    Uf  = U((Tini * m_ctr + 1):end,:);

    E   = hankel_matrix(ed, Tini + N);
    Ep  = E(1:Tini,:);
    Ef  = E((Tini + 1):end,:);

    Y   = hankel_matrix(yd, Tini + N);
    Yp  = Y(1:Tini * p_ctr,:);
    Yf  = Y((Tini * p_ctr + 1):end,:);
    
    C   = hankel_matrix(yd, Tini + N);


    % ------------------
    %  Reference trajectory
    % ------------------
    r      = zeros(p_ctr, total_time_step + N);            % stabilization

    % ---------------------------------------------------------------------
    %   Simulation starts here
    %----------------------------------------------------------------------
    tic

    % ------------------
    %  Initial trajectory
    % ------------------
    uini = zeros(m_ctr,Tini);
    eini = zeros(1,Tini);
    yini = zeros(p_ctr,Tini);

    for k = 1:Tini-1
        % Update acceleration for flow 1
        acel = HDV_dynamics(S(k,:,:),hdv_parameter) ...
             - acel_noise + 2*acel_noise*rand(n,1);

        S(k,1,3)           = 0;             % the head vehicle
        S(k,2:end,3)       = acel;          % all the vehicles using HDV model
        S(k,Omega_c+1,3)  = uini(:,k);      % the CAV

        S(k+1,:,2) = S(k,:,2) + Tstep*S(k,:,3);
        S(k+1,1,2) = eini(k) + v_star;      % the velocity of the head vehicle
        S(k+1,:,1) = S(k,:,1) + Tstep*S(k,:,2);

        yini(:,k) = measure_mixed_traffic(S(k,2:end,2),S(k,:,1),ID,v_star,s_star,measure_type);

    end

    k_end = k+1;
    yini(:,k_end) = measure_mixed_traffic(S(k_end,2:end,2),S(k_end,:,1),ID,v_star,s_star,measure_type);

    u(:,1:Tini) = uini;
    e(:,1:Tini) = eini;
    y(:,1:Tini) = yini;


    % ------------------
    %  simulation starts here
    % ------------------
    for k = Tini:total_time_step-1
        % Update acceleration
        tic
        acel = HDV_dynamics(S(k,:,:),hdv_parameter) ...
             - acel_noise + 2*acel_noise*rand(n,1);

        S(k,2:end,3) = acel;     % all the vehicles using HDV model

        [u_opt,y_opt,pr] = DeeP_LCC_2(Up,Yp,Uf,Yf,Ep,Ef,C,uini,yini,eini,Q,R,r(:,k:k+N-1),...
                           lambda_g,lambda_y,u_limit,s_limit);
        toc
        computation_time(k) = toc;

        % One-step formulation
        u(:,k) = u_opt(1:m_ctr,1);
        % Update accleration for the CAV
        S(k,Omega_c+1,3)   = u(:,k);
        % Judge whether SD system commands to brake
        brake_vehicle_ID = find(acel==dcel_max);                 % the vehicles that need to brake
        brake_cav_ID     = intersect(brake_vehicle_ID,Omega_c);  % the CAVs that need to brake
        if ~isempty(brake_cav_ID)
            S(k,brake_cav_ID+1,3) = dcel_max;
        end

        % Update state
        S(k+1,:,2) = S(k,:,2) + Tstep*S(k,:,3);

        % -------------
        % Perturbation for the head vehicle
        % -------------
        switch per_type
            case 1
                S(k+1,1,2) = v_star + sine_amp*sin(2*pi/(10/Tstep)*(k-Tini));
                S(k+1,:,1) = S(k,:,1) + Tstep*S(k,:,2);
            case 2
                if (k-Tini)*Tstep <= brake_amp/5
                    S(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= brake_amp/5+5
                    S(k+1,1,3) = 1;
                else
                    S(k+1,1,3) = 0;
                end
                S(k+1,:,2) = S(k,:,2) + Tstep*S(k,:,3);
                S(k+1,:,1) = S(k,:,1) + Tstep*S(k,:,2);
            case 3
                if (k-Tini)*Tstep <= brake_amp/5
                    S(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= brake_amp/5+3
                    S(k+1,1,3) = 0;
                elseif (k-Tini)*Tstep <= brake_amp/5+3+5
                    S(k+1,1,3) = 1;
                else
                    S(k+1,1,3) = 0;
                end
                S(k+1,:,2) = S(k,:,2) + Tstep*S(k,:,3);
                S(k+1,:,1) = S(k,:,1) + Tstep*S(k,:,2);
            case 4
                if (k-Tini)*Tstep <= 2
                    S(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= 2+5
                    S(k+1,1,3) = 2;
                end
                S(k+1,:,2) = S(k,:,2) + Tstep*S(k,:,3);
                S(k+1,:,1) = S(k,:,1) + Tstep*S(k,:,2);
            case 5
                if (k-Tini)*Tstep <= brake_amp/2
                    S(k,4,3) = -2;
                end
                S(k+1,:,2) = S(k,:,2) + Tstep*S(k,:,3);
                S(k+1,:,1) = S(k,:,1) + Tstep*S(k,:,2);
        end

        % Record output
        y(:,k) = measure_mixed_traffic(S(k,2:end,2),S(k,:,1),ID,v_star,s_star,measure_type);
        e(k)   = S(k,1,2) - v_star;

        % update past data in control process
        uini = u(:,k-Tini+1:k);
        yini = y(:,k-Tini+1:k);
        eini = S(k-Tini+1:k,1,2) - v_star;

        fprintf('Simulation number: %d  |  process... %2.2f%% \n',i_data,k/total_time_step*100);

        str=['Num:', num2str(i_data),' | Processing: ',num2str(k/total_time_step*100),'%'];
        waitbar(k/total_time_step,h_wait,str);
    end

    k_end = k+1;
    y(:,k_end) = measure_mixed_traffic(S(k_end,2:end,2),S(k_end,:,1),ID,v_star,s_star,measure_type);

    tsim = toc;

    fprintf('Simulation ends at %6.4f seconds \n', tsim);


    % -------------------------------------------------------------------------
    %   Plot Results
    %--------------------------------------------------------------------------

end
close(h_wait);
