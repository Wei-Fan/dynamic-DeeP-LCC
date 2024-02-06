% =========================================================================
%               Data collection for dynamic traffic flow
% =========================================================================

clc; close all; clear;
addpath('_fcn');
addpath('_data');

trial_name = 'one_flow';

% ------------------------------------------
% Simulation Parameters
% ------------------------------------------
% How many data sets to collect
data_total_number = 10;
h_wait = waitbar(0, 'please wait');

% Parameters in Simulation
total_time       = 40;              % Total Simulation Time
Tstep            = 0.05;            % Time Step
total_time_step  = total_time / Tstep;

% ------------------------------------------
% DeeP-LCC Parameters
% ------------------------------------------
% T       = 1200;      % length of data samples
T               = 3000;       % Second choise

Tini            = 20;        % length of past data
N               = 50;        % length of predicted horizon

weight_v        = 1;        % weight coefficient for velocity error
weight_s        = 0.5;      % weight coefficient for spacing error
weight_u        = 0.1;      % weight coefficient for control input

lambda_g        = 10;        % penalty on ||g||_2^2 in objective
lambda_y        = 1e4;      % penalty on ||sigma_y||_2^2 in objective


% ------------------------------------------
% Parameters in Mixed Traffic
% ------------------------------------------
n = 100;  % number of total vehicles
m = 10;   % number of CAVs
% Type for HDV car-following model
hdv_type        = 1;    % 1. OVM   2. IDM
% Uncertainty for HDV behavior
acel_noise      = 0.1;  % A white noise signal on HDV's original acceleration

ID = generate_mixed_traffic_flow(n, m);

Omega   = find(ID==0);
Omega_c = find(ID==1);

v_star      = 15;                   % Equilibrium velocity
s_star      = 20;                   % Equilibrium spacing for CAV

acel_max = 2;
dcel_max = -5;

save(['./_data/ID_',trial_name,'_n_',num2str(n),'_tmp.mat'],'ID','acel_noise','v_star','s_star');

% Random setup for HDV
hdv_parameter = generate_HDVs(n);

save(['./_data/hdv_',trial_name,'_n_',num2str(n),'_tmp.mat'],'hdv_type','hdv_parameter');

% ------------------
%  size in DeeP-LCC
% ------------------

n_ctr = 2 * n; % number of state variables
m_ctr = m;     % number of input variables
p_ctr = n + m; % number of output variables

% -------------------------------------------------------------------------
%   Scenario initialization
%--------------------------------------------------------------------------
for i_data = 1:data_total_number
    % -------------------------------------------------------------------------
    %   Scenario initialization
    %--------------------------------------------------------------------------

    % There are two identical head vehicle at the very beginning
    S           = zeros(total_time_step, n + 1, 3);
    S(1, 1, 1)  = 0;
    for i = 2 : n + 1
        S(1, i, 1) = S(1, i - 1, 1) - hdv_parameter.s_star(i - 1);
    end
    S(1, :, 2)  = v_star * ones(n + 1, 1);

    % -------------------------------------------------------------------------
    %   Data collection
    %--------------------------------------------------------------------------

    % ------------------
    %  persistently exciting input data
    % ------------------
    ud          = -1 + 2 * rand(m_ctr, T);
    ed          = -1 + 2 * rand(1, T);
    yd          = zeros(p_ctr, T);

    % ------------------
    %  generate output data
    % ------------------
    for k = 1:T-1
        % Update acceleration
        acel                = HDV_dynamics(S(k,:,:), hdv_parameter) ...
                             -acel_noise + 2 * acel_noise * rand(n,1);

        S(k,1,3)            = 0;         % the head vehicle
        S(k,2:end,3)        = acel;      % all the vehicles using HDV model
        S(k,Omega_c + 1,3) = ud(:,k);    % the CAVs

        S(k+1,:,2) = S(k,:,2) + Tstep * S(k,:,3);
        S(k+1,1,2) = ed(k) + v_star;   % the velocity of the head vehicle
        S(k+1,:,1) = S(k,:,1) + Tstep * S(k,:,2);

        yd(:,k) = measure_mixed_traffic(S(k,2:end,2), S(k,:,1), ID, v_star, s_star, 3);
    end
    k = k+1;
    yd(:,k) = measure_mixed_traffic(S(k,2:end,2), S(k,:,1), ID, v_star, s_star, 3);

    % ------------------
    %  data Hankel matrices
    % ------------------
    % For DeeP-LCC 1
    U   = hankel_matrix(ud, Tini + N);
    Up  = U(1:Tini * m_ctr,:);
    Uf  = U((Tini * m_ctr + 1):end,:);

    E   = hankel_matrix(ed, Tini + N);
    Ep  = E(1:Tini,:);
    Ef  = E((Tini + 1):end,:);

    Y   = hankel_matrix(yd, Tini + N);
    Yp  = Y(1:Tini * p_ctr,:);
    Yf  = Y((Tini * p_ctr + 1):end,:);

    str=['Processing...',num2str(i_data/data_total_number*100),'%'];
    waitbar(i_data/data_total_number, h_wait, str);

    save(['./_data/trajectory_data_collection/',trial_name, '_data',num2str(i_data),'_T_',num2str(T),'_',num2str(Tini),'_',num2str(N),'_noiseLevel_',num2str(acel_noise),'_tmp.mat'],...
        'ud','ed','yd','T','Tini','N','Tstep');

end

close(h_wait);
