% =========================================================================
%               Data collection for dynamic traffic flow
% =========================================================================

clc; close all; clear;
addpath('_fcn');
addpath('_data');

trial_name = 'two_flow';

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
T               = 600;       % Second choise

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
n1 = 100;  % number of total vehicles
m1 = 10;   % number of CAVs
n2 = 100;  % number of total vehicles
m2 = 10;   % number of CAVs
% Type for HDV car-following model
hdv_type        = 1;    % 1. OVM   2. IDM
% Uncertainty for HDV behavior
acel_noise      = 0.1;  % A white noise signal on HDV's original acceleration

ID1 = generate_mixed_traffic_flow(n1, m1);
ID2 = generate_mixed_traffic_flow(n2, m2);
Omega1   = find(ID1==0);
Omega1_c = find(ID1==1);
Omega2   = find(ID2==0);
Omega2_c = find(ID2==1);

v_star      = 15;                   % Equilibrium velocity
s_star      = 20;                   % Equilibrium spacing for CAV

acel_max = 2;
dcel_max = -5;

% Random setup for HDV
hdv_parameter1 = generate_HDVs(n1);
hdv_parameter2 = generate_HDVs(n2);


% ------------------
%  size in DeeP-LCC
% ------------------

n1_ctr = 2 * n1; % number of state variables
m1_ctr = m1;     % number of input variables
p1_ctr = n1 + m1; % number of output variables
n2_ctr = 2 * n2; % number of state variables
m2_ctr = m2;     % number of input variables
p2_ctr = n2 + m2; % number of output variables

% -------------------------------------------------------------------------
%   Scenario initialization
%--------------------------------------------------------------------------
for i_data = 1:data_total_number
    % -------------------------------------------------------------------------
    %   Scenario initialization
    %-------------------------------------------------------------------------- 
    
    % There are two identical head vehicle at the very beginning
    S1           = zeros(total_time_step, n1 + 1, 3);
    S1(1, 1, 1)  = 0;
    for i = 2 : n1 + 1
        S1(1, i, 1) = S1(1, i - 1, 1) - hdv_parameter1.s_star(i - 1);
    end
    S1(1, :, 2)  = v_star * ones(n1 + 1, 1);
    
    % There are two identical head vehicle at the very beginning
    S2           = zeros(total_time_step, n1 + 1, 3);
    S2(1, 1, 1)  = 0;
    for i = 2 : n1 + 1
        S2(1, i, 1) = S2(1, i - 1, 1) - hdv_parameter1.s_star(i - 1);
    end
    S2(1, :, 2)  = v_star * ones(n1 + 1, 1);

    % -------------------------------------------------------------------------
    %   Data collection
    %-------------------------------------------------------------------------- 
    
    % ------------------
    %  persistently exciting input data
    % ------------------
    ud1          = -1 + 2 * rand(m1_ctr, T);
    ed1          = -1 + 2 * rand(1, T);
    yd1          = zeros(p1_ctr, T);
    ud2          = -1 + 2 * rand(m2_ctr, T);
    ed2          = -1 + 2 * rand(1, T);
    yd2          = zeros(p2_ctr, T);
        
    % ------------------
    %  generate output data
    % ------------------
    for k = 1:T-1
        % Update acceleration
        acel1                = HDV_dynamics(S1(k,:,:), hdv_parameter1) ...
                             -acel_noise + 2 * acel_noise * rand(n1,1);
        
        S1(k,1,3)            = 0;         % the head vehicle
        S1(k,2:end,3)        = acel1;      % all the vehicles using HDV model
        S1(k,Omega1_c + 1,3) = ud1(:,k);   % the CAVs
        
        S1(k+1,:,2) = S1(k,:,2) + Tstep * S1(k,:,3);
        S1(k+1,1,2) = ed1(k) + v_star;   % the velocity of the head vehicle
        S1(k+1,:,1) = S1(k,:,1) + Tstep * S1(k,:,2);    
        
        yd1(:,k) = measure_mixed_traffic(S1(k,2:end,2), S1(k,:,1), ID1, v_star, s_star, 3);

        % Update acceleration
        acel2                = HDV_dynamics(S2(k,:,:), hdv_parameter2) ...
                             -acel_noise + 2 * acel_noise * rand(n2,1);
        
        S2(k,1,3)            = 0;         % the head vehicle
        S2(k,2:end,3)        = acel2;      % all the vehicles using HDV model
        S2(k,Omega2_c + 1,3) = ud2(:,k);   % the CAVs
        
        S2(k+1,:,2) = S2(k,:,2) + Tstep * S2(k,:,3);
        S2(k+1,1,2) = ed2(k) + v_star;   % the velocity of the head vehicle
        S2(k+1,:,1) = S2(k,:,1) + Tstep * S2(k,:,2);    
        
        yd2(:,k) = measure_mixed_traffic(S2(k,2:end,2), S2(k,:,1), ID2, v_star, s_star, 3);
    end
    k = k+1;
    yd1(:,k) = measure_mixed_traffic(S1(k,2:end,2), S1(k,:,1), ID1, v_star, s_star, 3);
    yd2(:,k) = measure_mixed_traffic(S2(k,2:end,2), S2(k,:,1), ID2, v_star, s_star, 3);
    
    
    % ------------------
    %  data Hankel matrices
    % ------------------
    % For DeeP-LCC 1
    U1   = hankel_matrix(ud1, Tini + N);
    U1p  = U1(1:Tini * m1_ctr,:);
    U1f  = U1((Tini * m1_ctr + 1):end,:);
    
    E1   = hankel_matrix(ed1, Tini + N);
    E1p  = E1(1:Tini,:);
    E1f  = E1((Tini + 1):end,:);
    
    Y1   = hankel_matrix(yd1, Tini + N);
    Y1p  = Y1(1:Tini * p1_ctr,:);
    Y1f  = Y1((Tini * p1_ctr + 1):end,:);

    % For DeeP-LCC 2
    U2   = hankel_matrix(ud2, Tini + N);
    U2p  = U2(1:Tini * m2_ctr,:);
    U2f  = U2((Tini * m2_ctr + 1):end,:);
    
    E2   = hankel_matrix(ed2, Tini + N);
    E2p  = E2(1:Tini,:);
    E2f  = E2((Tini + 1):end,:);
    
    Y2   = hankel_matrix(yd2, Tini + N);
    Y2p  = Y2(1:Tini * p2_ctr,:);
    Y2f  = Y2((Tini * p2_ctr + 1):end,:);
    
    
    str=['Processing...',num2str(i_data/data_total_number*100),'%'];
    waitbar(i_data/data_total_number, h_wait, str);

    save_data = {hdv_type,acel_noise,U1p,Y1p,U1f,Y1f,E1p,E1f,U2p,Y2p,U2f,Y2f,E2p,E2f,T,Tini,N,ID1,ID2,Tstep,v_star};
    save(['.\_data\trajectory_data_collection\',trial_name, '_data',num2str(i_data),'_T_',num2str(T),'_',num2str(Tini),'_',num2str(N),'_noiseLevel_',num2str(acel_noise),'_tmp.mat'],"save_data");

end

close(h_wait);