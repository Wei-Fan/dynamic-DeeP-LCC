% =========================================================================
%               Data collection for dynamic traffic flow
% =========================================================================

clc; close all; clear;
addpath('_fcn');
addpath('_data');

% How many data sets to collect
data_total_number = 100;
h_wait = waitbar(0, 'please wait');

% -------------------------------------------------------------------------
%   Parameter setup
% -------------------------------------------------------------------------

% Type for HDV car-following model
hdv_type        = 1;    % 1. OVM   2. IDM
% Uncertainty for HDV behavior
acel_noise      = 0.1;  % A white noise signal on HDV's original acceleration

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
load('ID_LargeScale.mat'); % record ID
ID_str      = num2str(ID);
ID_str(find(ID_str==' ')) = '';
n2          = 10;                   % number of vehicles from other lanes
ID1         = ID(1:(length(ID) - n2));
ID2         = ID((length(ID) - n2 + 1) : end);

pos1_cav    = find(ID1==1);          % position of CAVs
n1_vehicle  = length(ID) - n2;           % number of vehicles
n1_cav      = length(pos1_cav);      % number of CAVs
n1_hdv      = n1_vehicle - n1_cav;   % number of HDVs

pos2_cav    = find(ID2==1);          % position of CAVs
n2_vehicle  = n2;                    % number of vehicles
n2_cav      = length(pos2_cav);      % number of CAVs
n2_hdv      = n2_vehicle - n2_cav;   % number of HDVs

mix         = 1;                    % whether mixed traffic flow

v_star      = 15;                   % Equilibrium velocity
s_star      = 20;                   % Equilibrium spacing for CAV

acel_max = 2;
dcel_max = -5;

% Random setup for OVM
load(['_data/hdv_ovm_random_largescale.mat']);
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
            hdv_parameter1.(field) = data(1:(n - n2));
            hdv_parameter2.(field) = data((n - n2 + 1):end);
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

n1_ctr = 2 * n1_vehicle;    % number of state variables
m1_ctr = n1_cav;          % number of input variables
switch measure_type     % number of output variables
    case 1
        p1_ctr = n1_vehicle;
    case 2
        p1_ctr = 2 * n1_vehicle;
    case 3
        p1_ctr = n1_vehicle + n1_cav;
end

n2_ctr = 2 * n2_vehicle;    % number of state variables
m2_ctr = n2_cav;          % number of input variables
switch measure_type     % number of output variables
    case 1
        p2_ctr = n2_vehicle;
    case 2
        p2_ctr = 2 * n2_vehicle;
    case 3
        p2_ctr = n2_vehicle + n2_cav;
end

for i_data = 1:data_total_number
    % -------------------------------------------------------------------------
    %   Scenario initialization
    %-------------------------------------------------------------------------- 
    
    % There are two identical head vehicle at the very beginning
    S1           = zeros(total_time_step, n1_vehicle + 1, 3);
    S1(1, 1, 1)  = 0;
    for i = 2 : n1_vehicle + 1
        S1(1, i, 1) = S1(1, i - 1, 1) - hdv_parameter1.s_star(i - 1);
    end
    S1(1, :, 2)  = v_star * ones(n1_vehicle + 1, 1);
    
    S2           = zeros(total_time_step, n2_vehicle + 1, 3);
    S2(1, 1, 1)  = 0;
    for i = 2 : n2_vehicle + 1
        S2(1, i, 1) = S2(1, i - 1, 1) - hdv_parameter2.s_star(i - 1);
    end
    S2(1, :, 2)  = v_star * ones(n2_vehicle + 1, 1);
    
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
                             -acel_noise + 2 * acel_noise * rand(n1_vehicle,1);
        
        S1(k,1,3)            = 0;         % the head vehicle
        S1(k,2:end,3)        = acel1;      % all the vehicles using HDV model
        S1(k,pos1_cav + 1,3) = ud1(:,k);   % the CAVs
        
        S1(k+1,:,2) = S1(k,:,2) + Tstep * S1(k,:,3);
        S1(k+1,1,2) = ed1(k) + v_star;   % the velocity of the head vehicle
        S1(k+1,:,1) = S1(k,:,1) + Tstep * S1(k,:,2);    
        
        yd1(:,k) = measure_mixed_traffic(S1(k,2:end,2),S1(k,:,1),ID1,v_star,s_star,measure_type);
        
        acel2                = HDV_dynamics(S2(k,:,:), hdv_parameter2) ...
                             -acel_noise + 2 * acel_noise * rand(n2_vehicle,1);
        
        S2(k,1,3)            = 0;         % the head vehicle
        S2(k,2:end,3)        = acel2;      % all the vehicles using HDV model
        S2(k,pos2_cav + 1,3) = ud2(:,k);   % the CAVs
        
        S2(k+1,:,2) = S2(k,:,2) + Tstep * S2(k,:,3);
        S2(k+1,1,2) = ed2(k) + v_star;   % the velocity of the head vehicle
        S2(k+1,:,1) = S2(k,:,1) + Tstep * S2(k,:,2);    
        
        yd2(:,k) = measure_mixed_traffic(S2(k,2:end,2), S2(k,:,1), ID2, v_star, s_star, measure_type);
    
    
    end
    k = k+1;
    yd1(:,k) = measure_mixed_traffic(S1(k,2:end,2), S1(k,:,1), ID1, v_star, s_star, measure_type);
    yd2(:,k) = measure_mixed_traffic(S2(k,2:end,2), S2(k,:,1), ID2, v_star, s_star, measure_type);
    
    
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

    save_data = {hdv_type,acel_noise,U1p,Y1p,U1f,Y1f,E1p,E1f,T,Tini,N,ID1,ID2,Tstep,v_star,U2p,U2f,E2p,E2f,Y2p,Y2f};
    save(['.\_data\trajectory_data_collection\','data',num2str(i_data),'_T_',num2str(T),'_',num2str(Tini),'_',num2str(N),'_noiseLevel_',num2str(acel_noise),'.mat'],"save_data");

end

close(h_wait);
