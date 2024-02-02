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

%TODO:!!!
for i_data = 1:data_number
    % Load trajectory data
    load(['./_data/trajectory_data_collection/','data',num2str(i_data),'_T_',num2str(T),'_',num2str(Tini),'_',num2str(N),'_noiseLevel_',num2str(acel_noise),'.mat']);
    
    % ---------------------------------------
    %   Scenario initialization
    %----------------------------------------
    
    % There is one head vehicle at the very beginning
    S1          = zeros(total_time_step, n1_vehicle + 1, 3);
    S(1,1,1)    = 0;
    for i = 2 : n1_vehicle + 1
        S1(1,i,1) = S1(1,i-1,1) - hdv_parameter1.s_star(i-1);
    end
    S1(1,:,2)    = v_star * ones(n1_vehicle + 1,1);
end
