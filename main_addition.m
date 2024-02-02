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
warning off;


% whether output data
output_bool = 1;

% whether mixed traffic flow
mix         = 1;

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

% Constraints
acel_max        = 2;
dcel_max        = -5;
spacing_max     = 40;
spacing_min     = 5;
u_limit         = [dcel_max, acel_max];
s_limit         = [spacing_min, spacing_max] - s_star;

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
