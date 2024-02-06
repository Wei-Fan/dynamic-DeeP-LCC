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
per_type            = 0;    % 0. no perturbation 1. sinuoid perturbation 2. small brake perturbation
                            % 3. large brake perturbation 4. larger brake perturbation
                            % 5. Perturbation on a vehicle in the middle of the platoon
sine_amp            = 4; % amplitidue of sinuoid perturbation
brake_amp           = 5; % brake amplitude of brake perturbation

constraint_bool     = 1; % Whether there exist constraints

% Type for HDV car-following model
hdv_type            = 1;    % 1. OVM   2. IDM
% Uncertainty for HDV behavior
acel_noise          = 0.1;  % A white noise signal on HDV's original acceleration

% Parameters in Simulation
total_time          = 100;              % Total Simulation Time
Tstep               = 0.1;              % Time Step
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
Omega_c     = find(ID==1);

% Generate the first random integer
firstInteger = randi([2, m_total - 1], 1, 1);

% Initialize the second integer
secondInteger = firstInteger;

% Ensure the second integer is different from the first
while secondInteger == firstInteger
    secondInteger = randi([2, m_total], 1, 1);
end

% Sort the integers
m_separate = sort([firstInteger, secondInteger]);

m_separate  = sort(randi([2, m_total], 1, 2));
n_separate  = Omega_c(m_separate);

ID1         = ID(1:n_separate(1) - 1);
IDb         = ID(n_separate(1) : n_separate(2) - 1);
ID2         = ID(n_separate(2) : end);

Omega1_c    = find(ID1==1);          % position of CAVs
n1          = length(ID1);           % number of vehicles
m1          = length(Omega1_c);      % number of CAVs

Omegab_c    = find(IDb==1);          % position of CAVs
nb          = length(IDb);           % number of vehicles
mb          = length(Omegab_c);      % number of CAVs

Omega2_c    = find(ID2==1);          % position of CAVs
n2          = length(ID2);           % number of vehicles
m2          = length(Omega2_c);      % number of CAVs

IDa         = [ID1; ID2];
Omegaa_c    = find(IDa==1);          % position of CAVs
na          = length(IDa);           % number of vehicles
ma          = length(Omegaa_c);      % number of CAVs

% Constraints
acel_max        = 2;
dcel_max        = -5;
spacing_max     = 40;
spacing_min     = 5;
u_limit         = [dcel_max, acel_max];
s_limit         = [spacing_min, spacing_max] - s_star;

% Random setup for OVM
load(['_data/hdv_', trial_name, '_n_100_tmp.mat']);
% Initialize empty structs
hdv_parametera = struct();
hdv_parameterb = struct();

% Loop through the fields of the input struct
fields = fieldnames(hdv_parameter);
for i = 1:numel(fields)
    field = fields{i};
    data = hdv_parameter.(field);

    % Check if the field is a vector
    if length(data) > 1
        hdv_parametera.(field) = [data(1:n1);data((n1 + nb + 1): end)];
        hdv_parameterb.(field) = data((n1 + 1):(n1 + nb));
    else
        % If the field is not a vector, include it in both structs
        hdv_parametera.(field) = data;
        hdv_parameterb.(field) = data;
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

na_ctr = 2 * na;      % number of state variables
ma_ctr = ma;          % number of input variables
pa_ctr = na + ma;  % number of output variables

nb_ctr = 2 * nb;      % number of state variables
mb_ctr = mb;          % number of input variables
pb_ctr = nb + mb;     % number of output variables


for i_data = 1:data_number
    % Load trajectory data
    load(['./_data/trajectory_data_collection/', trial_name, '_data', num2str(i_data), '_T_',num2str(T),'_',num2str(Tini),'_',num2str(N),'_noiseLevel_',num2str(acel_noise),'_tmp.mat']);
    % ---------------------------------------
    %   Scenario initialization
    %----------------------------------------

    % There is one head vehicle at the very beginning
    Sa          = zeros(total_time_step, na + 1, 3);
    Sa(1,1,1)    = 0;
    for i = 2 : na + 1
        Sa(1,i,1) = Sa(1,i-1,1) - hdv_parametera.s_star(i-1);
    end
    Sa(1,:,2)    = v_star * ones(na + 1,1);

    Sb          = zeros(total_time_step, nb + 1, 3);
    Sb(1,1,1)    = 0;
    for i = 2 : nb + 1
        Sb(1,i,1) = Sb(1,i-1,1) - hdv_parameterb.s_star(i-1);
    end
    Sb(1,:,2)    = v_star * ones(nb + 1,1);

    % ------------------
    %  DeeP-LCC Formulation
    % ------------------
    Qa_v         = weight_v*eye(na);              % penalty for velocity error
    Qa_s         = weight_s*eye(pa_ctr-na);       % penalty for spacing error
    Qa           = blkdiag(Qa_v,Qa_s);            % penalty for trajectory error
    Ra           = weight_u*eye(ma_ctr);          % penalty for control input

    ua           = zeros(ma_ctr,total_time_step); % control input
    xa           = zeros(na_ctr,total_time_step); % state variables
    ya           = zeros(pa_ctr,total_time_step); % output variables
    pr_statusa   = zeros(total_time_step,1);      % problem status
    ea           = zeros(1,total_time_step);      % external input

    Qb_v         = weight_v*eye(nb);              % penalty for velocity error
    Qb_s         = weight_s*eye(pb_ctr-nb);       % penalty for spacing error
    Qb           = blkdiag(Qb_v,Qb_s);            % penalty for trajectory error
    Rb           = weight_u*eye(mb_ctr);          % penalty for control input

    ub           = zeros(mb_ctr,total_time_step); % control input
    xb           = zeros(nb_ctr,total_time_step); % state variables
    yb           = zeros(pb_ctr,total_time_step); % output variables
    pr_statusb   = zeros(total_time_step,1);      % problem status
    eb           = zeros(1,total_time_step);      % external input

    % ------------------
    %  Separate Hankel Matrix
    % ------------------
    uda = [ud(1:m1, :);ud(m1 + mb + 1: end, :)];
    udb = ud(m1 + 1:m1 + mb, :);
    eda = ed;
    edb = yd(n1, :);
    yda = [yd(1:n1, :); yd(n1+nb+1:n1+nb+n2, :);...
           yd(n1+nb+n2+1:n1+nb+n2+m1, :); yd(n1+nb+n2+m1+mb+1:end, :)];
    ydb = [yd(n1+1:n1+nb, :); yd(n1+nb+n2+m1+1:n1+nb+n2+m1+mb, :)];
    
    cd   = yd(n1+nb, :) - yd(n1, :);

    Ua   = hankel_matrix(uda, Tini + N);
    Uap  = Ua(1:Tini * ma_ctr,:);
    Uaf  = Ua((Tini * ma_ctr + 1):end,:);

    Ea   = hankel_matrix(eda, Tini + N);
    Eap  = Ea(1:Tini,:);
    Eaf  = Ea((Tini + 1):end,:);

    Ya   = hankel_matrix(yda, Tini + N);
    Yap  = Ya(1:Tini * pa_ctr,:);
    Yaf  = Ya((Tini * pa_ctr + 1):end,:);
    
    C   = hankel_matrix(cd, Tini + N);

    Ub   = hankel_matrix(udb, Tini + N);
    Ubp  = Ub(1:Tini * mb_ctr,:);
    Ubf  = Ub((Tini * mb_ctr + 1):end,:);

    Eb   = hankel_matrix(edb, Tini + N);
    Ebp  = Eb(1:Tini,:);
    Ebf  = Eb((Tini + 1):end,:);

    Yb   = hankel_matrix(ydb, Tini + N);
    Ybp  = Yb(1:Tini * pb_ctr,:);
    Ybf  = Yb((Tini * pb_ctr + 1):end,:);

    % ------------------
    %  Reference trajectory
    % ------------------
    ra       = zeros(pa_ctr, total_time_step + N);            % stabilization
    rb       = zeros(pb_ctr, total_time_step + N);            % stabilization

    % ---------------------------------------------------------------------
    %   Simulation starts here
    %----------------------------------------------------------------------
    tic

    % ------------------
    %  Initial trajectory
    % ------------------
    uaini = zeros(ma_ctr,Tini);
    eaini = zeros(1,Tini);
    yaini = zeros(pa_ctr,Tini);

    ubini = zeros(mb_ctr,Tini);
    ebini = zeros(1,Tini);
    ybini = zeros(pb_ctr,Tini);

    for k = 1:Tini-1
        % Update acceleration for flow 1
        acela = HDV_dynamics(Sa(k,:,:),hdv_parametera) ...
              - acel_noise + 2*acel_noise*rand(na,1);

        Sa(k,1,3)           = 0;               % the head vehicle
        Sa(k,2:end,3)       = acela;           % all the vehicles using HDV model
        Sa(k,Omegaa_c+1,3)  = uaini(:,k);      % the CAV

        Sa(k+1,:,2) = Sa(k,:,2) + Tstep*Sa(k,:,3);
        Sa(k+1,1,2) = eaini(k) + v_star;       % the velocity of the head vehicle
        Sa(k+1,:,1) = Sa(k,:,1) + Tstep*Sa(k,:,2);

        yaini(:,k) = measure_mixed_traffic(Sa(k,2:end,2),Sa(k,:,1),IDa,v_star,s_star,measure_type);

        % Update acceleration for flow 2
        acelb = HDV_dynamics(Sb(k,:,:),hdv_parameterb) ...
              - acel_noise + 2*acel_noise*rand(nb,1);

        Sb(k,1,3)           = 0;               % the head vehicle
        Sb(k,2:end,3)       = acelb;           % all the vehicles using HDV model
        Sb(k,Omegab_c+1,3)  = ubini(:,k);      % the CAV

        Sb(k+1,:,2) = Sb(k,:,2) + Tstep*Sb(k,:,3);
        Sb(k+1,1,2) = ebini(k) + v_star;       % the velocity of the head vehicle
        Sb(k+1,:,1) = Sb(k,:,1) + Tstep*Sb(k,:,2);

        ybini(:,k) = measure_mixed_traffic(Sb(k,2:end,2),Sb(k,:,1),IDb,v_star,s_star,measure_type);
    end

    k_end = k+1;
    yaini(:,k_end) = measure_mixed_traffic(Sa(k_end,2:end,2),Sa(k_end,:,1),IDa,v_star,s_star,measure_type);
    ybini(:,k_end) = measure_mixed_traffic(Sb(k_end,2:end,2),Sb(k_end,:,1),IDb,v_star,s_star,measure_type);

    ua(:,1:Tini) = uaini;
    ea(:,1:Tini) = eaini;
    ya(:,1:Tini) = yaini;

    ub(:,1:Tini) = ubini;
    eb(:,1:Tini) = ebini;
    yb(:,1:Tini) = ybini;

    % ------------------
    %  simulation starts here
    % ------------------
    for k = Tini:total_time_step-1
        % Update acceleration
        tic
        acela = HDV_dynamics(Sa(k,:,:),hdv_parametera) ...
              - acel_noise + 2*acel_noise*rand(na,1);
        acelb = HDV_dynamics(Sb(k,:,:),hdv_parameterb) ...
              - acel_noise + 2*acel_noise*rand(nb,1);

        Sa(k,2:end,3) = acela;     % all the vehicles using HDV model
        Sb(k,2:end,3) = acelb;     % all the vehicles using HDV model

        [ua_opt,ya_opt,pra] = DeeP_LCC_2(Uap,Yap,Uaf,Yaf,Eap,Eaf, C,uaini,yaini,eaini,Qa,Ra,ra(:,k:k+N-1),...
                           lambda_g,lambda_y,u_limit,s_limit);
        [ub_opt,yb_opt,prb] = DeeP_LCC(Ubp,Ybp,Ubf,Ybf,Ebp,Ebf,ubini,ybini,ebini,Qb,Rb,rb(:,k:k+N-1),...
                           lambda_g,lambda_y,u_limit,s_limit);
        toc
        computation_time(k) = toc;

        % One-step formulation
        ua(:,k) = ua_opt(1:ma_ctr,1);
        ub(:,k) = ub_opt(1:mb_ctr,1);
        % Update accleration for the CAV
        Sa(k,Omegaa_c+1,3)   = ua(:,k);
        Sb(k,Omegab_c+1,3)   = ub(:,k);
        % Judge whether SD system commands to brake
        brake_vehicle_ID1a = find(acela==dcel_max);                 % the vehicles that need to brake
        brake_cav_IDa     = intersect(brake_vehicle_IDa,Omegaa_c); % the CAVs that need to brake
        if ~isempty(brake_cav_IDa)
            Sa(k,brake_cav_IDa+1,3) = dcel_max;
        end

        brake_vehicle_IDb = find(acelb==dcel_max);                 % the vehicles that need to brake
        brake_cav_IDb     = intersect(brake_vehicle_IDb,Omegab_c); % the CAVs that need to brake
        if ~isempty(brake_cav_IDb)
            Sb(k,brake_cav_IDb+1,3) = dcel_max;
        end

        % Update state
        Sa(k+1,:,2) = Sa(k,:,2) + Tstep*Sa(k,:,3);
        Sb(k+1,:,2) = Sb(k,:,2) + Tstep*Sb(k,:,3);

        % -------------
        % Perturbation for the head vehicle
        % -------------
        switch per_type
            case 0
                Sa(k+1,1,2) = v_star;
                Sa(k+1,:,1) = Sa(k,:,1) + Tstep*Sa(k,:,2);

                Sb(k+1,1,2) = v_star;
                Sb(k+1,:,1) = Sb(k,:,1) + Tstep*Sb(k,:,2);
            case 1
                Sa(k+1,1,2) = v_star + sine_amp*sin(2*pi/(10/Tstep)*(k-Tini));
                Sa(k+1,:,1) = Sa(k,:,1) + Tstep*Sa(k,:,2);

                Sb(k+1,1,2) = v_star + sine_amp*sin(2*pi/(10/Tstep)*(k-Tini));
                Sb(k+1,:,1) = Sb(k,:,1) + Tstep*Sb(k,:,2);
            case 2
                if (k-Tini)*Tstep <= brake_amp/5
                    Sa(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= brake_amp/5+5
                    Sa(k+1,1,3) = 1;
                else
                    Sa(k+1,1,3) = 0;
                end
                Sa(k+1,:,2) = Sa(k,:,2) + Tstep*Sa(k,:,3);
                Sa(k+1,:,1) = Sa(k,:,1) + Tstep*Sa(k,:,2);

                if (k-Tini)*Tstep <= brake_amp/5
                    Sb(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= brake_amp/5+5
                    Sb(k+1,1,3) = 1;
                else
                    Sb(k+1,1,3) = 0;
                end
                Sb(k+1,:,2) = Sb(k,:,2) + Tstep*Sb(k,:,3);
                Sb(k+1,:,1) = Sb(k,:,1) + Tstep*Sb(k,:,2);
            case 3
                if (k-Tini)*Tstep <= brake_amp/5
                    Sa(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= brake_amp/5+3
                    Sa(k+1,1,3) = 0;
                elseif (k-Tini)*Tstep <= brake_amp/5+3+5
                    Sa(k+1,1,3) = 1;
                else
                    Sa(k+1,1,3) = 0;
                end
                Sa(k+1,:,2) = Sa(k,:,2) + Tstep*Sa(k,:,3);
                Sa(k+1,:,1) = Sa(k,:,1) + Tstep*Sa(k,:,2);

                if (k-Tini)*Tstep <= brake_amp/5
                    Sb(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= brake_amp/5+3
                    Sb(k+1,1,3) = 0;
                elseif (k-Tini)*Tstep <= brake_amp/5+3+5
                    Sb(k+1,1,3) = 1;
                else
                    Sb(k+1,1,3) = 0;
                end
                Sb(k+1,:,2) = Sb(k,:,2) + Tstep*Sb(k,:,3);
                Sb(k+1,:,1) = Sb(k,:,1) + Tstep*Sb(k,:,2);
            case 4
                if (k-Tini)*Tstep <= 2
                    Sa(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= 2+5
                    Sa(k+1,1,3) = 2;
                end
                Sa(k+1,:,2) = Sa(k,:,2) + Tstep*Sa(k,:,3);
                Sa(k+1,:,1) = Sa(k,:,1) + Tstep*Sa(k,:,2);

                if (k-Tini)*Tstep <= 2
                    Sb(k+1,1,3) = -5;
                elseif (k-Tini)*Tstep <= 2+5
                    Sb(k+1,1,3) = 2;
                end
                Sb(k+1,:,2) = Sb(k,:,2) + Tstep*Sb(k,:,3);
                Sb(k+1,:,1) = Sb(k,:,1) + Tstep*Sb(k,:,2);
            case 5
                if (k-Tini)*Tstep <= brake_amp/2
                    Sa(k,4,3) = -2;
                end
                Sa(k+1,:,2) = Sa(k,:,2) + Tstep*Sa(k,:,3);
                Sa(k+1,:,1) = Sa(k,:,1) + Tstep*Sa(k,:,2);

                if (k-Tini)*Tstep <= brake_amp/2
                    Sb(k,4,3) = -2;
                end
                Sb(k+1,:,2) = Sb(k,:,2) + Tstep*Sb(k,:,3);
                Sb(k+1,:,1) = Sb(k,:,1) + Tstep*Sb(k,:,2);
        end

        % Record output
        ya(:,k) = measure_mixed_traffic(Sa(k,2:end,2),Sa(k,:,1),IDa,v_star,s_star,measure_type);
        ea(k)   = Sa(k,1,2) - v_star;
        yb(:,k) = measure_mixed_traffic(Sb(k,2:end,2),Sb(k,:,1),IDb,v_star,s_star,measure_type);
        eb(k)   = Sb(k,1,2) - v_star;

        % update past data in control process
        uaini = ua(:,k-Tini+1:k);
        yaini = ya(:,k-Tini+1:k);
        eaini = Sa(k-Tini+1:k,1,2) - v_star;
        ubini = ub(:,k-Tini+1:k);
        ybini = yb(:,k-Tini+1:k);
        ebini = Sb(k-Tini+1:k,1,2) - v_star;

        fprintf('Simulation number: %d  |  process... %2.2f%% \n',i_data,k/total_time_step*100);

        str=['Num:', num2str(i_data),' | Processing: ',num2str(k/total_time_step*100),'%'];
        waitbar(k/total_time_step,h_wait,str);
    end

    k_end = k+1;
    ya(:,k_end) = measure_mixed_traffic(Sa(k_end,2:end,2),Sa(k_end,:,1),IDa,v_star,s_star,measure_type);
    yb(:,k_end) = measure_mixed_traffic(Sb(k_end,2:end,2),Sb(k_end,:,1),IDb,v_star,s_star,measure_type);

    tsim = toc;

    fprintf('Simulation ends at %6.4f seconds \n', tsim);


    % -------------------------------------------------------------------------
    %   Plot Results
    %--------------------------------------------------------------------------

end
close(h_wait);
