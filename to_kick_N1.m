clear all;clc;

addpath('path to data');
addpath("path to results of run_HOPF_genetic_AD_atrophy");

load('simulation_stage_II.mat');
SC_EP = load('all_STAGES_ep_sc.mat'); SC_EP = SC_EP.SC_EP_stages;
stage_II_SC = SC_EP{1,1};
SC = mean(stage_II_SC, 3);

[~,min_iterada] = min([out.ssimfinal]);
a = Rta(:,min_iterada); G = RtG(:,min_iterada);

Cfg.nodos_kick = 115; 
Cfg.parallel=0;
Cfg.amplitud= 20;
Cfg.Repe = 10;
Cfg.filt.bpass=0.04;
Cfg.nSub=1;
Cfg.Tmax=1500;
Cfg.TRsec =2.;

% SC_best = sign(FC_emp_W_P).*SC;

[val_best, ~] = resamplingID(Cfg.Tmax, Cfg.TRsec,0.1);

%to load the information of each basal state
%load('adapt_corridatotal._sujetoN1.mat')

w = obs.f_diff;

SC = SC/max(max(SC))*0.2; %SC_best;
aes = a; %Rta_best;
G = G; %RtG_best;
% omega = obs.f_diff; %w_best;

[FC_sim_kick, S, nodo_kicked]=fun_sim_patadas(aes,G,w,Cfg,SC,val_best);

save kick_II_bona.mat FC_sim_kick S

fprintf(['Stage II complete','\n']);

%%

clear all;clc;

addpath('path to data');
addpath("path to results of run_HOPF_genetic_AD_atrophy");

load('simulation_stage_IIIb.mat');
SC_EP = load('all_STAGES_ep_sc.mat'); SC_EP = SC_EP.SC_EP_stages;
stage_IIIb_SC = SC_EP{1,3};
SC = mean(stage_IIIb_SC, 3);

[~,min_iterada] = min([out.ssimfinal]);
a = Rta(:,min_iterada); G = RtG(:,min_iterada);

Cfg.nodos_kick = 115; 
Cfg.parallel=0;
Cfg.amplitud= 20;
Cfg.Repe = 10;
Cfg.filt.bpass=0.04;
Cfg.nSub=1;
Cfg.Tmax=1500;
Cfg.TRsec =2.;

% SC_best = sign(FC_emp_W_P).*SC;

[val_best, ~] = resamplingID(Cfg.Tmax, Cfg.TRsec,0.1);

%to load the information of each basal state
%load('adapt_corridatotal._sujetoN1.mat')

w = obs.f_diff;

SC = SC/max(max(SC))*0.2; %SC_best;
aes = a; %Rta_best;
G = G; %RtG_best;
% omega = obs.f_diff; %w_best;

[FC_sim_kick, S, nodo_kicked]=fun_sim_patadas(aes,G,w,Cfg,SC,val_best);

save kick_IIIb_bona.mat FC_sim_kick S

fprintf(['Stage IIIb complete','\n']);

%%

clear all;clc;

addpath('path to data');
addpath("path to results of run_HOPF_genetic_AD_atrophy");

load('simulation_stage_IIIc.mat');
SC_EP = load('all_STAGES_ep_sc.mat'); SC_EP = SC_EP.SC_EP_stages;
stage_IIIc_SC = SC_EP{1,4};
SC = mean(stage_IIIc_SC, 3);

[~,min_iterada] = min([out.ssimfinal]);
a = Rta(:,min_iterada); G = RtG(:,min_iterada);

Cfg.nodos_kick = 115; 
Cfg.parallel=0;
Cfg.amplitud= 20;
Cfg.Repe = 10;
Cfg.filt.bpass=0.04;
Cfg.nSub=1;
Cfg.Tmax=1500;
Cfg.TRsec =2.;

% SC_best = sign(FC_emp_W_P).*SC;

[val_best, ~] = resamplingID(Cfg.Tmax, Cfg.TRsec,0.1);

%to load the information of each basal state
%load('adapt_corridatotal._sujetoN1.mat')

w = obs.f_diff;

SC = SC/max(max(SC))*0.2; %SC_best;
aes = a; %Rta_best;
G = G; %RtG_best;
% omega = obs.f_diff; %w_best;

[FC_sim_kick, S, nodo_kicked]=fun_sim_patadas(aes,G,w,Cfg,SC,val_best);

save kick_IIIc_bona.mat FC_sim_kick S

fprintf(['Stage IIIc complete','\n']);