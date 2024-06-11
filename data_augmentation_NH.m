% We take the smallest ssim = best fit

clear all;clc;close all;

addpath('path to data');
addpath("path to results of run_HOPF_genetic_AD_atrophy.m");

SC_CNT = load('control_SC.mat'); SC_CNT = SC_CNT.SC_CNT;
SC_EP = load('EP_SC.mat'); SC_EP = SC_EP.SC_EP_stages;  

% different stage classification
stage_II_SC = SC_EP{1,1};
stage_IIIa_SC = SC_EP{1,2};
stage_IIIb_SC = SC_EP{1,3};
stage_IIIc_SC = SC_EP{1,4};
stage_IV_SC = SC_EP{1,5};
n = 1000;
fprintf(['Starting data augmentation...','\n']);

%% CONTROLS

load('simulation_CNT.mat');
SC = mean(SC_CNT, 3);

[~,min_iterada] = min([out.ssimfinal]);
a = Rta(:,min_iterada); G = RtG(:,min_iterada);
beta = zeros(115,1);
w = obs.f_diff;
TR = 2;

C = (SC/max(max(SC)))*0.2;
Tmax = 200;

for i = 1:n;
    [xs]=hopf_supersub_G_int(a,G,beta,C,Tmax,w, TR);
    FCsim_cnt(i,:,:) = corr(xs);
    FCsim_vector_cnt(i,:) = [0,reshape(FCsim_cnt(i,:,:),[1,13225])]; 
end

fprintf(['Controls Complete','\n']);

%% Stage II

load('simulation_stage_II.mat');
SC = mean(stage_II_SC, 3);

[~,min_iterada] = min([out.ssimfinal]);
a = Rta(:,min_iterada); G = RtG(:,min_iterada);
beta = zeros(115,1);
w = obs.f_diff;
TR = 2;

C = (SC/max(max(SC)))*0.2;
Tmax = 200;

for i = 1:n;
    [xs]=hopf_supersub_G_int(a,G,beta,C,Tmax,w, TR);
    FCsim_stage_II(i,:,:) = corr(xs);
    FCsim_vector_stage_II(i,:) = [1,reshape(FCsim_stage_II(i,:,:),[1,13225])]; 
end

fprintf(['Stage II Complete','\n']);

%% Stage IIIb

load('simulation_stage_IIIb.mat');
SC = mean(stage_IIIb_SC, 3);

[~,min_iterada] = min([out.ssimfinal]);
a = Rta(:,min_iterada); G = RtG(:,min_iterada);
beta = zeros(115,1);
w = obs.f_diff;
TR = 2;

C = (SC/max(max(SC)))*0.2;
Tmax = 200;

for i = 1:(n/2);
    [xs]=hopf_supersub_G_int(a,G,beta,C,Tmax,w, TR);
    FCsim_stage_IIIb(i,:,:) = corr(xs);
    FCsim_vector_stage_IIIb(i,:) = [3,reshape(FCsim_stage_IIIb(i,:,:),[1,13225])]; 
end

fprintf(['Stage IIIb complete','\n']);

%% Stage IIIc

load('simulation_stage_IIIc.mat');
SC = mean(stage_IIIc_SC, 3);

[~,min_iterada] = min([out.ssimfinal]);
a = Rta(:,min_iterada); G = RtG(:,min_iterada);
beta = zeros(115,1);
w = obs.f_diff;
TR = 2;

C = (SC/max(max(SC)))*0.2;
Tmax = 200;

for i = 1:(n/2);
    [xs]=hopf_supersub_G_int(a,G,beta,C,Tmax,w, TR);
    FCsim_stage_IIIc(i,:,:) = corr(xs);
    FCsim_vector_stage_IIIc(i,:) = [4,reshape(FCsim_stage_IIIc(i,:,:),[1,13225])]; 
end

fprintf(['Stage IIIc Complete','\n']);

%%

dataset_FCsim = [FCsim_vector_cnt; FCsim_vector_stage_II; FCsim_vector_stage_IIIb; FCsim_vector_stage_IIIc]; 
r = randperm(size(dataset_FCsim,1)); dataset_FCsim = dataset_FCsim(r,:);
labels = dataset_FCsim(:,1); dataset_FCsim = dataset_FCsim(:,[2:13226]);
r_2 = randperm(size(labels,1)); labels2 = labels(r_2,:);
csvwrite('dataset_FCsim_1_NH.csv', dataset_FCsim); csvwrite('labels_1_NH.csv', labels2);

fprintf(['Data Augmentation Complete (Null Hypothesis)','\n']);
