%
clear all; clc; close all;
addpath('path to data');

all_Stages_TS = load('all_STAGES_ts.mat');
for i = 1:37;
    stage_II_TS{1,i} = reshape(all_Stages_TS.Bold_EP_stages{1,1}(i,:,:),115,276);
end

for i = 1:19;
    stage_IIIa_TS{1,i} = reshape(all_Stages_TS.Bold_EP_stages{1,2}(i,:,:),115,276);
end

for i = 1:22;
    stage_IIIb_TS{1,i} = reshape(all_Stages_TS.Bold_EP_stages{1,3}(i,:,:),115,276);
end

for i = 1:9;
    stage_IIIc_TS{1,i} = reshape(all_Stages_TS.Bold_EP_stages{1,4}(i,:,:),115,276);
end

stage_IV_TS{1,1} = reshape(all_Stages_TS.Bold_EP_stages{1,5}(1,:,:),115,276);

SC_CNT = load('clean_cnt_sc.mat'); SC_CNT = SC_CNT.SC_CNT;

SC_EP = load('all_STAGES_ep_sc.mat'); SC_EP = SC_EP.SC_EP_stages; 

stage_II_SC = SC_EP{1,1};
stage_IIIa_SC = SC_EP{1,2};
stage_IIIb_SC = SC_EP{1,3};
stage_IIIc_SC = SC_EP{1,4};
stage_IV_SC = SC_EP{1,5};

Bold_CNT = load('cnt_ts_scale.mat'); Bold_CNT = Bold_CNT.Bold_CNT;
for i = 1:128
    CNT{1,i} = reshape(Bold_CNT(i,:,:),115,276);
end

% Bold_EP = load('ep_ts_scale1_clearNaN.mat'); Bold_EP = Bold_EP.Bold_EP; 
%%
csvwrite('FC_control_em.csv', CNT(1)); 
csvwrite('FC_II_em.csv', stage_II_TS(1));
csvwrite('FC_IIIb_em.csv', stage_IIIb_TS(1));
csvwrite('FC_IIIc_em.csv', stage_IIIc_TS(1));

%%
tseries = stage_IIIc_TS; 
Cnew = mean(stage_IIIc_SC, 3); % or SC_EP{1,n}

C = Cnew/max(max(Cnew))*0.2;

nNodes=size(C,1);
N = nNodes;
Isubdiag = find(tril(ones(N),-1));


% Parameters of the data
TR=2;  % Repetition Time (seconds)



% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.04;                    % lowpass frequency of filter (Hz) Original: 0.008
fhi = 0.07;                    % highpass; ORIGINAL 0.08
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter


% load ts compute FC, KoP and SF HighOrder and then the node freq

    xs = tseries; % or xs=tseries(:,cond);
    NSUB=size(find(~cellfun(@isempty,xs)),2);
    
 kk=1;   
for sub=1:NSUB
    %sub   
    Tmax = size(xs{sub},2);
    fce1=zeros(NSUB,N,N);
    clear signal_filt Phases ts;
    ts=zeros(N,Tmax);
    sub;

    ts= double(xs{sub}); % or ts= xs{sub};
    ts = ts(:,1:Tmax);
    Isubdiag = find(tril(ones(N),-1));
    
    for seed=1:N
        ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
        Xanalytic = hilbert(demean(signal_filt(seed,:)));
        Phases(seed,:) = angle(Xanalytic);
    end
    FC_spatime(:,:,1)=corr(signal_filt',signal_filt','rows','complete');
    FC_spatime(:,:,2)=corr(signal_filt(:,1:end-1)',signal_filt(:,2:end)','rows','complete');
    FC_spatime(:,:,3)=corr(signal_filt(:,1:end-2)',signal_filt(:,3:end)','rows','complete');
    FC_spatime(:,:,4)=corr(signal_filt(:,1:end-3)',signal_filt(:,4:end)','rows','complete');
    FC_spatime(:,:,5)=corr(signal_filt(:,1:end-4)',signal_filt(:,5:end)','rows','complete');
    FC_spatime(:,:,6)=corr(signal_filt(:,1:end-5)',signal_filt(:,6:end)','rows','complete');
    FC_spatime(:,:,7)=corr(signal_filt(:,1:end-6)',signal_filt(:,7:end)','rows','complete');
    
    GBC_spatime_empsub(:,:,sub)=squeeze(nanmean(FC_spatime,2));
        
 % FCD      
    ii2=1;
 
    for t=1:15:Tmax-30
        jj2=1;
        cc=corrcoef((signal_filt(:,t:t+30))');
        for t2=1:15:Tmax-30
            cc2=corrcoef((signal_filt(:,t2:t2+30))');
            ca=corrcoef(cc(Isubdiag),cc2(Isubdiag));
            if jj2>ii2
                cotsamplingdata(kk)=ca(2);   %% this accumulate all elements of the FCD empirical
                kk=kk+1;
            end
            jj2=jj2+1;
        end
        ii2=ii2+1;
    end
    
    KoP=abs(nansum(complex(cos(Phases),sin(Phases)),1))/N;
    KuramotoOrderParameter2(sub)=mean(KoP);
    FCemp2(sub,:,:)=corrcoef(signal_filt');
    
    [Ns, Tmaxred]=size(signal_filt);
    TT=Tmaxred;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    nfreqs=length(freq);
    clear PowSpect;
    for seed=1:N
        pw = abs(fft(signal_filt(seed,:)));
        PowSpect(:,seed,sub) = pw(1:floor(TT/2)).^2/(TT/TR);
    end
end

Power_Areas=squeeze(mean(PowSpect,3));
for seed=1:N
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
end

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);

obs.FCemp=squeeze(mean(FCemp2,1)); % or obs.FCemp=squeeze(mean(FCemp2));
obs.KuramotoOrderParameter=squeeze(mean(KuramotoOrderParameter2));

obs.GBCspatime = squeeze(nanmean(GBC_spatime_empsub(:,:,:),3));
obs.cotsam = cotsamplingdata;
obs.f_diff = f_diff;
obs.SC = C;

group = xlsread('grouping_p.xlsx');
prior = 1;

Cfgpob = 10;
Cfggens = 150;

% some modular parameters
% the number of parameters that will optimizeobs
Cfg.dimpar = 3;  % if 1-> only a_i, if 2--> ai and G full; if 3-> a_i and G_i, if 3-> a_i, G_i and b_i
% the target funcion of the Genetic algorithm

Cfg.fitting = 1;
Cfg.verbose = 0; % 0 is faster and display in command window with 1, is lower and create a plot
Cfg.parallel=1; % 1 use parallel process, otherwise dont

%Numbers of subjects in each simulation 
Cfg.NsubSim = 1; %cda simulación se hace 5 veces
Tmaxsim = Tmax;
 
tiradas_Max = 20 
for i=1:tiradas_Max
    istring=num2str(i);
    
    fprintf(['Tirada #', istring ,' de ' , num2str(tiradas_Max),'\n']);
    
    [FC_master(:,:,i),Rta(:,i),RtG(:,i),RtB(:,i),out(i)]= hopf_genetic_AD(C,obs,nNodes,Cfgpob,Cfggens,Tmaxsim,group, obs.f_diff,TR,Cfg);
    fprintf(['1-SSIM alcanzada :', num2str(out(i).ssimfinal),'\n']); 
    fprintf(['log Frob alcanzada :', num2str(out(i).fc_dist),'\n']); 
end
cond = 1;
save(sprintf('simulation_stage_IIIc_CORR_dimpar3_new_filt10_%03d_%03d.mat',cond,prior), 'FC_master', 'Rta', 'RtG', 'out', 'obs', 'Cfg')