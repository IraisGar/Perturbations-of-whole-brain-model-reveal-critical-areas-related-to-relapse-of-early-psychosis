% to iterate
function [FC_sim_kick, S, nodo_kicked]=fun_sim_patadas(a,G,w,Cfg,SC,val);

Tmax=1500;
TR =2.;
beta = zeros(115,1);
nodo_kick = zeros(115,1);
if Cfg.parallel==1
    for ii=1:Cfg.nodos_kick; % goes though pairs of kicks
        nodo_kick = zeros(115,1);
        if Cfg.nodos_kick==115
            nodo_kick(ii)=1;
        else
            nodo_kick(ii)=1;
            nodo_kick(116-ii)=1 % change depending on the number of nodes
        end
        ii
        for j=1:Cfg.amplitud; % swips the value of the amplitude of (S)
            S(j)= 0+0.1*(j-1);
            parfor i=1:Cfg.Repe % iteration for every amplitude value, statistic for the noise
                tic
                %xs = simulAG_transtion_loop_G(aes,G,S(j),SC,Cfg.nSub,Cfg.Tmax,omega,Cfg.TRsec,nodo_kick,val);
                xs = hopf_supersub_G_int_pert(a,G,beta,SC,Tmax,w,TR,S(j),nodo_kick,Cfg.nSub);
                if Cfg.filt.bpass==1; xs = filtroign(xs,Cfg.TRsec,Cfg.filt.lb,Cfg.filt.ub); end
                
                toc
                FC_sim_2(:,:,i,j,ii)=corr(xs');
                
            end
            FC_sim_av_2(:,:,j,ii) = mean(FC_sim_2(:,:,:,j,ii),3);
            %figure
            %imagesc(FC_sim_av(:,:,j,ii));title(num2str(ii));
            %S(j) = S;
        end
        %FC_sim_av_nod(:,:,ii) = FC_sim_av(:,:,j);
        nodo_kicked(:,ii) = nodo_kick;
    end
    FC_sim_kick = FC_sim_av_2;
    
    % without paralelaising
    
else
    
    for ii=1:Cfg.nodos_kick; 
        nodo_kick = zeros(115,1);
        if Cfg.nodos_kick==115
            nodo_kick(ii)=1;
        else
            nodo_kick(ii)=1;
            nodo_kick(116-ii)=1; % change depending on the number of nodes
        end
        ii
        for j=1:Cfg.amplitud; 
            S(j)= 0+0.1*(j-1);
            for i=1:Cfg.Repe 
                tic
                %xs = simulAG_transtion_loop_G(aes,G,S(j),SC,Cfg.nSub,Cfg.Tmax,omega,Cfg.TRsec,nodo_kick,val);
                xs = hopf_supersub_G_int_pert(a,G,beta,SC,Tmax,w,TR,S(j),nodo_kick,Cfg.nSub);
                if Cfg.filt.bpass==1; xs = filtroign(xs,Cfg.TRsec,Cfg.filt.lb,Cfg.filt.ub); end
                
                toc
                FC_sim_2(:,:,i,j,ii)=corr(xs');
                
            end
            FC_sim_av_2(:,:,j,ii) = mean(FC_sim_2(:,:,:,j,ii),3);
            %figure
            %imagesc(FC_sim_av(:,:,j,ii));title(num2str(ii));
            %S(j) = S;
        end
        %FC_sim_av_nod(:,:,ii) = FC_sim_av(:,:,j);
        nodo_kicked(:,ii) = nodo_kick;
    end
    FC_sim_kick = FC_sim_av_2;
end