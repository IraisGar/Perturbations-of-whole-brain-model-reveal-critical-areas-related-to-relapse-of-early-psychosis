
function [xs]=simulAG_transtion_loop_G(a,G,S,C,nSubs,Tmax,omega,TRsec,kick,val)

    
    % a: (n,1)
    % G: (n,1)--> it is computed internaly
    % S: amplitude of the kick
    % w: (n,2) --> comes from "saca_w" where each column has the same but 
    % with opposite signs (to give ODES a complex form)
    % kick: (n,1) --> 1 in the nodes we want to kick.
    % n : size of the used parcellation, in this example n = 90
    
    dt = 0.1;
    sig = 0.04; 
    dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step
    %fix G = 0.5
    %G = ones(90,1)*0.5;
    nNodes=length(a);
    
    a=repmat(a,1,2);
    
    G=repmat(G,1,2);
    S_in = repmat(S,1,2);
    
    kick = repmat(kick,1,2);
    
    omega_pos = omega;
    omega_pos(:,1) = omega(:,2);
        
    wC = C; 

    xs=zeros(nNodes,Tmax*nSubs);
    z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)

    sumC = repmat(sum(wC,2),1,2);

   % [start simulating with "a" (remember that neither "a" nor "G" are
   % homogeneous from this point on
     
          nn=1; % ends the thing
        for t=1:dt:1000 
            %suma =  wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + G.* (wC*z - sumC.*z)+kick.*S_in.*cos(omega_pos.*t)) + dsig*randn(nNodes,2);
        end
 
        t=1:dt:Tmax*TRsec*nSubs;
        
        
          for i=1:length(t) %JVS: was 15000, now faster
            %suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + G.* (wC*z - sumC.*z)+kick.*S_in.*cos(omega_pos.*t(i))) + dsig*randn(nNodes,2);
% 
       
       

            if t(i)==val(nn)
                
                %REMEMBER time series have to have a row based shape
                xs(:,nn)=z(:,1);
                nn=nn+1;
                
            end
        
        
        end
        
