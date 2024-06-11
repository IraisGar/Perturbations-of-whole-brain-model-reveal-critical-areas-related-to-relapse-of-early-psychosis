function [xs]=hopf_supersub_G_int_pert(a,G,beta,C,Tmax,w, TR,S,kick, nSubs)

% integrates the model's equations
sig = 0.01; 
dt=0.1*TR/2;
dsig = sqrt(dt)*sig; 

nNodes=length(a);
N = nNodes;

a=repmat(a,1,2);
wC = C ; 



wC = G.*C;
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
S_in = repmat(S,1,2);

kick = repmat(kick,1,2);

omega1 = repmat(2*pi*w',1,2); omega1(:,1) = -omega1(:,1);
omega2=beta.*ones(nNodes,2);
omega2(:,1) = -omega2(:,1);
omega=omega1+omega2;

omega_pos = omega;
omega_pos(:,1) = omega(:,2);

%xs=zeros(Tmax,N);
xs=zeros(nNodes,Tmax*nSubs);
z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
nn=0;
% discard first 2000 time steps
for t=0:dt:2000
    zz = z(:,end:-1:1);
    z = z + dt*(a.*z + zz.*omega - (z+omega2.*zz).*(z.*z+zz.*zz) + wC*(z) - sumC.*z + kick.*S_in.*cos(omega_pos.*t)) +dsig*randn(N,2);%G.* (wC*z - sumC.*z)) + dsig*randn(N,2);
end
for t=0:dt:((Tmax-1)*TR)
    zz = z(:,end:-1:1);
    z = z + dt*(a.*z + zz.*omega - (z+omega2.*zz).*(z.*z+zz.*zz) + wC*(z) - sumC.*z + kick.*S_in.*cos(omega_pos.*t)) +dsig*randn(N,2);
    if abs(mod(t,TR))<0.01
        nn=nn+1;
        xs(:,nn)=z(:,1)';
    end
end
