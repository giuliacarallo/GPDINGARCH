function [xpost2, ypost2, alpha0Gibbs2, alphaGibbs2, betaGibbs2, phiGibbs2, lambdaGibbs2, ZfSeq2, Vf2, mut_f2, mutstar_f2] = gibbsadap(niter,nburn, wn, H, p, q, alpha0in, alphain, betain,phiin,lambdain,xpostin,ypostin,Z,Z0,V)

mut_f2=zeros(niter,wn+H);
mutstar_f2=zeros(niter,wn+H);
sigma2t_f2=zeros(niter,wn+H);
sigma2tstar_f2=zeros(niter,wn+H);
mut_fSeq2=zeros(niter,wn+H);
mutstar_fSeq2=zeros(niter,wn+H);
sigma2t_fSeq2=zeros(niter,wn+H);
sigma2tstar_fSeq2=zeros(niter,wn+H);
ZfSeq2=zeros(niter,wn+H);
Vf2 = zeros(niter,H);

rho_pl = zeros(1, niter+nburn);
sig_pl = zeros(2,2,niter+nburn);
mu_pl = zeros(2,niter+nburn);
lamb_pl = zeros(1,niter+nburn);

rho_ab = zeros(1, niter+nburn);
sig_ab = zeros(2,2,niter+nburn);
mu_ab = zeros(2,niter+nburn);
lamb_ab = zeros(1,niter+nburn);

rho_a = zeros(1, niter+nburn);
sig_a = zeros(1,niter+nburn);
mu_a = zeros(1,niter+nburn);
lamb_a = zeros(1,niter+nburn);

rho_l = zeros(1, niter+nburn);
sig_l = zeros(1,niter+nburn);
mu_l = zeros(1,niter+nburn);
lamb_l = zeros(1,niter+nburn);


% Gamma prior for Phi
a = 5; % anche con a e b = 10 funziona
b = 5;
% Dirichlet prior
sc=0.5;
c = [3*sc*ones(1,p), 4*sc*ones(1,q)];
c0 = 3*sc;
eps=0.90;

% Gibbs initialization


alpha0Gibbs2=zeros(1,niter);
alphaGibbs2=zeros(p,niter);
betaGibbs2=zeros(q,niter);
lambdaGibbs2=zeros(1,niter);
phiGibbs2=zeros(1,niter);
xpost2=zeros(wn,niter);
ypost2=zeros(wn,niter);
acalpha02=zeros(1,niter);

%%% step 1: simulate varphipost_i from q(varphi_i|-)
i=1;
xpost2(:,i)=xpostin;
ypost2(:,i)=ypostin;
phiGibbs2(1,i)=phiin;
alphaGibbs2(:,i)=alphain;
betaGibbs2(:,i)=betain;
alpha0Gibbs2(1,i)=alpha0in;
lambdaGibbs2(1,i)= lambdain;

alp_ab = 0.6;
gam_ab =  0.8*((1:niter+nburn).^(-alp_ab));
sig_ab(:,:,i) = diag([0.2; 0.2]);
mu_ab(:,i) = [0.9*alphaGibbs2(:,i); 0.9*betaGibbs2(:,i)];
lamb_ab(1,i) = 0.001; %0.01;

alp_pl = 0.8;
gam_pl =  0.8*((1:niter+nburn).^(-alp_pl));
sig_pl(:,:,i) = diag([0.2; 0.2]);
mu_pl(:,i) = [0.9*phiGibbs2(:,i); 0.9*lambdaGibbs2(:,i)];
lamb_pl(1,i) = 0.01; %0.1; %0.01

%%% For the nested models

alp_a = 0.6;
gam_a =  0.8*((1:niter+nburn).^(-alp_a));
sig_a(:,i) = 0.01;
mu_a(:,i) = 0.9*alphaGibbs2(:,i);
lamb_a(1,i) = 0.1; %0.01;

alp_l = 0.8;
gam_l =  0.8*((1:niter+nburn).^(-alp_ab));
sig_l(:,i) = 0.2;
mu_l(:,i) = 0.9*lambdaGibbs2(:,i);
lamb_l(1,i) = 0.01; %0.01;





for i = 2:niter
    
    % LATENTS X AND Y
    [xpost2(:,i),ypost2(:,i),ac2(:,i)]=xySampl2New(Z,alpha0Gibbs2(1,i-1), alphaGibbs2(:,i-1), betaGibbs2(:,i-1), lambdaGibbs2(1,i-1), phiGibbs2(1,i-1),xpost2(:,i-1),ypost2(:,i-1));
    %  xpost(:,i)=X_t;
    %  ypost(:,i)=Y_t;
    %
    lklStar2 = lkl(Z,alpha0Gibbs2(:,i-1), alphaGibbs2(:,i-1), betaGibbs2(:,i-1), lambdaGibbs2(:,i-1),phiGibbs2(:,i-1), xpost2(:,i),ypost2(:,i));
    lklOld2 = lkl(Z,alpha0Gibbs2(:,i-1), alphaGibbs2(:,i-1), betaGibbs2(:,i-1), lambdaGibbs2(:,i-1),phiGibbs2(:,i-1), xpost2(:,i-1),ypost2(:,i-1));
    acpost = lklStar2-lklOld2;
    u = rand(1,1);
    if log(u) < min(0,acpost)
        xpost2(:,i) = xpost2(:,i);
        ypost2(:,i) = ypost2(:,i);
    else
        xpost2(:,i) = xpost2(:,i-1);
        ypost2(:,i) = ypost2(:,i-1);
    end
    
    %%%%%%%%%% alpha0
    %  alpha0Gibbs(1,i)=alpha0;
    [alpha0Gibbs2(1,i),acalpha02(i)]= alpha0Sampl(Z,alpha0Gibbs2(1,i-1), alphaGibbs2(:,i-1), betaGibbs2(:,i-1), lambdaGibbs2(1,i-1), phiGibbs2(1,i-1),xpost2(:,i),ypost2(:,i));
    
    
    % ALPHA AND BETA (INGARCH PARAMETERS)
   % [alphaGibbs2(:,i),betaGibbs2(:,i), rho_ab(1,i),mu_ab(:,i), sig_ab(:,:,i), lamb_ab(1,i)]= varphiAdap(Z,alpha0Gibbs2(1,i), alphaGibbs2(:,i-1), betaGibbs2(:,i-1), lambdaGibbs2(1,i-1),phiGibbs2(1,i-1),xpost2(:,i),ypost2(:,i),c,c0,eps,gam_ab(:,i-1),mu_ab(:,i-1), sig_ab(:,:,i-1), lamb_ab(:,i-1));
    % alphaGibbs(:,i)=alpha;
      betaGibbs2(:,i)=betain;
     [alphaGibbs2(:,i),rho_a(1,i),mu_a(:,i), sig_a(:,i), lamb_a(1,i)]= alphaAdap(Z,alpha0Gibbs2(1,i), alphaGibbs2(:,i-1), betaGibbs2(:,i), lambdaGibbs2(1,i-1),phiGibbs2(1,i-1),xpost2(:,i),ypost2(:,i),c,c0,eps,gam_a(:,i-1),mu_a(:,i-1), sig_a(:,i-1), lamb_a(:,i-1));
    
    % LAMBDA AND PHI (GPD PARAMETERS)
   % [phiGibbs2(1,i), lambdaGibbs2(1,i),rho_pl(1,i),mu_pl(:,i), sig_pl(:,:,i), lamb_pl(1,i)] = philambdaAdap(Z,alpha0Gibbs2(1,i), alphaGibbs2(:,i), betaGibbs2(:,i), lambdaGibbs2(1,i-1), phiGibbs2(1,i-1),xpost2(:,i),ypost2(:,i),a,b,gam_pl(:,i-1),mu_pl(:,i-1), sig_pl(:,:,i-1), lamb_pl(:,i-1));
      lambdaGibbs2(1,i) = lambdain;
      [phiGibbs2(:,i),rho_l(1,i),mu_l(:,i), sig_l(:,i), lamb_l(1,i)]= phiAdap(Z,alpha0Gibbs2(1,i), alphaGibbs2(:,i), betaGibbs2(:,i), lambdaGibbs2(1,i),phiGibbs2(1,i-1),xpost2(:,i),ypost2(:,i),a,b,gam_l(:,i-1),mu_l(:,i-1), sig_l(:,i-1), lamb_l(:,i-1));
    
    
    %%%%% Forecasting
    %    [mut_f(i,:),mutstar_f(i,:), sigma2t_f(i,:),sigma2tstar_f(i,:), Zf(i,:),...
    %       mut_fDyn(i,:),mutstar_fDyn(i,:), sigma2t_fDyn(i,:),sigma2tstar_fDyn(i,:), ZfDyn(i,:)] = mut_forec(Z,Zout,alpha0Gibbs(i,1), alphaGibbs(i,:), betaGibbs(i,:), lambdaGibbs(i,1), phiGibbs(i,1), H);
    %
    
    [mut_f2(i,:),mutstar_f2(i,:), sigma2t_f2(i,:),sigma2tstar_f2(i,:),...
        mut_fSeq2(i,:),mutstar_fSeq2(i,:), sigma2t_fSeq2(i,:),sigma2tstar_fSeq2(i,:), ZfSeq2(i,:), Vf2(i,:)] = mut_forec2(Z0,Z,V,alpha0Gibbs2(:,i),alphaGibbs2(:,i), betaGibbs2(:,i), lambdaGibbs2(1,i), phiGibbs2(1,i),H);
    %% 
    %% 
    
    
    %
    %     message   = sprintf('MCMC iteration: %d/%d\n', i, niter);
    %     delString = repmat(sprintf('\b'), 1, length(message));
    %     fprintf([delString, message])
end