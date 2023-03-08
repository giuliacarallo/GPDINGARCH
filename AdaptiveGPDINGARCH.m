clear all
clc
close all


%%%% Schiphol accidents in 2001 dataset
%save('/Users/giuliacarallo/Documents/GIBBS/Dataset/Schiphol.mat')
load('/Users/giuliacarallo/Desktop/GIBBS/Dataset/Schiphol.mat')
% mu = mean(Z);
% sigma2 = var(Z);

%%
T=length(Z);
n=1;
p=1;
q=1;
%lambda=0;
lambda=0.4;
alpha0 = -0.002;
alpha = 0.30;
%beta = 0;
beta=0.25;
phi = 1/(1-lambda)^2+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set prior Hyperparameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Gamma prior for Phi
a = 5; % anche con a e b = 10 funziona
b = 5;
% Dirichlet prior
sc=0.5;
c = [3*sc*ones(1,p), 4*sc*ones(1,q)];
c0 = 3*sc;
eps=0.90;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Gibbs parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
niter=100000;
nburn=10000;
thin=10;
% Gibbs initialization
alpha0in = alpha0;
alphain = alpha;
betain = beta;
lambdain = lambda;
phiin = phi;

%% Adaptive parameters

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

%% ADAPTIVE SAMPLING
% theta = (alpha0, alpha1, beta1, lambda, phi)

alpha0Gibbs=zeros(1,niter+nburn);
alphaGibbs=zeros(p,niter+nburn);
betaGibbs=zeros(q,niter+nburn);
lambdaGibbs=zeros(1,niter+nburn);
phiGibbs=zeros(1,niter+nburn);
xpost=zeros(T,niter+nburn);
ypost=zeros(T,niter+nburn);
ac=zeros(T,niter+nburn);
% acphi=zeros(nGibbs+nburn,1);
% acvarphi=zeros(nGibbs+nburn,1);
acalpha0=zeros(niter+nburn,1);
% aclambda=zeros(nGibbs+nburn,1);
% acphilambda=zeros(nGibbs+nburn,1);


%%% step 1: simulate varphipost_i from q(varphi_i|-)
i = 1;
xpost(:,i)=ones(T,1);
ypost(:,i)=2*ones(T,1);
phiGibbs(1,i)=phiin;
alphaGibbs(:,i)=alphain;
betaGibbs(:,i)=betain;
alpha0Gibbs(1,i)=alpha0in*0.9;
lambdaGibbs(1,i)= lambdain;
%lkltot(i) = lkl(Z,alpha0Gibbs(i,1),alphaGibbs(i,1),betaGibbs(i,1),lambdaGibbs(i,1),phiGibbs(i,1),xpost(i,:)',ypost(i,:)');

alp_ab = 0.6;
gam_ab =  0.8*((1:niter+nburn).^(-alp_ab));
sig_ab(:,:,i) = diag([0.02; 0.02]);
mu_ab(:,i) = [0.9*alphaGibbs(:,i); 0.9*betaGibbs(:,i)];
lamb_ab(1,i) = 0.001; %0.01;

alp_pl = 0.8;
gam_pl =  0.8*((1:niter+nburn).^(-alp_pl));
sig_pl(:,:,i) = diag([0.2; 0.2]);
mu_pl(:,i) = [0.9*phiGibbs(:,i); 0.9*lambdaGibbs(:,i)];
lamb_pl(1,i) = 0.01; %0.1; %0.01

%%% For the nested models

alp_a = 0.6;
gam_a =  0.8*((1:niter+nburn).^(-alp_ab));
sig_a(:,i) = 0.01;
mu_a(:,i) = 0.9*alphaGibbs(:,i);
lamb_a(1,i) = 0.1; %0.01;

alp_l = 0.8;
gam_l =  0.8*((1:niter+nburn).^(-alp_ab));
sig_l(:,i) = 0.2;
mu_l(:,i) = 0.9*lambdaGibbs(:,i);
lamb_l(1,i) = 0.01; %0.01;



for i = 2:niter+nburn
    % LATENTS X AND Y
    [xpost(:,i),ypost(:,i),ac(:,i)]=xySampl2New(Z,alpha0Gibbs(1,i-1), alphaGibbs(:,i-1), betaGibbs(:,i-1), lambdaGibbs(1,i-1), phiGibbs(1,i-1),xpost(:,i-1),ypost(:,i-1));
  %  xpost(:,i)=X_t;
  %  ypost(:,i)=Y_t;
%     
    lklStar = lkl(Z,alpha0Gibbs(:,i-1), alphaGibbs(:,i-1), betaGibbs(:,i-1), lambdaGibbs(:,i-1),phiGibbs(:,i-1), xpost(:,i),ypost(:,i));
    lklOld = lkl(Z,alpha0Gibbs(:,i-1), alphaGibbs(:,i-1), betaGibbs(:,i-1), lambdaGibbs(:,i-1),phiGibbs(:,i-1), xpost(:,i-1),ypost(:,i-1));
    acpost = lklStar-lklOld;
    u = rand(1,1);
    if log(u) < min(0,acpost)
        xpost(:,i) = xpost(:,i);
        ypost(:,i) = ypost(:,i);
    else
        xpost(:,i) = xpost(:,i-1);
        ypost(:,i) = ypost(:,i-1);
    end
    
    %%%%%%%%%% alpha0
    %  alpha0Gibbs(1,i)=alpha0;
    [alpha0Gibbs(1,i),acalpha0(i)]= alpha0Sampl(Z,alpha0Gibbs(1,i-1), alphaGibbs(:,i-1), betaGibbs(:,i-1), lambdaGibbs(1,i-1), phiGibbs(1,i-1),xpost(:,i),ypost(:,i));
    
    
    % ALPHA AND BETA (INGARCH PARAMETERS)
     [alphaGibbs(:,i),betaGibbs(:,i), rho_ab(1,i),mu_ab(:,i), sig_ab(:,:,i), lamb_ab(1,i)]= varphiAdap(Z,alpha0Gibbs(1,i), alphaGibbs(:,i-1), betaGibbs(:,i-1), lambdaGibbs(1,i-1),phiGibbs(1,i-1),xpost(:,i),ypost(:,i),c,c0,eps,gam_ab(:,i-1),mu_ab(:,i-1), sig_ab(:,:,i-1), lamb_ab(:,i-1));
    % alphaGibbs(:,i)=alpha;
   %  betaGibbs(:,i)=beta;
   % [alphaGibbs(:,i),rho_a(1,i),mu_a(:,i), sig_a(:,i), lamb_a(1,i)]= alphaAdap(Z,alpha0Gibbs(1,i), alphaGibbs(:,i-1), betaGibbs(:,i), lambdaGibbs(1,i-1),phiGibbs(1,i-1),xpost(:,i),ypost(:,i),c,c0,eps,gam_a(:,i-1),mu_a(:,i-1), sig_a(:,i-1), lamb_a(:,i-1));

    
    % LAMBDA AND PHI (GPD PARAMETERS)
    [phiGibbs(1,i), lambdaGibbs(1,i),rho_pl(1,i),mu_pl(:,i), sig_pl(:,:,i), lamb_pl(1,i)] = philambdaAdap(Z,alpha0Gibbs(1,i), alphaGibbs(:,i), betaGibbs(:,i), lambdaGibbs(1,i-1), phiGibbs(1,i-1),xpost(:,i),ypost(:,i),a,b,gam_pl(:,i-1),mu_pl(:,i-1), sig_pl(:,:,i-1), lamb_pl(:,i-1));
   % lambdaGibbs(1,i) = lambda;
   % [phiGibbs(:,i),rho_l(1,i),mu_l(:,i), sig_l(:,i), lamb_l(1,i)]= phiAdap(Z,alpha0Gibbs(1,i), alphaGibbs(:,i), betaGibbs(:,i), lambdaGibbs(1,i),phiGibbs(1,i-1),xpost(:,i),ypost(:,i),a,b,gam_l(:,i-1),mu_l(:,i-1), sig_l(:,i-1), lamb_l(:,i-1));

    
    
    
        message   = sprintf('MCMC iteration: %d/%d\n', i, niter+nburn);
        delString = repmat(sprintf('\b'), 1, length(message));
        fprintf([delString, message])
end

if i > nburn
    xpostNew = xpost(:,nburn+1:thin:end);
    ypostNew = ypost(:,nburn+1:thin:end);
    alpha0GibbsNew = alpha0Gibbs(1,nburn+1:thin:end);
    alphaGibbsNew = alphaGibbs(:,nburn+1:thin:end);
    betaGibbsNew = betaGibbs(:,nburn+1:thin:end);
    phiGibbsNew = phiGibbs(1,nburn+1:thin:end);
    lambdaGibbsNew = lambdaGibbs(1,nburn+1:thin:end);
end
%%
malpha = mean(alphaGibbs);
mbeta = mean(betaGibbs);
mlambda = mean(lambdaGibbs);
mphi = mean(phiGibbs);
SEalpha = std(alphaGibbs)/sqrt(niter+nburn);
SEbeta = std(betaGibbs)/sqrt(niter+nburn);
SElambda = std(lambdaGibbs)/sqrt(niter+nburn);
SEphi = std(phiGibbs)/sqrt(niter+nburn);
% % var = [alphaGibbs,alpha0Gibbs]
% % save('var')


%save('/Users/giuliacarallo/Documents/GIBBS/PermitsmonthlyEst2.mat')
%save('/home/carallo/.matlab/R2019a/REMOTO/CyberAttacksEst.mat')
%save('CyberAttacksEst.mat')


%% Estimation
%
%
% % for i = 1:nGibbs+nburn
% %    if (alphaGibbs(i)<0 || betaGibbs(i)<0) || (alphaGibbs(i)+betaGibbs(i)>=1)
% %       alphaGibbs(i)=[];
% %       betaGibbs(i)=[];
% %       alpha0Gibbs(i)=[];
% %       phiGibbs(i)=[];
% %       lambdaGibbs(i)=[];
% %    end
% % end
%
if i > nburn
    xpostNew = xpost(:,nburn+1:thin:end);
    ypostNew = ypost(:,nburn+1:thin:end);
    alpha0GibbsNew = alpha0Gibbs(1,nburn+1:thin:end);
    alphaGibbsNew = alphaGibbs(:,nburn+1:thin:end);
    betaGibbsNew = betaGibbs(:,nburn+1:thin:end);
    phiGibbsNew = phiGibbs(1,nburn+1:thin:end);
    lambdaGibbsNew = lambdaGibbs(1,nburn+1:thin:end);
end
% % sl=(betaGibbsNew2>0.18);
% % alphaGibbsNew2=alphaGibbsNew2(sl);
% % betaGibbsNew2=betaGibbsNew2(sl);
% % phiGibbsNew2=phiGibbsNew2(sl);
% % lambdaGibbsNew2=lambdaGibbsNew2(sl);
% %
% % for i = 1:length(betaGibbsNew2)
% %    if betaGibbsNew2(i) < 0.18
% %       alphaGibbsNew2(i)=[];
% %       betaGibbsNew2(i)=[];
% %       phiGibbsNew2(i)=[];
% %       lambdaGibbsNew2(i)=[];
% %    end
% % end
% 
 %%%%%%%%%% mean of the parameters
 malpha0N = mean(alpha0GibbsNew);
malphaN = mean(alphaGibbsNew);
mbetaN = mean(betaGibbsNew);
mlambdaN = mean(lambdaGibbsNew);
mphiN = mean(phiGibbsNew);

%%%%%%%%%% standard deviation of the parameters
SDalphaN = std(alphaGibbsNew);
SDalpha0N = std(alpha0GibbsNew);
SDbetaN = std(betaGibbsNew);
SDlambdaN = std(lambdaGibbsNew);
SDphiN = std(phiGibbsNew);
% 
% %%%%%%%%%% credible intervals for the parameters
quantlambdaN = quantile(lambdaGibbsNew, [0.025 0.975]);
quantalphaN = quantile(alphaGibbsNew, [0.025 0.975]);
quantbetaN = quantile(betaGibbsNew, [0.025 0.975]);
quantphiN = quantile(phiGibbsNew, [0.025 0.975]);
quantalpha0N = quantile(alpha0GibbsNew, [0.025 0.975]);

% tsa = tinv([0.025  0.975],length(alphaGibbsNew)-1);
% tsb = tinv([0.025  0.975],length(betaGibbsNew)-1);
% tsl = tinv([0.025  0.975],length(lambdaGibbsNew)-1);
% tsp = tinv([0.025  0.975],length(phiGibbsNew)-1);
% %
% %
% %
% %
% % CIalpha=[malphaN-1.96*SDalphaN, malphaN+1.96*SDalphaN];
% % CIbeta=[mbetaN-1.96*SDbetaN, mbetaN+1.96*SDbetaN];
% % CIlambda=[mlambdaN-1.96*SDlambdaN, mlambdaN+1.96*SDlambdaN];
% % CIphi=[mphiN-1.96*SDphiN, mphiN+1.96*SDphiN];
% %
% % % %
%save('/Users/giuliacarallo/Desktop/GIBBS/SchipholAdaptive.mat')
%save('/Users/giuliacarallo/Desktop/GIBBS/Code_Cyber/SchipholAdaptiveGPIN.mat')
%save('/Users/giuliacarallo/Desktop/GIBBS/Code_Cyber/SchipholAdaptivePDIN.mat')
%save('/Users/giuliacarallo/Desktop/GIBBS/Code_Cyber/SchipholAdaptivePD.mat')
