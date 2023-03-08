
clear all
close all
clc

% % Cyber Attacks dataset - 2017 e 2018 daily + difference
% poolobj = parpool(30);
% addAttachedFiles(poolobj);
% updateAttachedFiles(poolobj);
%addpath('/home/giulia.carallo/.matlab/R2022b');
%5load('/home/giulia.carallo/.matlab/R2022b/CyberAttacks2017_2018.mat','Z_old','TOTAL');
%load('/home/giulia.carallo/Schiphol.mat')
%load('CyberAttacks2017_2018.mat','Z_old','TOTAL')
%load('/Users/giuliacarallo/Desktop/GIBBS/Dataset/Schiphol.mat')


%load('/home/giulia.carallo/.matlab/R2022b/Schiphol.mat','Z','OBSV');
%%
Z_old=Z; % for schiphol
TOTAL = OBSV; % for schiphol
TT=size(Z_old,1);
%wn=365; % window size
wn = 160; % window size for schiphol
niter=110000;
niter2=5000;
nburn=500;
thin = 20;
th = TT-wn+1;  % number rolling windows
H = 7;
V_old =  TOTAL;
%% Forecast variables
mut_f=zeros(niter,wn+H);
mutstar_f=zeros(niter,wn+H);
sigma2t_f=zeros(niter,wn+H);
sigma2tstar_f=zeros(niter,wn+H);
mut_fSeq=zeros(niter,wn+H);
mutstar_fSeq=zeros(niter,wn+H);
sigma2t_fSeq=zeros(niter,wn+H);
sigma2tstar_fSeq=zeros(niter,wn+H);
ZfSeq=zeros(niter,wn+H);
Vf = zeros(niter,H);

muCum=zeros(th,wn+H);
mustarCum=zeros(th,wn+H);
sigma2Cum=zeros(th,1);
sigma2starCum=zeros(th,1);

alpha0Cum = zeros(th,1);
alphaCum = zeros(th,1);
betaCum = zeros(th,1);
lambdaCum = zeros(th,1);
phiCum = zeros(th,1);
ZfCum = zeros(th,wn+H);
VCum = zeros(th,H);
xpostCum = zeros(th,wn);
ypostCum = zeros(th,wn);

muCumQ1=zeros(th,wn+H);
mustarCumQ1=zeros(th,wn+H);
ZfCumQ1=zeros(th,wn+H);
VfCumQ1=zeros(th,H);
muCumQ2=zeros(th,wn+H);
mustarCumQ2=zeros(th,wn+H);
ZfCumQ2=zeros(th,wn+H);
VfCumQ2=zeros(th,H);

%% First window
Z = Z_old(1:wn);
V = V_old(1:wn);

t=1;
n=1;
p=1;
q=1;

Z0=zeros(p+q,1);

lambda=0.4;
alpha = 0.40;
%beta = 0.25;
beta=0;
alpha0 = mean(Z)*(1-alpha-beta);
phi = 1/(1-lambda)^2+60;

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

% Gibbs initialization
alpha0in = alpha0;
alphain = alpha;
betain = beta;
lambdain = lambda;
%phiin = 1/((1-lambda)^2)+2;
phiin = phi;
xpostin=ones(wn,1);
ypostin=2*ones(wn,1);


%%
% theta = (alpha0, alpha1, beta1, lambda, phi)

alpha0Gibbs=zeros(1,niter+nburn);
alphaGibbs=zeros(p,niter+nburn);
betaGibbs=zeros(q,niter+nburn);
lambdaGibbs=zeros(1,niter+nburn);
phiGibbs=zeros(1,niter+nburn);
xpost=zeros(wn,niter+nburn);
ypost=zeros(wn,niter+nburn);
ac=zeros(wn,niter+nburn);
% acphi=zeros(niter+nburn,1);
% acvarphi=zeros(niter+nburn,1);
acalpha0=zeros(niter+nburn,1);
% aclambda=zeros(niter+nburn,1);
% acphilambda=zeros(niter+nburn,1);

%%% step 1: simulate varphipost_i from q(varphi_i|-)
i = 1;
xpost(:,i)=xpostin;
ypost(:,i)=ypostin;
phiGibbs(1,i)=phiin;
alphaGibbs(:,i)=alphain;
betaGibbs(:,i)=betain;
alpha0Gibbs(1,i)=alpha0in*0.9;
lambdaGibbs(1,i)= lambdain;

[xpost, ypost, alpha0Gibbs, alphaGibbs, betaGibbs, phiGibbs, lambdaGibbs, ZfSeq, Vf, mut_f, mutstar_f] = gibbsadap(niter, nburn, wn, H, p, q, alpha0in, alphain, betain,phiin,lambdain,xpostin,ypostin,Z,Z0,V);


alpha0Cum(t,1) = mean(alpha0Gibbs(1,nburn+1:thin:end));
alphaCum(t,1) = mean(alphaGibbs(:,nburn+1:thin:end));
betaCum(t,1) = mean(betaGibbs(:,nburn+1:thin:end));
phiCum(t,1) = mean(phiGibbs(1,nburn+1:thin:end));
lambdaCum(t,1) = mean(lambdaGibbs(1,nburn+1:thin:end));
ZfCum(t,:) = mean(ZfSeq(nburn+1:thin:end,:));
VCum(t,:) = mean(Vf(nburn+1:thin:end,:));
muCum(t,:)=mean(mut_f(nburn+1:thin:end,:));
mustarCum(t,:)=mean(mutstar_f(nburn+1:thin:end,:));

muCumQ1(t,:)=quantile(mut_f(nburn+1:thin:end,:),0.025);
muCumQ2(t,:)=quantile(mut_f(nburn+1:thin:end,:),0.975);
mustarCumQ1(t,:)=quantile(mutstar_f(nburn+1:thin:end,:),0.025);
mustarCumQ2(t,:)=quantile(mutstar_f(nburn+1:thin:end,:),0.975);
ZfCumQ1(t,:)=quantile(ZfSeq(nburn+1:thin:end,:),0.025);
ZfCumQ2(t,:)=quantile(ZfSeq(nburn+1:thin:end,:),0.975);
VfCumQ1(t,:)=quantile(Vf(nburn+1:thin:end,:),0.025);
VfCumQ2(t,:)=quantile(Vf(nburn+1:thin:end,:),0.975);

xpostCum(t,:) = floor(mean(xpost(:,nburn+1:thin:end),2));
ypostCum(t,:) = floor(mean(ypost(:,nburn+1:thin:end),2));



%% Other windows
alpha0in = alpha0Cum(t,1);
alphain = alphaCum(t,1);
betain = betaCum(t,1);
lambdain = lambdaCum(t,1);
phiin = phiCum(t,1);
xpostin = xpostCum(t,:);
ypostin = ypostCum(t,:);

parfor t=2:th
    n=1;
    p=1;
    q=1;
    Z0=zeros(p+q,1);
    Z = Z_old(t:t-1+wn);
    V = V_old(t:t-1+wn);
    
    
    
    [xpost2, ypost2, alpha0Gibbs2, alphaGibbs2, betaGibbs2, phiGibbs2, lambdaGibbs2, ZfSeq2, Vf2, mut_f2, mutstar_f2] = gibbsadap(niter2, nburn, wn, H, p, q, alpha0in, alphain, betain,phiin,lambdain,xpostin,ypostin,Z,Z0,V);
    
    
 
    
    xpostCum(t,:) = floor(mean(xpost2(:,nburn+1:thin:end),2));
    ypostCum(t,:) = floor(mean(ypost2(:,nburn+1:thin:end),2));
    alpha0Cum(t,1) = mean(alpha0Gibbs2(1,nburn+1:thin:end));
    alphaCum(t,1) = mean(alphaGibbs2(:,nburn+1:thin:end));
    betaCum(t,1) = mean(betaGibbs2(:,nburn+1:thin:end));
    phiCum(t,1) = mean(phiGibbs2(1,nburn+1:thin:end));
    lambdaCum(t,1) = mean(lambdaGibbs2(1,nburn+1:thin:end));
    ZfCum(t,:) = mean(ZfSeq2(nburn+1:thin:end,:));
    VCum(t,:) = mean(Vf2(nburn+1:thin:end,:));
    muCum(t,:)=mean(mut_f2(nburn+1:thin:end,:));
    mustarCum(t,:)=mean(mutstar_f2(nburn+1:thin:end,:));
    
    muCumQ1(t,:)=quantile(mut_f2(nburn+1:thin:end,:),0.025);
    muCumQ2(t,:)=quantile(mut_f2(nburn+1:thin:end,:),0.975);
    mustarCumQ1(t,:)=quantile(mutstar_f2(nburn+1:thin:end,:),0.025);
    mustarCumQ2(t,:)=quantile(mutstar_f2(nburn+1:thin:end,:),0.975);
    ZfCumQ1(t,:)=quantile(ZfSeq2(nburn+1:thin:end,:),0.025);
    ZfCumQ2(t,:)=quantile(ZfSeq2(nburn+1:thin:end,:),0.975);
    VfCumQ1(t,:)=quantile(Vf2(nburn+1:thin:end,:),0.025);
    VfCumQ2(t,:)=quantile(Vf2(nburn+1:thin:end,:),0.975);
    
    
    message2   = sprintf('Forecast Window: %d/%d\n', t, th);
    delString2 = repmat(sprintf('\b'), 1, length(message2));
    fprintf([delString2, message2])
end

figure(1)
plot(1:wn,Z_old(1:wn));
hold on ;
plot(wn+1:TT-1,Z_old(wn+1:end-1),'k-');
plot((wn+1):TT-1,ZfCum(1:(end-2),wn+1),'r--');

V_old_1 = V_old(2:end);
figure(2)
plot(1:wn,V_old_1(1:wn));
hold on ;
plot(wn+1:TT-1,V_old_1(wn+1:end-1),'k-');
plot((wn+1):TT-1,VCum(1:(end-2),1),'r--');


%save('/Users/giuliacarallo/Dropbox/GPDingarchFOREC/CyberSeqAdapForecGPD.mat')
%save('/Users/giuliacarallo/Dropbox/GPDingarchFOREC/SchipSeqAdapForecGPD.mat')

%save('/home/giulia.carallo/.matlab/R2022b/CyberSeqAdapForecGPD.mat')
%save('/home/giulia.carallo/SchipSeqAdapForecGPD.mat')

%save('/home/giulia.carallo/.matlab/R2022b/SchipSeqAdapForecGPDIN.mat')
