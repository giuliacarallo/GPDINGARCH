function [mut_f,mutstar_f, sigma2t_f,sigma2tstar_f,mut_fDyn,mutstar_fDyn, sigma2t_fDyn,sigma2tstar_fDyn, ZtildeDyn, Vf] = mut_forec2(Z0,Z,V,alpha0, alpha, beta, lambda, phi, H)
%[mut_f,mutstar_f, sigma2t_f,sigma2tstar_f,mut_fDyn,mutstar_fDyn, sigma2t_fDyn,sigma2tstar_fDyn, ZtildeDyn, Vf, VDyn]
T = size(Z,1);
p = size(alpha,1);
q = size(beta,1);
mut_f=zeros(p+q+T+H,1);
Zf =zeros(H,1);
Vf=zeros(H,1);
VDyn=zeros(H,1);

%Ztilde=[Z0;Z;Zf];
ZtildeDyn=[Z0;Z;Zf];
%ZtildeAll=[Z0;Z;Z];
V0 = V(T,1);
sigma2t_f=zeros(p+q+T+H,1);
sigma2tstar_f=zeros(p+q+T+H,1);  % AGGIUNTO
mutstar_f=zeros(p+q+T+H,1);
sigma2t_fDyn=zeros(p+q+T+H,1);
sigma2tstar_fDyn=zeros(p+q+T+H,1);  % AGGIUNTO
mutstar_fDyn=zeros(p+q+T+H,1);
mut_fDyn=zeros(p+q+T+H,1);

for t=p+q+1:p+q+T
   mut_f(t,1)=alpha0+alpha*ZtildeDyn(t-(1:p))+beta*mut_f(t-(1:q));
   mutstar_f(t,1)=mut_f(t,1)*(1-lambda);
   sigma2t_f(t,1)=phi*abs(mut_f(t,1));
   sigma2tstar_f(t,1) =sigma2t_f(t,1)*(1-lambda)^3;   % AGGIUNTO
   mut_fDyn(t,1)=mut_f(t,1);
   mutstar_fDyn(t,1)=mutstar_f(t,1);
   sigma2t_fDyn(t,1)=sigma2t_f(t,1);
   sigma2tstar_fDyn(t,1)=sigma2tstar_f(t,1);
   %    thetaxt(t,1)= (sigma2tstar(t,1)+mutstar(t,1))/2;
   %    thetayt(t,1)= (sigma2tstar(t,1)-mutstar(t,1))/2;
end
mut_f = mut_f(p+q+1:end,1);
mutstar_f = mutstar_f(p+q+1:end,1);
sigma2t_f = sigma2t_f(p+q+1:end,1);
sigma2tstar_f = sigma2tstar_f(p+q+1:end,1);
%Ztilde=Ztilde(p+q+1:end,1);

mut_fDyn = mut_fDyn(p+q+1:end,1);
mutstar_fDyn = mutstar_fDyn(p+q+1:end,1);
sigma2t_fDyn = sigma2t_fDyn(p+q+1:end,1);
sigma2tstar_fDyn = sigma2tstar_fDyn(p+q+1:end,1);
ZtildeDyn=ZtildeDyn(p+q+1:end,1);
%ZtildeAll=ZtildeAll(p+q+1:end,1);

%primo step Forecast
% h = 1;
% mut_f(T+h,1)=alpha0+alpha*Ztilde(T+h-(1:p))+beta*mut_f(T+h-(1:q));
% mutstar_f(T+h,1)=mut_f(T+h,1)*(1-lambda);
% sigma2t_f(T+h,1)=phi*abs(mut_f(T+h,1));
% sigma2tstar_f(T+h,1) =sigma2t_f(T+h,1)*(1-lambda)^3;   % AGGIUNTO
% thetaxt= (sigma2tstar_f(T+h,1)+mutstar_f(T+h,1))/2;
% thetayt= (sigma2tstar_f(T+h,1)-mutstar_f(T+h,1))/2;
% X_f = gpbranching(thetaxt,lambda,1)';
% Y_f = gpbranching(thetayt,lambda,1)';
% Ztilde(T+h,1)= X_f-Y_f;


for h=1:H
  
   mut_f(T+h,1)=alpha0+alpha*ZtildeDyn(T+h-(1:p))+beta*mut_f(T+h-(1:q));
   mutstar_f(T+h,1)=mut_f(T+h,1)*(1-lambda);
   sigma2t_f(T+h,1)=phi*abs(mut_f(T+h,1));
   sigma2tstar_f(T+h,1) =sigma2t_f(T+h,1)*(1-lambda)^3;   % AGGIUNTO
   thetaxt= (sigma2tstar_f(T+h,1)+mutstar_f(T+h,1))/2;
   thetayt= (sigma2tstar_f(T+h,1)-mutstar_f(T+h,1))/2;
   cnd = 1;
   while cnd
   X_f = gpbranching(thetaxt,lambda,1);
   Y_f = gpbranching(thetayt,lambda,1);
   ZtildeDyn(T+h,1)= X_f-Y_f;
   
   if h == 1
     Vf(h,1) = V0 + ZtildeDyn(T+h,1);
   else
   Vf(h,1) = V(h-1,1) + ZtildeDyn(T+h,1);
   end
   cnd = Vf(h,1) < 0;
   end
   %%
%    mut_fDyn(T+h,1)=alpha0+alpha*ZtildeAll(T+h-(1:p))+beta*mut_fDyn(T+h-(1:q));
%    mutstar_fDyn(T+h,1)=mut_fDyn(T+h,1)*(1-lambda);
%    sigma2t_fDyn(T+h,1)=phi*abs(mut_fDyn(T+h,1));
%    sigma2tstar_fDyn(T+h,1) =sigma2t_fDyn(T+h,1)*(1-lambda)^3;   % AGGIUNTO
%    thetaxt= (sigma2tstar_fDyn(T+h,1)+mutstar_fDyn(T+h,1))/2;
%    thetayt= (sigma2tstar_fDyn(T+h,1)-mutstar_fDyn(T+h,1))/2;
%    X_f = gpbranching(thetaxt,lambda,1);
%    Y_f = gpbranching(thetayt,lambda,1);
%    ZtildeDyn(T+h,1)= X_f-Y_f;  
%    if h == 1
%       VDyn(h,1) = V0 + ZtildeDyn(T+h,1);
%    else
%    VDyn(h,1) = V(h-1,1) + ZtildeDyn(T+h,1);
%    end
end
mut_f = mut_f';
mutstar_f = mutstar_f';
sigma2t_f = sigma2t_f';
sigma2tstar_f = sigma2tstar_f';
%Ztilde=Ztilde';
%
% mut_fDyn = mut_fDyn';
% mutstar_fDyn = mutstar_fDyn';
% sigma2t_fDyn = sigma2t_fDyn';
% sigma2tstar_fDyn = sigma2tstar_fDyn';
% ZtildeDyn=ZtildeDyn';
Vf = Vf';
% VDyn = VDyn';

% mut_f = mut_f(p+q+1:end,1);
% mutstar_f = mutstar_f(p+q+1:end,1);
% sigma2t_f = sigma2t_f(p+q+1:end,1);
% sigma2tstar_f = sigma2tstar_f(p+q+1:end,1);  % AGGIUNTO
% % thetaxt_f = thetaxt_f(p+q+1:end,1);
% % thetayt_f = thetayt_f(p+q+1:end,1);
% end