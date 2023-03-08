function [mut,mutstar, sigma2t,sigma2tstar, thetaxt, thetayt] = mut(Z,alpha0, alpha, beta, lambda, phi)
T = size(Z,1);
p = size(alpha,1);
q = size(beta,1);
mut=zeros(p+q+T,1);
Z0=zeros(p+q,1);
Ztilde=[Z0;Z];
sigma2t=zeros(p+q+T,1);
sigma2tstar=zeros(p+q+T,1);  % AGGIUNTO
thetaxt=zeros(p+q+T,1);
thetayt=zeros(p+q+T,1);
mutstar=zeros(p+q+T,1);

for t=p+q+1:p+q+T
   mut(t,1)=alpha0+alpha*Ztilde(t-(1:p))+beta*mut(t-(1:q));
   mutstar(t,1)=mut(t,1)*(1-lambda);
   sigma2t(t,1)=phi*abs(mut(t,1));
   sigma2tstar(t,1) = sigma2t(t,1)*(1-lambda)^3;   % AGGIUNTO 
   thetaxt(t,1)= abs((sigma2tstar(t,1)+mutstar(t,1))/2);
   thetayt(t,1)= abs((sigma2tstar(t,1)-mutstar(t,1))/2);
end

mut = mut(p+q+1:end,1);
mutstar = mutstar(p+q+1:end,1);
sigma2t = sigma2t(p+q+1:end,1);
sigma2tstar = sigma2tstar(p+q+1:end,1);  % AGGIUNTO
thetaxt = thetaxt(p+q+1:end,1);
thetayt = thetayt(p+q+1:end,1);

end

