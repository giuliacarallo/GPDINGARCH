function [phiGibbsnew, lambdaGibbsnew,rho_pl,mu_plnew, sig_plnew, lamb_plnew] = philambdaAdap(Z,alpha0Gibbs, alphaGibbs, betaGibbs, lambdaGibbs, phiGibbs,xpost,ypost,a,b,gam_pl,mu_plold,sig_plold,lamb_plold)

philamGibbs = [phiGibbs;  lambdaGibbs];

cnd = true;
while cnd 
philamStar = philamGibbs + sqrt(lamb_plold)*chol(sig_plold,'lower')*randn(2,1);

phiStar = philamStar(1,:);
lambdaStar = philamStar(2,:);

k = 1/(1-lambdaStar)^2;
cnd = phiStar < k || (0 > lambdaStar) || (lambdaStar > 1);
end

s=100;
qstar=log(betapdf(lambdaStar,lambdaGibbs*s,(1-lambdaGibbs)*s));
qGibbs=log(betapdf(lambdaGibbs, lambdaStar*s,(1-lambdaStar)*s));

c = 1/(1-lambdaGibbs)^2;
prstar=log(pdf('gamma',phiStar-k,a,b));
prGibbs=log(pdf('gamma',phiGibbs-c,a,b));

lklstar = lkl(Z,alpha0Gibbs, alphaGibbs, betaGibbs, lambdaStar,phiStar, xpost,ypost);
lklold = lkl(Z,alpha0Gibbs, alphaGibbs, betaGibbs, lambdaGibbs,phiGibbs, xpost,ypost);

ac = lklstar-lklold+prstar-prGibbs+qGibbs-qstar;
rho_pl  = exp(min(0,ac));


u = rand(1,1);
if log(u) < min(0,ac) 
    philamGibbsnew = philamStar;
else
    philamGibbsnew = philamGibbs;
    rho_pl = 0;
end

phiGibbsnew = philamGibbsnew(1,:);
lambdaGibbsnew = philamGibbsnew(2,:);



lamb_plnew = exp(log(lamb_plold) + gam_pl * (rho_pl - 0.34));
mu_plnew = mu_plold + gam_pl*(philamGibbsnew - mu_plold);
sig_plnew = sig_plold + gam_pl*((philamGibbsnew-mu_plold)*(philamGibbsnew-mu_plold)'-sig_plold);

end
