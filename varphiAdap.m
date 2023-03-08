function [alphaGibbsnew, betaGibbsnew,rho_ab,mu_abnew, sig_abnew, lamb_abnew] = varphiAdap(Z,alpha0Gibbs, alphaGibbs, betaGibbs, lambdaGibbs, phiGibbs,xpost,ypost,c,c0,eps,gam_ab,mu_abold,sig_abold,lamb_abold)

varphiGibbs = [alphaGibbs;  betaGibbs];

cnd = true;
while cnd 
varphiStar = varphiGibbs + sqrt(lamb_abold)*chol(sig_abold,'lower')*randn(2,1);

alphaStar = varphiStar(1,:);
betaStar = varphiStar(2,:);

cnd = alphaStar+betaStar > eps || (0.01 > alphaStar) || (0.01 > betaStar);
end

prstar=log(dirchletpdf(varphiStar', c, c0));
prold=log(dirchletpdf(varphiGibbs', c, c0));

lklstar = lkl(Z,alpha0Gibbs, alphaStar, betaStar, lambdaGibbs,phiGibbs, xpost,ypost);
lklold = lkl(Z,alpha0Gibbs, alphaGibbs, betaGibbs, lambdaGibbs,phiGibbs, xpost,ypost);

ac = lklstar-lklold+prstar-prold;
rho_ab  = exp(min(0,ac));


u = rand(1,1);
if log(u) < min(0,ac) 
    varphiGibbsnew = varphiStar;
else
    varphiGibbsnew = varphiGibbs;
    rho_ab = 0;
end

alphaGibbsnew = varphiGibbsnew(1,:);
betaGibbsnew = varphiGibbsnew(2,:);


lamb_abnew = exp(log(lamb_abold) + gam_ab * (rho_ab - 0.34));
mu_abnew = mu_abold + gam_ab*(varphiGibbsnew - mu_abold);
sig_abnew = sig_abold + gam_ab*((varphiGibbsnew-mu_abold)*(varphiGibbsnew-mu_abold)'-sig_abold);

end