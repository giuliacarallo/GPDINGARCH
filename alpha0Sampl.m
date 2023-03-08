function [alpha0Gibbsnew, ac] = alpha0Sampl(Z,alpha0Gibbs, alphaGibbs, betaGibbs, lambdaGibbs, phiGibbs,xpost,ypost)

sc = 0.01;
alpha0star = alpha0Gibbs + randn(1,1)*sc;
%alpha0star = mean(Z) * (1-alphaGibbs-betaGibbs) + randn(1,1)*sc;

lklstar = lkl(Z,alpha0star, alphaGibbs, betaGibbs, lambdaGibbs,phiGibbs, xpost,ypost);
lklold = lkl(Z,alpha0Gibbs, alphaGibbs, betaGibbs, lambdaGibbs,phiGibbs, xpost,ypost);

ac = exp(lklstar-lklold);
%ac = exp(lklstar-lklold+log(pdf('normal',alpha0Gibbs,mean(Z) * (1-alphaGibbs-betaGibbs),sc))-log(pdf('normal',alpha0star,mean(Z) * (1-alphaGibbs-betaGibbs),sc)));

u = rand(1,1);
if u < min(1,ac)
   alpha0Gibbsnew = alpha0star;
else
   alpha0Gibbsnew = alpha0Gibbs;
end


end 