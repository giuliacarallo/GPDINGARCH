
function likelihood= lkl(Z,alpha0, alpha, beta, lambda,phi, xpost,ypost)
T = size(Z,1);
[~,mutstar,~, sigma2tstar, thetaxt, thetayt] = mut(Z,alpha0, alpha, beta, lambda, phi);
P = zeros(T,1);
for t=1:T
   
   P(t,1) = log(thetaxt(t,1)) + (xpost(t,1)-1)*log(thetaxt(t,1)+lambda*xpost(t,1))-lambda*xpost(t,1)-thetaxt(t,1)-gammaln(xpost(t,1)+1)+...
      log(thetayt(t,1)) + (ypost(t,1)-1)*log(thetayt(t,1)+lambda*ypost(t,1))-lambda*ypost(t,1)-thetayt(t,1)-gammaln(ypost(t,1)+1);
   
%    P(t,1) = -sigma2tstar(t,1)+log(((sigma2tstar(t,1)^(2)- mutstar(t,1)^(2))/4))+(xpost(t,1)-1)*log((sigma2tstar(t,1)+mutstar(t,1))/2 + xpost(t,1)*lambda)...
%       +(ypost(t,1)-1)*log((sigma2tstar(t,1)-mutstar(t,1))/2 + ypost(t,1)*lambda)...
%       -lambda*(xpost(t,1)+ypost(t,1));%-log(factorial(xpost(t,1))*factorial(ypost(t,1)));
  

end

likelihood = sum(P);

end
