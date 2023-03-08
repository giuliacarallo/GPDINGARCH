function [xpost,ypost,ac,qstar,qold]=xySampl2New(Z,alpha0Gibbs, alphaGibbs, betaGibbs, lambdaGibbs, phiGibbs, x,y)
T = size(Z,1);
xpost = zeros(T,1);
ypost = zeros(T,1);
ac = zeros(T,1);
[~,~,~,~,thetaxt, thetayt] = mut(Z,alpha0Gibbs, alphaGibbs, betaGibbs, lambdaGibbs, phiGibbs);

%nu = mean((Z-mean(Z))>0);
nu = mean(Z>0);
%nu = (mean(Z>0)>0.5)*0.9+(mean(Z>0)<=0.5)*0.1;
for t=1:T
   u = rand(1,1);
   %thetaxtstar=thetaxt(t);
   %thetaytstar=thetayt(t);
   % thetaxtstar=min([thetaxt(t),abs(mean(abs(Z((Z-mean(Z))>0))))]);
   % thetaytstar=min([thetayt(t),abs(mean(abs(Z((Z-mean(Z))<0))))]);
     thetaxtstar=abs(mean(Z((Z-mean(Z))>0)));
     thetaytstar=abs(mean(Z((Z-mean(Z))<0)));

   %thetaxtstar=abs(mean(Z((Z-mean(Z))>0)));
   %thetaytstar=abs(mean(Z((Z-mean(Z))<0)));
   lambdastar = 0.6;
   
   if u <= nu
      %xstar = -10^5;
      %while xstar - Z(t) < 0
    xstar = gpbranching(thetaxtstar,lambdastar,1);
    ystar = xstar - Z(t,1);
      %end
   else
      %xstar = -10^5;
      %while xstar < 0
   ystar = gpbranching(thetaytstar,lambdastar,1);
      xstar = ystar + Z(t,1);
      %end
   end
   
   %    if ystar < 0
   %       ac(t,1) = 0;
   %       if u < min(1,ac(t,1))
   %          xpost(t,1) = xstar;
   %          ypost(t,1) = xstar - Z(t,1);
   %       else
   %       end
   %    else
   if (xstar - Z(t) < 0) || (xstar < 0)
      xpost(t,1)=x(t);
      ypost(t,1)=y(t);
   else
      
%      qstar = nu*((thetaxt(t,1)*(thetaxt(t,1)+lambdaGibbs*xstar)^(xstar-1)*exp(-xstar*lambdaGibbs-thetaxt(t,1)))/factorial(xstar))+ ...
%         (1-nu)*((thetayt(t,1)*(thetayt(t,1)+lambdaGibbs*(xstar-Z(t,1)))^(xstar-Z(t,1)-1)*exp(-(xstar-Z(t,1))*lambdaGibbs-thetaxt(t,1)))/factorial(xstar-Z(t,1)));
      
%      qold = nu*((thetaxt(t,1)*(thetaxt(t,1)+lambdaGibbs*x(t,1))^(x(t,1)-1)*exp(-x(t,1)*lambdaGibbs-thetaxt(t,1)))/factorial(x(t,1)))+ ...
%         (1-nu)*((thetayt(t,1)*(thetayt(t,1)+lambdaGibbs*(y(t)))^(y(t)-1)*exp(-(y(t))*lambdaGibbs-thetaxt(t,1)))/factorial(y(t)));

  %    qstar = nu*((thetaxtstar*(thetaxtstar+lambdastar*xstar)^(xstar-1)*exp(-xstar*lambdastar-thetaxtstar))/factorial(xstar))+ ...
  %      (1-nu)*((thetaytstar*(thetaytstar+lambdastar*(xstar-Z(t,1)))^(xstar-Z(t,1)-1)*exp(-(xstar-Z(t,1))*lambdastar-thetaxtstar))...
  %       /factorial(xstar-Z(t,1)));
     
       qstar =  nu*exp((log(thetaxtstar)+(xstar-1)*log(thetaxtstar+lambdastar*xstar)-thetaxtstar-lambdastar*xstar - gammaln(xstar+1)))+ ...
         (1-nu)*exp((log(thetaytstar)+(ystar-1)*log(thetaytstar+lambdastar*ystar)-thetaytstar-lambdastar*ystar - gammaln(ystar+1)));
     
    %  qold = nu*((thetaxtstar*(thetaxtstar+lambdastar*x(t,1))^(x(t,1)-1)*exp(-x(t,1)*lambdastar-thetaxtstar))/factorial(x(t,1)))+ ...
     %    (1-nu)*((thetaytstar*(thetaytstar+lambdastar*(y(t)))^(y(t)-1)*exp(-(y(t))*lambdastar-thetaxtstar))/factorial(y(t)));
     
        qold =  nu*exp((log(thetaxtstar)+(x(t,1)-1)*log(thetaxtstar+lambdastar*x(t,1))-thetaxtstar-lambdastar*x(t,1) - gammaln(x(t,1)+1)))+ ...
         (1-nu)*exp((log(thetaytstar)+(y(t,1)-1)*log(thetaytstar+lambdastar*y(t,1))-thetaytstar-lambdastar*y(t,1) - gammaln(y(t,1)+1)));

     lklstar = log(thetayt(t,1))+(xstar-Z(t,1)-1)*log(thetayt(t,1)+lambdaGibbs*(xstar-Z(t,1)))-thetayt(t,1)-lambdaGibbs*(xstar-Z(t,1)) - gammaln(xstar-Z(t,1)+1);
     lklold = log(thetayt(t,1))+(y(t,1)-1)*log(thetayt(t,1)+lambdaGibbs*y(t,1))-thetayt(t,1)-lambdaGibbs*y(t,1) - gammaln(y(t,1)+1); 
     
      ac(t) = exp(qold-qstar+lklstar-lklold);
      
     
     
      if rand(1) < min(1,ac(t))
         xpost(t,1) = xstar;
         ypost(t,1) = xstar - Z(t);
      else
         xpost(t,1)=x(t);
         ypost(t,1)=y(t);
      end
   end
end
end

