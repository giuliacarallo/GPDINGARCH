function mynumb = gpbranching(theta,lambda,n)
mynumb = zeros(1,n);
for i=1:n;
   index=0;
   y=poissrnd(theta);
   x=y;
   if y<=0 
      mynumb(i)= x;
   else while index==0;
         z=poissrnd(lambda*y);
         x=x+z;
         y=z;
         index=1*(y<=0);
         mynumb(i)=x;
      end
   end
   mynumb;
end