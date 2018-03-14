%%%%%%%%%%%%%% Question4 %%%%%%%%%%%%%%
tic;
phi = 0.9;
c = 0.00001;
alpha = 1;

A = eye(2);
x = [-1.2,1]';

dfdx1 = 200*(x(2)-x(1)^2)*(-2*x(1)) - 2*(1-x(1));
dfdx2 = 100*(x(2)-x(1)^2);
g = [dfdx1,dfdx2]';
d = -g;
f = 100*(x(2)-x(1)^2)^2 + (1 - x(1))^2;


k = 1;

while (1 == 1)
  [retcode, x, f, g, alpha] = linesearch(x, f, g, d, alpha);
  if (norm(g, inf) < 0.00001) break; end;
  d = -g; 
 % p = -A^(-1)*g;
 % funold = 100*(x(2)-x(1)^2)^2 + (1 - x(1))^2 + c*alpha*g'*phi*p;
 % xnew = x + alpha*d;
 % funnew = 100*(xnew(2)-xnew(1)^2)^2 + (1 - xnew(1))^2;
 % x = xnew;
  sprintf('k=%d fval=%f err=%f alpha=%f\n',k,funnew,norm(g, inf),alpha) 
  k = k + 1;
end
time = toc;