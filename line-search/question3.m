%%%%%%%%%%%%%% Question3 %%%%%%%%%%%%%%
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

p = -A^(-1)*g;

funold = 100*(x(2)-x(1)^2)^2 + (1 - x(1))^2 + c*alpha*g'*p;
xnew = x + alpha*d;
funnew = 100*(xnew(2)-xnew(1)^2)^2 + (1 - xnew(1))^2;
x = xnew;

k = 1;

while (funnew > funold)
  xnew = x + alpha*d;
  dfdx1 = 200*(xnew(2)-xnew(1)^2)*(-2*xnew(1)) - 2*(1-xnew(1));
  dfdx2 = 100*(xnew(2)-xnew(1)^2);
  g = [dfdx1,dfdx2]';
  if (norm(g, inf) < 0.00001) break; end;
  d = -g; 
  A = eye(2);
  p = -A^(-1)*g;
  funold = 100*(x(2)-x(1)^2)^2 + (1 - x(1))^2 + c*alpha*g'*phi*p;
  xnew = x + alpha*d;
  funnew = 100*(xnew(2)-xnew(1)^2)^2 + (1 - xnew(1))^2;
  x = xnew;
  sprintf('k=%d fval=%f err=%f alpha=%f\n',k,funnew,norm(g, inf),alpha) 
  k = k + 1;
end
time = toc;