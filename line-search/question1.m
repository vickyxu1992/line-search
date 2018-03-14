%%%%%%%%%%%%%% Question1 %%%%%%%%%%%%%%
tic;
n = 20;
x = rand(n,1);
f = fun(x); g = grad(x);
d = -g;

n = 20;
A = zeros(n,n);
for i=1:n
  for j=1:n
    A(i,j) = 1.0/(i+j-1);
  end
end
alpha0 = (g'*g)/(g'*A*g);

k = 1;

while ( 1== 1) 
  [retcode, x, f, g, alpha] = linesearch(x, f, g, d, alpha0);
  if (norm(g, inf) < 1e-2) break; end;
  d = -g; k = k + 1;
  sprintf('k=%d fval=%f err=%f alpha=%f\n',k,f,norm(g, inf),alpha) 
end
time = toc;
