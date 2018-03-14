%%%%%%%%%%%%%% Question2 %%%%%%%%%%%%%%
tic;
alpha0 = 1.0;
n = 20;
x = rand(n,1);
f = fun(x); g = grad(x);
d = -g;
k = 1;

while ( 1== 1) 
  [retcode, x, f, g, alpha] = linesearch(x, f, g, d, alpha0);
  if (norm(g, inf) < 1e-2) break; end;
  d = -g; k = k + 1;
  sprintf('k=%d fval=%f err=%f alpha=%f\n',k,f,norm(g, inf),alpha) 
end
time = toc;

