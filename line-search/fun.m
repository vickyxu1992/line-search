function fval = fun(x)

n = 20;
A = zeros(n,n);
for i=1:n
  for j=1:n
    A(i,j) = 1.0/(i+j-1);
  end
end
b = ones(n,1);

fval = x'*(A*x)/2 - b'*x;

%tmp = inv(A)*b;
%val = tmp'*A*tmp/2 - b'*tmp
