function [coeff,norm2_res]=polminquad(x,y,n)
x = x';
y = y';
m = size(x);
m = m(1);
error = 'error';
if m ~= size(y)
    disp(error)
end

A = ones(m,1);
for i = 2:n
    A = [A,x.^(i-1)];
end

M = A;
R = zeros(n,n);
Q = zeros(m,n);
for k = 1:n
    R(k,k) = norm(A(:,k));
    Q(:,k) = (1/R(k,k))*A(:,k);
    for s = k+1:n
        R(k,s) = dot(Q(:,k),A(:,s));
        A(:,s) = A(:,s) - R(k,s)*Q(:,k);
    end
    A(:,k) = Q(:,k);
end
Residu = M - Q*R;
MSG = 'Error, la norma de un dels residu és molt gran!';

if norm(Residu) > 0.05 
    disp(MSG)
end

sol = zeros(n);
Qy = Q'*y;
sol = R\Qy;
Residu = M'*M*sol - M'*y;
if norm(Residu) > 0.05 
    disp(MSG)
end
coeff = sol;
norm2_res = norm(Residu);
end