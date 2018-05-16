% Gauss Seidel
function [x,iter]=seidel(A,b,x)
format long;	
tol=10^-12;
maxit=500;
D=diag(diag(A));
BGS=-inv(tril(A))*(triu(A,1));
normB=norm(BGS,2);
X1=sprintf('Norma sub-2 de B_GS: %f', normB);
disp(X1)
if (normB>=1) 
  disp('Radi espectral de B mes gran o igual que 1. Gauss Seidel no convergeix')
  return; 
end;
cGS=inv(tril(A))*b;
iter=0;
errabs=1;
x0=x;xnext=x0;
while ((abs(errabs)>tol) && (iter<maxit))
    x=xnext;		 
    iter=iter+1;
    xnext=BGS*x+cGS;
    errabs=norm(xnext-x,inf); 
end

XD=sprintf('\nEstimacio error absolut a la solucio iterativa: %e',errabs);
disp(XD) 
XD1=sprintf('\nIteracions usades: %d \n',iter);
disp(XD1) 
disp('Solucio aproximada usant Gauss Seidel:');
disp(xnext)
end


% Exemple de crida des de la finestra de comandes:
% clear all
% A = [1,0,0;0,1,0;1,0,131]
% x = [1,1,2]'
% b = [11,13,12]'
% [x,iter]=seidel(A,b,x)
