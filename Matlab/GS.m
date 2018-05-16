clear all

x = linspace(0,10,100); %Vector X
y = exp(x);            %Vector Y
n = 3; % Grau -1 del polinomi que aproxima
[coeff,norma] = polminquad(x,y,n);
X = linspace(0,10);
Y = 0;
for i = 1:n
    Y = Y + coeff(i)*X.^(i-1);
end
plot(x,y,'r',X,Y);
