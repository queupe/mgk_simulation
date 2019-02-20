clear all

A = csvread('./output/100k_samples/B_factor_0.000500_alpha_0.9900_rho_0.8000_round_1.csv');
B = csvread('./output/100k_samples/B_factor_0.000500_alpha_0.9900_rho_0.8000_round_2.csv');
C = csvread('./output/100k_samples/B_factor_0.000500_alpha_0.9900_rho_0.8000_round_3.csv');
D = csvread('./output/100k_samples/B_factor_0.000500_alpha_0.9900_rho_0.8000_round_4.csv');
E = csvread('./output/100k_samples/B_factor_0.000500_alpha_0.9900_rho_0.8000_round_5.csv');

tam = 5;
M = (A(:,2)+B(:,2)+C(:,2)+D(:,2)+E(:,2))./tam;
S = (A(:,2).^2+B(:,2).^2+C(:,2).^2+D(:,2).^2+E(:,2).^2)./tam - M.^2;
S = sqrt(S);
error = (1.96 .* S) ./ sqrt(tam);

Mp = M + error;
Mn = M - error;


semilogy (A(:,1), M, 'b');
hold on
%loglog (A(:,1), Mn, 'y');
%hold on
%loglog (A(:,1), Mp, 'y');
%hold on
semilogy (A(:,1), A(:,4), 'r--');
hold on

X = [A(:,1)', A(:,1)']';
Y = [Mp', Mn']';

fill(X, Y, 'b')
alpha(0.25);
length(A(:,1))
length(Mp)
