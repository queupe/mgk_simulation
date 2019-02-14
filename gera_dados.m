clear all

E_s = 54.13;
E_l_vec = [E_s/0.05, E_s/0.005, E_s/0.0005];
alpha_vec = [0.6, 0.8, 0.99];
rho_vec = [0.95, 0.8, 0.5];

d = 0.1; % tempo médio entre chegadas

% change this variable to choose the number of samples
% Take careful: this can take a lot of time and space.
n_samples = 500;

for run = 1:5
    for E_l = E_l_vec
        for rho = rho_vec
            for alpha = alpha_vec
                
                lamb = rho / (alpha*E_s + (1-alpha)*E_l);
                m = [[1:n_samples]',cumsum(exprnd(1/lamb,n_samples,1)) , ((rand(n_samples,1) > alpha )+ones(n_samples,1))];
                str_file = strrep(sprintf('input_data/B_factor_%8.6f_alpha_%0.4f_rho_%0.4f_round_%d.csv', E_s/E_l, alpha, rho, run), '.','_');
        
                dlmwrite(str_file, m, 'precision', '%16.6f')

            end
        end
    end
end
