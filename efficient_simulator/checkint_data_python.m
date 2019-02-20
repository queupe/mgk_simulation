clear all

colors_mu    = ['r', 'g', 'b'];

INPUT_FOLDER_CVS = './output/1M_samples/';
OUTPUT_FOLDER_PDF = './img/1M_samples_python/';

E_s = 54.13;
E_l_vec = [E_s/0.05, E_s/0.005, E_s/0.0005];
alpha_vec = [0.6, 0.8, 0.99];
rho_vec = [0.95, 0.8, 0.5];


%A = csvread('/home/vilc/Downloads/MDmk/output/100k_samples/fig_Bfactor_0.0050_alpha_0.99_rho_0.80.csv');

for E_l = E_l_vec
    for alpha1 = alpha_vec
        j = 0;
        for rho = rho_vec
            j = 1 + j;
            str_file = sprintf('%sfig_Bfactor_%6.4f_alpha_%4.2f_rho_%4.2f.csv', INPUT_FOLDER_CVS, E_s/E_l, alpha1, rho);

            disp(str_file)
            A = csvread(str_file);
            K = (A(:,1));

            M  = (A(:,2));
            Mp = (A(:,3));
            Mn = (A(:,4));

            S =  (A(:,5));
            Sp = (A(:,6));
            Serro = A(:,6) -A(:,5);
            Sn = A(:,5) - Serro;

            idx_aux = 1 + mod(j-1, length(colors_mu));
            str_aux = strcat(colors_mu(idx_aux), '-');
            semilogy (K, M, str_aux);
            hold on
            x_m = [K', flipud(K)'];
            y_m = [Mp', flipud(Mn)'];
            fill(x_m,y_m, 'b', 'edgecolor', 'none', 'facealpha', 0.25)
            hold on

            idx_aux = 1 + mod(j-1, length(colors_mu));
            str_aux = strcat(colors_mu(idx_aux), '--');
            semilogy (K, S, 'r--');
            hold on
            x_s = [K', flipud(K)'];
            y_s = [Sp', flipud(Sn)'];
            fill(x_s,y_s, 'r', 'edgecolor', 'none', 'facealpha', 0.25)
            hold on
            
            [MinM, iMinM] = min(M);
            [MinS, iMinS] = min(S);
            
            str = sprintf('Es/El:%7.5f, alpha:%4.2f, rho:%4.2f, min(E_t):%7.5f, Kmin(E_t):%3d, min(S_t):%7.5f, Kmin(S_t):%3d\n', E_s/E_l, alpha1, rho, MinM, iMinM, MinS, iMinS);
            disp(str)
        end
        ylabel('response time');
        xlabel('number of servers (K)');
        str_file = sprintf('B_factor_%8.6f__alpha_%4.2f__rho_%2.2fx%4.2f', E_s/E_l, alpha1, rho_vec(length(rho_vec)), rho_vec(1));
        str_file = strrep(str_file , '.','_');
        str_file_std = sprintf('%s%sb.pdf', OUTPUT_FOLDER_PDF ,str_file);
        fig = gcf;
        fig.PaperPositionMode = 'auto';
        fig_pos = fig.PaperPosition;
        fig.PaperSize = [fig_pos(3) fig_pos(4)];
        print (fig,str_file_std,'-dpdf');

        hold off
    end
end



