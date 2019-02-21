%/hpn/home/molinero/TCPSwitch/SJF-PS-simulations/get_average_service_time.pl for PARAM=0.995, rho=0.25, TYPE=BiModal, FAMILY=Exp_BMod_A_150
%A=150, B=270150, m1=1500, m2=3.6e+08, m3=98579119500000.1
clear all
hold off 

OUTPUT_FOLDER_IMAGES = './img/';

savepdf = true;
identifier   = ['^', 'o', 's'];
colors_mu    = ['r', 'g', 'b'];
colors_sigma = ['y', 'm', 'c'];
line_mu      = ['-' ];
line_sigma   = ['--'];
i = 1;
legd_txt =['', '', '','', '', ''];

% A=150; % mean service time of short jobs
% B=270150; % mean service time of long jobs (long tail)
A=54.13; % mean service time of short jobs
B=95.2; % mean service time of long jobs (long tail)
B_vec = [B];
alpha=0.995; % fraction of short jobs 
alpha_vec = [0.6, 0.8, 0.99];
rho_25_vec = [0.95, 0.8, 0.5]; % rho = lambda E(X) = system utilization
rho_ident_mu    = zeros(1,length(rho_25_vec));
rho_ident_sigma = zeros(1,length(rho_25_vec));
x_mu = 1;
x_sigma = 1;
k=[1:20];

txt = sprintf('\\begin{tabular}{|c|c|c|c|c|}'); disp(txt);
txt = sprintf('\\hline \n$E(X_S)/E(X_L)$ \t& $\\alpha$ \t& $\\rho$ \t& $min \\: k(\\mu)$  \t& $min \\: k(\\sigma)$ \t\\\\');
disp(txt);

for B = B_vec

    for alpha = alpha_vec

        f = figure('visible','off');

        j = 0;
        for rho = rho_25_vec
            j = j + 1;
            
            % Effect of min of average and standard deviation
            moment1 = alpha * A   + (1-alpha) * B  ; % mean service time  ============> E(X)
            moment2 = alpha * A^2 + (1-alpha) * B^2; % second moment of service time => E(X^2)
            moment3 = alpha * A^3 + (1-alpha) * B^3; % third moment of service time ==> E(X^3)

            
            rho_B= rho*(1-alpha)*B/moment1; 
            rho_A= rho*(alpha)  *A/moment1;
            
            Pblock  =1-poisscdf(k-2,rho_B.*k);
            Pblock_2=1-poisscdf(k-2,rho.^2.*k);
            Pblock_3=1-poisscdf(floor(k.*(1-rho_A)-1),rho_B.*k);
            %%%%%

            % equation (1)
            T  =Pblock  .*(rho./(1-rho).*(moment2)./2./(moment1))+(moment1).*k;
            %loglog(k,T_25  , 'DisplayName', sprintf('E(T) for \\rho=%3.2f',rho_25))
            idx_aux = 1 + mod(j-1, length(colors_mu));
            str_aux = strcat(colors_mu(idx_aux), '-');
            loglog(k, T, str_aux);

            axis auto;
            hold on;
            
            y = logspace(log10(T(x_mu)), log10(T(x_mu+1)), 10);
            rho_ident_mu(j) = y(3);
                

            %plot([x_mu, x_mu+1],[T_25(x_mu), T_25(x_mu+1)])
            %%%%%

            % equation (6), adapted to yield standard deviation rather than second moment
            T2  = Pblock  .*sqrt(rho./(1-rho).*moment3./(3.*moment1)) + sqrt(moment2).*k;
            %MT2  = Pblock  .* ((rho/(1-rho)) * (moment3 /(3 * moment1)) ) + moment2.* (k.^2);
            %T2 = MT2 - (T.^2);
            %loglog(k,T2  ,'--' , 'DisplayName', sprintf('\\sigma(T) for \\rho=%3.2f',rho))
            idx_aux = 1 + mod(j-1, length(colors_mu));
            str_aux = strcat(colors_mu(idx_aux), '--');
            loglog(k, T2, str_aux);
            
            y = logspace ( log10(T2(x_sigma)), log10(T2(x_sigma+1)), 10);
            rho_ident_sigma(j) = y(6);

            %plot([x_sigma, x_sigma+1],[T2_25(x_sigma), T2_25(x_sigma+1)])
            
            % print the table
            [T_min, T_idx] = min(T);
            [T2_min, T2_idx] = min(T2);
            k_min = k(T_idx);
            k2_min = k(T2_idx);
            
            %txt = sprintf('\\hline \n$%6.4f$ & $%4.2f$ \t& $%4.2f$ \t& $%3d$ \t& $%3d$ \t& $%3d$ \t\\\\',A/B, alpha, rho, k_min, k2_min, T2_min);
            txt = sprintf('\\hline \n$%6.4f$ & $%4.2f$ \t& $%4.2f$ \t& $%3d$ \t& $%3d$ \t\\\\',A/B, alpha, rho, k_min, k2_min);
            disp(txt);
            end

        str_file = strrep(sprintf('B_factor_%8.6f__alpha_%0.4f__rho_%0.4f_%0.4f', A/B, alpha, rho_25_vec(length(rho_25_vec)), rho_25_vec(1)), '.','_');
        str_title = sprintf('\\fontsize{10} \\fontname{Courier} \\alpha=%4.2f; \\rho=%4.2f, %4.2f, %4.2f;  (E(X_{s})/E(X_{l}))=%6.4f',alpha, rho_25_vec(1), rho_25_vec(2), rho_25_vec(3), A/B);
        %title([str_title],'Color','k')
        
        
        %lgd = legend('Location','southeast');
        for j = 0:(length(rho_25_vec) - 1)
            idx_aux = 1 + mod(j, length(colors_mu));
            str_aux = strcat(colors_mu(idx_aux), line_mu(1), identifier(idx_aux));
            test = 2 * j + 1;
            legd_plt(test) = plot(nan,nan,str_aux, 'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor',[0.01, 0.01, 0.01]);
            str_aux2 = sprintf('E(T) for \\rho=%3.2f',rho_25_vec(j+1));
            legd_txt{test} = str_aux2;
            
            idx_aux = 1 + mod(j, length(colors_mu));
            str_aux = strcat(colors_mu(idx_aux), '--', identifier(idx_aux));
            
            legd_plt(test+1) = plot(nan,nan,str_aux, 'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor',[0.01, 0.01, 0.01]);
            str_aux3 = sprintf('\\sigma(T) for \\rho=%3.2f',rho_25_vec(j+1));
            legd_txt{test+1} = str_aux3;
        end
        lgd = legend(legd_plt,legd_txt,'FontSize',14,'Location','northwest', 'Orientation','horizontal','AutoUpdate','off');
        lgd.NumColumns = 2;
        for i = [1:length(rho_25_vec)]
            i_aux = length(identifier);
            sz = 30;
            market1 = identifier(mod(i-1, length(identifier)) +1 );
            s1 = scatter(x_mu+0.2   , rho_ident_mu(i)   , sz, market1, 'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor',[0.01, 0.01, 0.01]);
            s2 = scatter(x_sigma+0.5, rho_ident_sigma(i), sz, market1, 'MarkerEdgeColor',[0, 0, 0], 'MarkerFaceColor',[0.01, 0.01, 0.01]);
            
            %s1.mkr = identifier(mod(i-1, length(identifier)) +1 );
            %s2.mkr = identifier(mod(i-1, length(identifier)) +1 );
        end
        
        %lgd.Title.String = 'average(cont.) and std deviation(dashed)';
        ylabel('response time');
        xlabel('number of servers (K)');

        %mkdir('tmp', 'Exp_BMod');
        file_name_prfx = sprintf('%s%s', OUTPUT_FOLDER_IMAGES, 'fit_v12_');
        if savepdf
            str_file_std = sprintf('%s%s.pdf', file_name_prfx ,str_file);
            fig = gcf;
            fig.PaperPositionMode = 'auto';
            fig_pos = fig.PaperPosition;
            fig.PaperSize = [fig_pos(3) fig_pos(4)];
            print (fig,str_file_std,'-dpdf');
        else
            str_file_std = sprintf('%s%s.png', file_name_prfx ,str_file);
            print (str_file_std,'-dpng');
        end
        hold off; 

    end
end
txt = sprintf('\\hline \n\\end{tabular}'); disp(txt);
