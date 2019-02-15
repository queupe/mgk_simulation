clear all
hold off 

savepdf = true;
identifier   = ['^', 'o', 's'];
colors_an    = ['r', 'g', 'b']; % colors for analytical data
colors_si    = ['y', 'm', 'c']; % colors for simulation data
line_mu      = ['-' ];
line_sigma   = ['--'];
i = 1;
legd_txt =['', '', '','', '', ''];

PREFIX_CVS_FILES = '/home/vilc/Downloads/MDmk/output/100k_samples/';
PREFIX_IMG_FILES = '/home/vilc/Documents/UFRJ/git/microservice/img/analytical_model/fig_v12_';

Es        = 54.13; % mean service time of short jobs
El         = 95.2; % mean service time of long jobs (long tail)
El_vec     = [Es/0.05, Es/0.005, Es/0.0005];
alpha1    = 0.995; % fraction of short jobs 
alpha_vec = [0.6, 0.8, 0.99];
rho_vec   = [0.95, 0.8, 0.5]; % rho = lambda E(X) = system utilization
k         = [1:120]; % Number of cores in the processor

rho_ident_mu    = zeros(1,length(rho_vec));
rho_ident_sigma = zeros(1,length(rho_vec));
x_mu            = 1;
x_sigma         = 1;

%% Header of Latex table 
txt = sprintf('\\begin{tabular}{|c|c|c|c|c|c|c|}'); disp(txt);
txt = sprintf('\\hline \n$E(X_S)/E(X_L)$ \t& $\\alpha$ \t& $\\rho$ \t& $min \\: k(\\mu)$  \t& $min \\: k(\\sigma)$  \t& $min \\: k(\\mu)$  \t& $min \\: k(\\sigma)$ \t\\\\');
disp(txt);

%% Loops 
for El = El_vec

    for alpha1 = alpha_vec

        f = figure('visible','off');

        j = 0;
        for rho = rho_vec
            j = j + 1;
            
            % Effect of min of average and standard deviation
            moment1 = alpha1 * Es   + (1-alpha1) * El  ; % mean service time  ============> E(X)
            moment2 = alpha1 * Es^2 + (1-alpha1) * El^2; % second moment of service time => E(X^2)
            moment3 = alpha1 * Es^3 + (1-alpha1) * El^3; % third moment of service time ==> E(X^3)

            rho_l= rho*(1-alpha1)*El/moment1; 
            rho_s= rho*(alpha1)  *Es/moment1;
            
            Pblock  =1-poisscdf(k-2,rho_l.*k);
            %Pblock_2=1-poisscdf(k-2,rho.^2.*k);
            %Pblock_3=1-poisscdf(floor(k.*(1-rho_A)-1),rho_B.*k);
            %%%%%

            %% Analytical plotting
            % Mean
            % equation (1)
            T  = Pblock  .*(rho./(1-rho).*(moment2)./2./(moment1))+(moment1).*k;
            %loglog(k,T_25  , 'DisplayName', sprintf('E(T) for \\rho=%3.2f',rho_25))

            semilogy(k, T, 'c', 'DisplayName', 'Analytical E(T)');

            axis auto;
            hold on;

            % Standard deviation
            % equation (6), adapted to yield standard deviation rather than second moment
            T2  = Pblock  .*sqrt(rho./(1-rho).*moment3./(3.*moment1)) + sqrt(moment2).*k;
            %loglog(k,T2  ,'--' , 'DisplayName', sprintf('\\sigma(T) for \\rho=%3.2f',rho))
            loglog(k, T2, 'm--', 'DisplayName', 'Analytical \sigma(T)');
            hold on


            %% Simulation plotting
            file_simul = sprintf('fig_Bfactor_%6.4f_alpha_%4.2f_rho_%4.2f', Es/El, alpha1, rho);
            
            %disp(strcat(PREFIX_CVS_FILES,file_simul, '.csv'));
            A = csvread(strcat(PREFIX_CVS_FILES,file_simul, '.csv'));
            K = (A(:,1));

            M  = (A(:,2));
            Mp = (A(:,3));
            Mn = (A(:,4));

            S =  (A(:,5));
            Sp = (A(:,6));
            Serro = A(:,6) -A(:,5);
            Sn = A(:,5) - Serro;

            semilogy (K, M, 'b', 'DisplayName', 'Simulation E(T)');
            hold on

            semilogy (K, S, 'r--', 'DisplayName', 'Simulation \sigma(T)');
            hold on

            lgd = legend('FontSize',10,'Location','southeast', 'Orientation','horizontal','AutoUpdate','off');
            lgd.NumColumns = 1;
            
            x_m = [K', flipud(K)'];
            y_m = [Mp', flipud(Mn)'];
            fill(x_m,y_m, 'b', 'edgecolor', 'none', 'facealpha', 0.25)
            hold on
            
            x_s = [K', flipud(K)'];
            y_s = [Sp', flipud(Sn)'];
            fill(x_s,y_s, 'r', 'edgecolor', 'none', 'facealpha', 0.25)
            hold on
            

            ylabel('response time');
            xlabel('number of servers (K)');
            
            
            %% Save the figure
            str_file_pdf = sprintf('B_factor_%6.4f_alpha_%4.2f_rho_%4.2f', Es/El, alpha1, rho);
            %disp(str_file_pdf);
            str_title = sprintf('\\fontsize{10} \\fontname{Courier} \\alpha=%4.2f; \\rho=%4.2f;  (E(X_{s})/E(X_{l}))=%6.4f',alpha1, rho, Es/El);

            if savepdf
                str_file_std = strcat(PREFIX_IMG_FILES , file_simul, '.pdf');
                %disp(str_file_std)
                fig = gcf;
                fig.PaperPositionMode = 'auto';
                fig_pos = fig.PaperPosition;
                fig.PaperSize = [fig_pos(3) fig_pos(4)];
                print (fig,str_file_std,'-dpdf');
            else
                str_file_std = sprintf('%s%s.png', PREFIX_IMG_FILES ,str_file_pdf);
                print (str_file_std,'-dpng');
            end
            hold off;
            
            
            %% Latex table:
            [T_min, T_idx]   = min(T);
            [T2_min, T2_idx] = min(T2);
            k_min            = k(T_idx);
            k2_min           = k(T2_idx);
            
            [MinM, iMinM]    = min(M);
            [MinS, iMinS]    = min(S);
            
            txt = sprintf('\\hline \n$%6.4f$ & $%4.2f$ \t& $%4.2f$ \t& $%3d$ \t& $%3d$  \t& $%3d$ \t& $%3d$ \t\\\\',Es/El, alpha1, rho, k_min, k2_min, iMinM, iMinS);
            disp(txt);            

        end
    end
end
txt = sprintf('\\hline \n\\end{tabular}'); disp(txt);