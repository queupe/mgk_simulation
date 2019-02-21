clear all
hold off 

savepdf = true;
table_values = true;
table_values_complete = true;
table_figures = false;

identifier   = ['^', 'o', 's'];
colors_an    = ['r', 'g', 'b']; % colors for analytical data
colors_si    = ['y', 'm', 'c']; % colors for simulation data
line_mu      = ['-' ];
line_sigma   = ['--'];
i = 1;
legd_txt =['', '', '','', '', ''];
index = 0;
value = 0;
str2 = '';

INPUT_FOLDER_CVS = './output/1M_samples/';
OUTPUT_FOLDER_PDF = './img/comp_analytical_simul/';


Es         = 54.13; % mean service time of short jobs
El         = 95.2;  % mean service time of long jobs (long tail)
El_vec     = [Es/0.0005, Es/0.005, Es/0.05];
%El_vec     = [Es/0.05];
alpha1    = 0.995; % fraction of short jobs 
alpha_vec = [0.99, 0.8, 0.6];
%alpha_vec = [0.99];
rho_vec   = [0.95, 0.8, 0.5]; % rho = lambda E(X) = system utilization
%rho_vec   = [0.95, 0.5]; % rho = lambda E(X) = system utilization
k         = [1:110]; % Number of cores in the processor

rho_ident_mu    = zeros(1,length(rho_vec));
rho_ident_sigma = zeros(1,length(rho_vec));
x_mu            = 1;
x_sigma         = 1;

%% Header of Latex table 
if table_values

    % Columns:  Es/El, \alpha, \rho, K_{\mu}^{\star}(simulation),
    %           K_{\sigma}^{\star}(simulation), \mu^{\star}, \sigma^{\star},
    %           \mu_{k=1}, \sigma_{k=1}
    
    if table_values_complete
        txt = sprintf('\\begin{tabular}{|c|c|c|c|c|rrr|rrr|r|r|}'); disp(txt);
        txt = sprintf('\\hline'); disp(txt);
        %txt = '$E_{S}/$ &          &        & \multicolumn{8}{c|}{analytical (simulation)} \\';
        %disp(txt);
        txt = '$E_{S}/$ &          &        & \multicolumn{6}{c|}{analytical (simulation)} & \multicolumn{2}{c|}{M/G/1}\\';
        disp(txt);
        txt = '$E_{L}$  & $\alpha$ & $\rho$ & $K^{\star}_{\mu}$ & $K^{\star}_{\sigma}$ & \multicolumn{3}{c|}{$\mu^{\star}$} & \multicolumn{3}{c|}{$\sigma^{\star}$} & \multicolumn{1}{c|}{$\mu_{K=1}$} & \multicolumn{1}{c|}{$\sigma_{K=1}$} \\';
        disp(txt);
    else
        txt = sprintf('\\begin{tabular}{|c|c|c|c|c|c|c|}'); disp(txt);
        txt = sprintf('\\hline'); disp(txt);
        txt = '$E_{S}/$ &          &        & \multicolumn{4}{c|}{analytical (simulation)} \\';
        disp(txt);
        txt = '$E_{L}$  & $\alpha$ & $\rho$ & $K^{\star}_{\mu}$ & $K^{\star}_{\sigma}$ & $\mu^{\star}$ & $\sigma^{\star}$ \\';
        disp(txt);        
    end


end

if table_figures
    txt_fig = '';
    stra_fig = '';
    strb_fig = '';
    strc_fig = '';
    txt = '\begin{figure*}[ht!]'; disp(txt);
    txt = '\begin{center}'; disp(txt);
    txt = '\begin{tabular}{@{\hspace{-1ex}}c @{\hspace{-3ex}}c@{\hspace{-3ex}}c}'; disp(txt);
end

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
            
            %disp(strcat(INPUT_FOLDER_CVS,file_simul, '.csv'));
            A = csvread(strcat(INPUT_FOLDER_CVS,file_simul, '.csv'));
            K = (A(:,1));

            M  = (A(:,2));
            Mp = (A(:,3));
            Mn = (A(:,4));
            Merro = A(:,3) - A(:,2);

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
            str_file_pdf = strrep(sprintf('fig_analy_simul_El_factor_%6.4f_alpha_%4.2f_rho_%4.2f', Es/El, alpha1, rho), '.','_');
            %disp(str_file_pdf);
            str_title = sprintf('\\fontsize{10} \\fontname{Courier} \\alpha=%4.2f; \\rho=%4.2f;  (E(X_{s})/E(X_{l}))=%6.4f',alpha1, rho, Es/El);
            %disp(str_file_pdf);
            if savepdf
                str_file_std = strcat(OUTPUT_FOLDER_PDF , strrep(str_file_pdf,'.','_'), '.pdf');
                %disp(str_file_std)
                fig = gcf;
                fig.PaperPositionMode = 'auto';
                fig_pos = fig.PaperPosition;
                fig.PaperSize = [fig_pos(3) fig_pos(4)];
                print (fig,str_file_std,'-dpdf');
            else
                str_file_std = sprintf('%s%s.png', OUTPUT_FOLDER_PDF ,str_file_pdf);
                print (str_file_std,'-dpng');
            end
            hold off;
            
            
            %% Latex table of values:
            if table_values
                [T_min, T_idx]   = min(T);
                [T2_min, T2_idx] = min(T2);
                k_min            = k(T_idx);
                k2_min           = k(T2_idx);

                [MinM, iMinM]    = min(M);
                [MinS, iMinS]    = min(S);
                
                if table_values_complete
                    % Columns: Es/El, \alpha, \rho,
                    % K^{\star}_{\mu}(simulation), K^{\star}_{\sigma}
                    % (simul), \mu^{\star}(simul), \sigma^{\star}(simul), 
                    % \mu_{K=1}, \sigma_{K=1}
                    txt = sprintf('\t\t\\hline'); disp(txt);
                    txt = sprintf('\t\t $%6.4f$ \t& $%4.2f$ \t& $%4.2f$ \t& $%d(%d)$ \t& $%d(%d)$ \t& $%7.2f($ \t& $%7.2f\\pm$ & $%0.2f)$ \t& $%7.2f($ \t& $%7.2f\\pm$ & $%0.2f)$ \t& %7.2f \t& %7.2f  \\\\', Es/El, alpha1, rho, k_min, iMinM, k2_min, iMinS, T_min, MinM, Merro(iMinM), T2_min, MinS, Serro(iMinS), T(1), T2(1));
                    disp(txt);
                    
                else
                    % Columns:  figure letter, \rho, K_{\mu}^{\star}(simulation),
                    %           K_{\sigma}^{\star}(simulation),
                    %           (\mu^{\star}), (\sigma^{\star})
                    [num, let] = get_figure(Es/El, alpha1);

                    txt = sprintf('\t\t\\hline'); disp(txt);
                    txt = sprintf('\t\t %s & $%4.2f$ \t& $%d(%d)$ \t& $%d$(%d) \t& $(%7.2f)$    \t& $(%7.2f)$ \t\\\\', let, rho, k_min, iMinM, k2_min, iMinS, MinM, MinS);
                    disp(txt);
                    %txt = sprintf('\t\t    &         \t&          \t&          \t& $(%7.2f)$ \t& $(%7.2f)$ \t\\\\', T_min, T2_min); 
                    %disp(txt);
                end


            end
            
            %% Latex table of figures
            if table_figures
                txt1 = sprintf('\t\t \\includegraphics[width=2.5in]{img/%s.pdf}', str_file_pdf);
                
                if index == 122
                    str2 = 'a';
                    value = - 25;
                end
                str1a = sprintf('\t\t (%s%c) $\\alpha=%4.2f$, $\\rho=%4.2f$',str2 , 97+index-value, alpha1, rho);
                str1b = sprintf('\t\t $E(X_S)/E(X_L)=%6.4f$, ' , Es/El);

                
                if mod(index,3) == 2
                    txt2 = sprintf('\\\\ \n');
                    txt_fig  = strcat(txt_fig,  txt1,  txt2);
                    stra_fig = strcat(stra_fig, str1a, txt2);
                    strb_fig = strcat(strb_fig, str1b, txt2);
                    disp(txt_fig);
                    disp(stra_fig); disp(strb_fig); disp(strc_fig);
                    txt_fig = '';
                    stra_fig = '';
                    strb_fig = '';
                    if mod(index,9) == 8
                        txt = '\end{tabular}'; disp(txt);
                        txt = '\end{center}'; disp(txt);
                        txt = '\caption{Response times as a function of  system parameters, for analytical and simulation analises. }'; disp(txt);
                        txt = '\label{fig:response_time}'; disp(txt);
                        txt = '\end{figure*}'; disp(txt);
                        %txt = sprintf('\n \n \%==============================\n');disp(txt);
                        
                        txt = '\begin{figure*}[ht!]'; disp(txt);
                        txt = '\begin{center}'; disp(txt);
                        txt = '\begin{tabular}{@{\hspace{-1ex}}c @{\hspace{-3ex}}c@{\hspace{-3ex}}c}'; disp(txt);
                        index = -1;
                    end
                else
                    txt2 = sprintf('& \n');
                    txt_fig = strcat(txt_fig, txt1, txt2);
                    stra_fig = strcat(stra_fig, str1a, txt2);
                    strb_fig = strcat(strb_fig, str1b, txt2);
                end
                index = index + 1;
            end

        end
    end
end

if table_values
    txt = sprintf('\\hline \n\\end{tabular}'); disp(txt);
end

if table_figures
    txt = '\end{tabular}'; disp(txt);
    txt = '\end{center}'; disp(txt);
    txt = '\caption{Response times as a function of  system parameters.  The number of servers that reduces the standard deviation of response times is typically smaller than or equal to the number of servers that reduces mean response times, indicating that to simplify intrusion detection one must trade between security and performance. }'; disp(txt);
    txt = '\label{fig:response_time}'; disp(txt);
    txt = '\end{figure*}'; disp(txt);
end

