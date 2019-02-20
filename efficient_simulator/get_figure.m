%% Declared function to get the number of figure from paper
function  [num, letter] = get_figure(factor_El, alpha1)
    num = '\ref{fig:response_time}';
    letter ='';
    if abs(factor_El - 0.0005) < 0.0001
        if alpha1 == 0.99
            letter = 'a';
        elseif alpha1 == 0.8
            letter = 'b';
        elseif alpha1 == 0.6
            letter = 'c';
        end
    elseif abs(factor_El - 0.005) < 0.0001
        if alpha1 == 0.99
            letter = 'd';
        elseif alpha1 == 0.8
            letter = 'e';
        elseif alpha1 == 0.6
            letter = 'f';
        end
    elseif abs(factor_El - 0.05) < 0.0001
        if alpha1 == 0.99
            letter = 'g';
        elseif alpha1 == 0.8
            letter = 'h';
        elseif alpha1 == 0.6
            letter = 'i';
        end
    end

end
