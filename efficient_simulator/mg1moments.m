for myratio=[0.0005,0.005,0.05]
    for alpha=[0.99,0.8,0.6]
        for rho=[0.95,0.8,0.5]
            exs=54.13;
            exl=exs/myratio;
            % alpha=0.99;
            % rho=alpha*lambda*exs+(1-alpha)*lambda*exl
            lambda =rho/(alpha*exs+(1-alpha)*exl); 
            ex=alpha*exs+(1-alpha)*exl;
            mu=1/ex;
            ex2=alpha*exs*exs+(1-alpha)*exl*exl;
            ex3=alpha*exs*exs*exs+(1-alpha)*exl*exl*exl;
            ew=lambda*ex2/(2*(1-rho));
            ew2=2*ew*ew+lambda*ex3/(3*(1-rho));
            vt=ew2+ex2+2*ew*ex - (ew+ex)^2;
            fprintf('ratio: %f alpha: %f rho: %f E(T): %f sigma(T): %f\n',myratio, alpha,rho,ew+ex,sqrt(vt))
        end
    end
end

