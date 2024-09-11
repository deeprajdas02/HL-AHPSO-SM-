function HL_AHPSO_SW(Runs,fhd,C,D,funcs,MaxFEs,NP,optimum)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%  HL-AHPSO-SM
    %%% ( Hybrid LShade Adaptive Heterogeneous Particle Swarm Optimisation with switiching mechanism )
    %%% Algorithm on CEC 2021 Benchmark Problems
    
    %%  Authors: Biraj Mahanta, Deepraj Das & Raj Roshan Singh
    %%  Under guidance : Somnath Mukhopadhaya
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %Rand_Seeds=load('input_data\Rand_Seeds.txt');

    Alg_Name=[ 'HL-AHPSO-SW_(' num2str(C(1)) num2str(C(2)) num2str(C(3)) ')']; 

    fprintf('Running %s algorithm on D= %d\n',Alg_Name, D)

    for n=0:15
        RecordFEsFactor(n+1) = round(D^((n/5)-3)*MaxFEs); 
    end

    progress = numel(RecordFEsFactor);
    
    for func_no = funcs

        % Combined paramters
        Tmax  = 40;
        lu = [-100 * ones(1, D); 100 * ones(1, D)];    % Max position
        
        % L-Shade parameters 
        H = 6;
        Co=0.5;
        Fo=0.5;
        
        % AHPSO parameters
        vmin = -1;      % Min velocity
        vmax = 1;       % Max velocity
        
        wmax = 0.9;     % Maximum inertia weight
        wmin = 0.2;     % Minimum inertia weight
        c = 0.08;        % Constant for adjusting inertia weight
        t_divisions = 10;

        fprintf('\n-------------------------------------------------------\n')
        fprintf('Function = %d, Dimension size = %d\n', func_no, D)
        allerrorvals = zeros(progress, Runs);

        for run_id = 1 : Runs

            FEs = 0;
            %run_seed=Rand_Seeds(1+mod(D*func_no*Runs+run_id-Runs,length(Rand_Seeds)));
            %rng(run_seed,'twister');
            %disp(run_seed);
            Run_RecordFEsFactor=RecordFEsFactor;
            run_funcvals = [];

            % initialize population and velocity for AHPSO
            X = repmat(lu(1, :), NP, 1) + rand(NP, D) .* (repmat(lu(2, :) - lu(1, :), NP, 1));
            V = repmat(vmin, NP, 1) + rand(NP, D) .* (repmat(vmax - vmin, NP, 1));
            FEs = NP;

            gbest_val = inf;
            gbest = X(1, :);
            
            phase = "AHPSO";

            while FEs < MaxFEs

                if phase == "AHPSO"
                    [gbest, gbest_val, FEs, X, pbest, phase] = AP   SO_SW(D, NP, Tmax, MaxFEs , FEs , X , V, lu(1, :), lu(2, :) , vmin, vmax ,wmin, wmax,t_divisions,c,phase,gbest,gbest_val,fhd, func_no, C);
                    
                    % Sort pbest based on the values in the last column (fitness values)
                    [~, sorted_indices] = sortrows(pbest, D + 1);
                    sorted_pbest = pbest(sorted_indices, :);
                    sorted_X = X(sorted_indices, :);
                    pbest = sorted_pbest;
                    X = sorted_X;
                end
                
                
                if phase == "LSHADE"    
                    % Select best performing subpopulation from HCLPSO
                    [subpopulation , N2] = selectSubpopulation(pbest(:,1:end-1), pbest(:,end));
                    
                    [gbest,gbest_val,FEs,X2,discarded_mem,phase] = L_SHADE(D,subpopulation,pbest(1:N2,end),Tmax,MaxFEs,FEs,N2,lu(1, :), lu(2, :),H,Co,Fo,phase,gbest,gbest_val, fhd, func_no, C);
                    
                    % caluclate number of subpopulation required 
                    numRequired = N2 - size(discarded_mem,1);
                    if FEs > MaxFEs / 2
                        % Use gbest mutation to generate the required number of solutions
                        mutatedSolutions = gbestMutation(gbest, D, numRequired);
                        X2 = [X2; mutatedSolutions];
                    else
                        % Repair discarded memory using gradient repair and add to the population
                        repairedMem = gradientRepair(discarded_mem, D, numRequired);
                        X2 = [X2; repairedMem];
                    end
                end
                X(1:N2,:) = X2(1:end,:);
                
                % Record performance
                if FEs >= Run_RecordFEsFactor(1)
                    run_funcvals = [run_funcvals; gbest_val];
                    Run_RecordFEsFactor(1) = [];
                end
            end

            if (C(1) == 1)
                run_funcvals = run_funcvals - optimum(func_no);
            end

            run_funcvals(run_funcvals < 1e-8) = 0;

            fprintf('%d th run, best-so-far error value = %1.8e\n', run_id, run_funcvals(end));
            allerrorvals(:, run_id) = run_funcvals;
        end

        fprintf('min_funvals:\t%e\n',min(allerrorvals(end,:)));
        fprintf('median_funvals:\t%e\n',median(allerrorvals(end,:)));
        fprintf('mean_funvals:\t%e\n',mean(allerrorvals(end,:)));
        fprintf('max_funvals:\t%e\n',max(allerrorvals(end,:)));
        fprintf('std_funvals:\t%e\n',std(allerrorvals(end,:)));
        
        file_name=sprintf('Results\\%s_%s_%s.txt',Alg_Name,int2str(func_no),int2str(D));
        save(file_name, 'allerrorvals', '-ascii');
    end
        

end