% %% HPSC Assignment 2 by 26AE91R03
% % ---- main Function -----
% clc;
% clear;
% close all;
% 
% % ============================================================
% % READ MATRIX
% % ============================================================
% 
% data = load('Kmat_1000.txt');
% 
% % Determine actual matrix size
% N = sqrt(length(data));
% 
% if N ~= floor(N)
%     error('Kmat does not contain N^2 entries.');
% end
% 
% N = round(N);
% 
% A = reshape(data,N,N);
% %A= 0.5.*(A+A');
% 
% % ============================================================
% % READ RHS
% % ============================================================
% 
% b = load('Fvec_1000.txt');
% b = b(:);
% 
% % ============================================================
% % CHECK
% % ============================================================
% 
% fprintf('\nMatrix size = %d x %d\n',N,N);
% fprintf('RHS size    = %d x %d\n',length(b),1);
% fprintf('NNZ         = %d\n',nnz(A));
% 
% % ============================================================
% % SOLVER PARAMETERS
% % ============================================================
% 
% tol = 1e-10;
% maxIter = 200000;
% 
% % ============================================================
% % SOR
% % ============================================================
% 
% [x_sor,iter_sor,elapsed_sor,res_sor,omega_sor] = ...
%     mySOR(A,b,tol,maxIter);
% 
% % RESULTS
% % ============================================================
% 
% fprintf('\n====================================\n');
% fprintf('SOR RESULTS\n');
% fprintf('====================================\n');
% 
% fprintf('Matrix size       = %d x %d\n',N,N);
% fprintf('Omega optimum     = %.8f\n',omega_sor);
% fprintf('Iterations        = %d\n',iter_sor);
% fprintf('Computing time    = %.6f seconds\n',elapsed_sor);
% fprintf('Final residual    = %.6e\n',res_sor);
% 
% %%  CG METHOD
% [xCG,iterCG,timeCG,resCG,resvecCG] = ...
%     CG(A,b,tol,maxIter);
% 
% fprintf('\n====================================\n');
% fprintf('CG RESULTS\n');
% fprintf('====================================\n');
% 
% fprintf('Iterations      = %d\n',iterCG);
% fprintf('Time            = %.8f seconds\n',timeCG);
% fprintf('Final residual  = %.12e\n',resCG);
% 
% 
% 
% 
% 
% 
% %%  FOR MR
% [xMR,iterMR,timeMR,resMR,resvecMR] = ...
%     MR(A,b,tol,maxIter);
% 
% fprintf('\n====================================\n');
% fprintf('MR RESULTS\n');
% fprintf('====================================\n');
% 
% fprintf('Iterations     = %d\n',iterMR);
% fprintf('Time            = %.8f seconds\n',timeMR);
% fprintf('Final residual  = %.12e\n',resMR);
% 
% 
% 
% 
% %% Steepest Descent SOLVE
% % ============================================================
% 
% [x_sd,iter_sd,elapsed_sd,res_sd,resvec_sd] = ...
%     steepest_descent(A,b,tol,maxIter);
% 
% % RESULTS
% % ============================================================
% fprintf('\n====================================\n');
% fprintf('Steepest Descent RESULTS\n');
% fprintf('====================================\n');
% 
% fprintf('\nIterations     = %d\n',iter_sd);
% fprintf('Time           = %.6f seconds\n',elapsed_sd);
% fprintf('Final residual = %.6e\n',res_sd);
% 
% % ============================================================
% % PLOT CONVERGENCE
% % ============================================================
% 
% figure;
% 
% semilogy(0:iter_sd,resvec_sd,'LineWidth',1.5);
% 
% grid on;
% 
% xlabel('Iteration');
% ylabel('Relative Residual');
% 
% title('Steepest Descent Convergence');
% 
% 
% %% BICGSTAB -----
% [xBi,iterBi,timeBi,resBi,resvecBi] = ...
%     BiCGSTAB(A,b,tol,maxIter);
% 
% fprintf('\n====================================\n');
% fprintf('BICGSTAB RESULTS\n');
% fprintf('====================================\n');
% 
% fprintf('Iterations      = %d\n',iterBi);
% fprintf('Time            = %.8f seconds\n',timeBi);
% fprintf('Final residual  = %.12e\n',resBi);
% 
% 
% 
% %% ****************************  HELPER FUNCTION ***************
% 
% %% 
% %--- CG METHOD 
% function [x,iter,elapsed,res,resvec] = CG(A,b,tol,maxIter)
% 
% % ============================================================
% % CONJUGATE GRADIENT METHOD
% %
% % Solves:
% %              A*x = b
% %
% % CG requires A to be symmetric positive definite.
% % ============================================================
% 
% N = length(b);
% 
% % Initial guess
% x = zeros(N,1);
% 
% % Initial residual
% r = b - A*x;
% 
% % Initial search direction
% p = r;
% 
% % Normalization
% normb = max(norm(b),1);
% 
% % Initial residual
% res = norm(r)/normb;
% 
% % Residual history
% resvec = zeros(maxIter+1,1);
% resvec(1) = res;
% 
% fprintf('\nInitial residual = %.12e\n',res);
% 
% % Already converged?
% if res <= tol
% 
%     iter = 0;
%     elapsed = 0;
%     resvec = resvec(1);
% 
%     return;
% 
% end
% 
% tic;
% 
% % ============================================================
% % CG ITERATIONS
% % ============================================================
% 
% for iter = 1:maxIter
% 
%     % A*p
%     Ap = A*p;
% 
%     % p^T A p
%     denominator = p' * Ap;
% 
%     % Breakdown / non-SPD check
%     if denominator <= 0 || abs(denominator) < 1e-30
% 
%         warning(['CG breakdown: A may not be ' ...
%                  'symmetric positive definite.']);
% 
%         break;
% 
%     end
% 
%     % r^T r
%     rr_old = r' * r;
% 
%     % Step size
%     alpha = rr_old / denominator;
% 
%     % Update solution
%     x = x + alpha*p;
% 
%     % Update residual
%     r_new = r - alpha*Ap;
% 
%     % Relative residual
%     res = norm(r_new)/normb;
% 
%     resvec(iter+1) = res;
% 
%     % Print first few iterations
%     if iter <= 10
% 
%         fprintf(['Iteration %4d : alpha = %.12e, ' ...
%                  'residual = %.12e\n'], ...
%                  iter,alpha,res);
% 
%     end
% 
%     % Convergence
%     if res <= tol
% 
%         r = r_new;
% 
%         fprintf('\nCG converged at iteration %d\n',iter);
% 
%         break;
% 
%     end
% 
%     % Beta
%     rr_new = r_new' * r_new;
% 
%     beta = rr_new / rr_old;
% 
%     % New search direction
%     p = r_new + beta*p;
% 
%     % Update residual
%     r = r_new;
% 
% end
% 
% elapsed = toc;
% 
% resvec = resvec(1:iter+1);
% 
% end
% 
% 
% 
% %-- SOR function ----
% function [x,iter,elapsed,res,omega] = mySOR(A,b,tol,maxIter)
% 
% % ============================================================
% % SOR SOLVER
% % ============================================================
% 
% N = length(b);
% 
% % Initial guess
% x = zeros(N,1);
% 
% % ------------------------------------------------------------
% % Jacobi iteration matrix
% % ------------------------------------------------------------
% 
% D = diag(diag(A));
% 
% BJ = -D\(A-D);
% 
% rhoJ = max(abs(eig(BJ)));
% 
% % ------------------------------------------------------------
% % Optimum omega
% % ------------------------------------------------------------
% 
% if rhoJ < 1
% 
%     omega = 2/(1 + sqrt(1-rhoJ^2));
% 
% else
% 
%     warning('Jacobi spectral radius >= 1');
% 
%     omega = 1.0;
% 
% end
% 
% fprintf('Jacobi spectral radius = %.8f\n',rhoJ);
% fprintf('Optimum omega = %.8f\n',omega);
% 
% % ------------------------------------------------------------
% % Initial residual
% % ------------------------------------------------------------
% 
% normb = max(norm(b),1);
% 
% tic;
% 
% % ------------------------------------------------------------
% % SOR ITERATION
% % ------------------------------------------------------------
% 
% for iter = 1:maxIter
% 
%     for i = 1:N
% 
%         % Contribution from already updated values
%         if i > 1
%             sigma1 = A(i,1:i-1)*x(1:i-1);
%         else
%             sigma1 = 0;
%         end
% 
%         % Contribution from old values
%         if i < N
%             sigma2 = A(i,i+1:N)*x(i+1:N);
%         else
%             sigma2 = 0;
%         end
% 
%         % Gauss-Seidel value
%         x_GS = (b(i)-sigma1-sigma2)/A(i,i);
% 
%         % SOR update
%         x(i) = (1-omega)*x(i) + omega*x_GS;
% 
%     end
% 
%     % Residual
%     r = b-A*x;
% 
%     res = norm(r)/normb;
% 
%     % Convergence
%     if res <= tol
%         break;
%     end
% 
% end
% 
% elapsed = toc;
% 
% end
% 
% 
% % --- STEEPEST DESCENT ------%
% function [x,iter,elapsed,res,resvec] = steepest_descent(A,b,tol,maxIter)
% 
% N = length(b);
% 
% % Initial guess
% x = zeros(N,1);
% 
% % Normalization
% normb = max(norm(b),1);
% 
% % Residual
% r = b - A*x;
% 
% % Store residual history
% resvec = zeros(maxIter+1,1);
% 
% resvec(1) = norm(r)/normb;
% 
% tic;
% 
% % ------------------------------------------------------------
% % Iteration
% % ------------------------------------------------------------
% 
% for iter = 1:maxIter
% 
%     % A*r
%     Ar = A*r;
% 
%     % Denominator
%     denominator = r' * Ar;
% 
%     % Check breakdown
%     if abs(denominator) < 1e-30
%         warning('Steepest Descent breakdown.');
%         break;
%     end
% 
%     % Optimum step length
%     alpha = (r' * r) / denominator;
% 
%     % Update solution
%     x = x + alpha*r;
% 
%     % New residual
%     r = b - A*x;
% 
%     % Relative residual
%     res = norm(r)/normb;
% 
%     resvec(iter+1) = res;
% 
%     % Convergence
%     if res <= tol
%         res
%         break;
%     end
% 
% end
% 
% elapsed = toc;
% 
% % Remove unused entries
% resvec = resvec(1:iter+1);
% 
% end
% 
% % ------------ MIM RESIDUAL METHOD --------
% function [x,iter,elapsed,res,resvec] = MR(A,b,tol,maxIter)
% 
% % ============================================================
% % MINIMUM RESIDUAL (MR) METHOD
% %
% % Solves A*x = b
% % ============================================================
% 
% N = length(b);
% 
% % Initial guess
% x = zeros(N,1);
% 
% % Initial residual
% r = b - A*x;
% 
% % Residual normalization
% normb = max(norm(b),1);
% 
% res = norm(r)/normb;
% 
% % Store residual history
% resvec = zeros(maxIter+1,1);
% resvec(1) = res;
% 
% fprintf('\nInitial residual = %.12e\n',res);
% 
% % ------------------------------------------------------------
% % If already converged
% % ------------------------------------------------------------
% 
% if res <= tol
% 
%     iter = 0;
%     elapsed = 0;
%     resvec = resvec(1);
% 
%     return;
% 
% end
% 
% tic;
% 
% % ============================================================
% % MR ITERATIONS
% % ============================================================
% 
% for iter = 1:maxIter
% 
%     % q = A*r
%     q = A*r;
% 
%     % Numerator
%     numerator = q' * r;
% 
%     % Denominator
%     denominator = q' * q;
% 
%     % Check breakdown
%     if abs(denominator) < 1e-30
% 
%         warning('MR breakdown.');
% 
%         break;
% 
%     end
% 
%     % MR step length
%     alpha = numerator / denominator;
% 
%     % Update solution
%     x = x + alpha*r;
% 
%     % Update residual
%     r = r - alpha*q;
% 
%     % Relative residual
%     res = norm(r)/normb;
% 
%     % Store
%     resvec(iter+1) = res;
% 
%     % Display first few iterations
%     if iter <= 10
% 
%         fprintf(['Iteration %4d : alpha = %.12e, ' ...
%                  'residual = %.12e\n'], ...
%                  iter,alpha,res);
% 
%     end
% 
%     % Convergence
%     if res <= tol
% 
%         fprintf('\nMR converged at iteration %d\n',iter);
% 
%         break;
% 
%     end
% 
% end
% 
% elapsed = toc;
% 
% % Remove unused values
% resvec = resvec(1:iter+1);
% 
% end
% 
% 
% % -- BICGSTAB ----
% function [x,iter,elapsed,res,resvec] = BiCGSTAB(A,b,tol,maxIter)
% 
% % ============================================================
% % BICGSTAB METHOD
% %
% % Solves:
% %              A*x = b
% %
% % Suitable for nonsymmetric matrices.
% % ============================================================
% 
% N = length(b);
% 
% % ------------------------------------------------------------
% % Initial guess
% % ------------------------------------------------------------
% 
% x = zeros(N,1);
% 
% % Initial residual
% r = b - A*x;
% 
% % Shadow residual
% rhat = r;
% 
% % Normalization
% normb = max(norm(b),1);
% 
% % Initial residual
% res = norm(r)/normb;
% 
% % Residual history
% resvec = zeros(maxIter+1,1);
% resvec(1) = res;
% 
% fprintf('\nInitial residual = %.12e\n',res);
% 
% if res <= tol
% 
%     iter = 0;
%     elapsed = 0;
%     resvec = resvec(1);
% 
%     return;
% 
% end
% 
% % ------------------------------------------------------------
% % Initial variables
% % ------------------------------------------------------------
% 
% p = zeros(N,1);
% v = zeros(N,1);
% 
% rho_old = 1;
% alpha = 1;
% omega = 1;
% 
% tic;
% 
% % ============================================================
% % BICGSTAB ITERATIONS
% % ============================================================
% 
% for iter = 1:maxIter
% 
%     % --------------------------------------------------------
%     % rho
%     % --------------------------------------------------------
% 
%     rho = rhat' * r;
% 
%     if abs(rho) < 1e-30
% 
%         warning('BiCGSTAB breakdown: rho = 0');
% 
%         break;
% 
%     end
% 
%     % --------------------------------------------------------
%     % beta
%     % --------------------------------------------------------
% 
%     if iter == 1
% 
%         beta = 0;
% 
%     else
% 
%         beta = (rho/rho_old)*(alpha/omega);
% 
%     end
% 
%     % --------------------------------------------------------
%     % Search direction
%     % --------------------------------------------------------
% 
%     p = r + beta*(p - omega*v);
% 
%     % --------------------------------------------------------
%     % v = A*p
%     % --------------------------------------------------------
% 
%     v = A*p;
% 
%     % --------------------------------------------------------
%     % alpha
%     % --------------------------------------------------------
% 
%     denominator = rhat' * v;
% 
%     if abs(denominator) < 1e-30
% 
%         warning('BiCGSTAB breakdown: denominator = 0');
% 
%         break;
% 
%     end
% 
%     alpha = rho/denominator;
% 
%     % --------------------------------------------------------
%     % Intermediate residual
%     % --------------------------------------------------------
% 
%     s = r - alpha*v;
% 
%     % --------------------------------------------------------
%     % Check convergence before second step
%     % --------------------------------------------------------
% 
%     res_s = norm(s)/normb;
% 
%     if res_s <= tol
% 
%         x = x + alpha*p;
% 
%         res = res_s;
% 
%         resvec(iter+1) = res;
% 
%         fprintf('\nBiCGSTAB converged at iteration %d\n',iter);
% 
%         break;
% 
%     end
% 
%     % --------------------------------------------------------
%     % t = A*s
%     % --------------------------------------------------------
% 
%     t = A*s;
% 
%     % --------------------------------------------------------
%     % omega
%     % --------------------------------------------------------
% 
%     tt = t' * t;
% 
%     if abs(tt) < 1e-30
% 
%         warning('BiCGSTAB breakdown: t^T*t = 0');
% 
%         break;
% 
%     end
% 
%     omega = (t' * s)/tt;
% 
%     % --------------------------------------------------------
%     % Update solution
%     % --------------------------------------------------------
% 
%     x = x + alpha*p + omega*s;
% 
%     % --------------------------------------------------------
%     % Update residual
%     % --------------------------------------------------------
% 
%     r = s - omega*t;
% 
%     % --------------------------------------------------------
%     % Relative residual
%     % --------------------------------------------------------
% 
%     res = norm(r)/normb;
% 
%     resvec(iter+1) = res;
% 
%     % --------------------------------------------------------
%     % Print first 10 iterations
%     % --------------------------------------------------------
% 
%     if iter <= 10
% 
%         fprintf(['Iteration %4d : alpha = %.12e, ' ...
%                  'omega = %.12e, residual = %.12e\n'], ...
%                  iter,alpha,omega,res);
% 
%     end
% 
%     % --------------------------------------------------------
%     % Convergence
%     % --------------------------------------------------------
% 
%     if res <= tol
% 
%         fprintf('\nBiCGSTAB converged at iteration %d\n',iter);
% 
%         break;
% 
%     end
% 
%     % --------------------------------------------------------
%     % Prepare next iteration
%     % --------------------------------------------------------
% 
%     rho_old = rho;
% 
% end
% 
% elapsed = toc;
% 
% resvec = resvec(1:iter+1);
% 
% end



%% PART II FOR THE ASSIGNMENT 3 
clc; clear all; 
%% PARAMETERS

matrixSizes = [100, 1000, 10000];
numWorkers  = [1, 2, 4, 6];

tol = 1e-10;
maxIter = 200000;
numRuns = 3;

results = [];

%% ============================================================
% MAIN EXPERIMENT
%% ============================================================

for ns = 1:length(matrixSizes)

    N_expected = matrixSizes(ns);

    fprintf('\n\n========================================\n');
    fprintf('MATRIX ORDER = %d\n',N_expected);
    fprintf('========================================\n');

    %% READ MATRIX

    matrixFile = sprintf('Kmat_%d.txt',N_expected);
    rhsFile = sprintf('Fvec_%d.txt',N_expected);

    data = load(matrixFile);

    N = sqrt(length(data));

    % if N ~= floor(N)
    %     error('Invalid matrix data.');
    % end
    % 
    % N = round(N);
    % 
    % if N ~= N_expected
    %     error('Unexpected matrix dimension.');
    % end

    % MATLAB reshape uses column-major ordering
    A = reshape(data,N,N);

    % Construct structured symmetric matrix
    A = 0.5*(A+A');

    % Convert to sparse matrix
    A = sparse(A);

    % Read RHS
    b = load(rhsFile);
    b = b(:);

    if length(b) ~= N
        error('Matrix and RHS dimensions mismatch.');
    end

    fprintf('Matrix size = %d x %d\n',N,N);
    fprintf('Number of nonzeros = %d\n',nnz(A));

    %% MATRIX CHECKS

    symmetryError = norm(A-A','fro') / ...
                    max(norm(A,'fro'),eps);

    fprintf('Symmetry error = %.6e\n',symmetryError);

    if any(abs(diag(A)) < 1e-30)
        error('Matrix has a zero diagonal element.');
    end

    %% DIAGONAL EXTRACTION

    D = diag(A);

    %% STORAGE FOR SPEEDUP

    timeCG = zeros(length(numWorkers),1);
    timeJacobi = zeros(length(numWorkers),1);

    %% ========================================================
    % LOOP OVER NUMBER OF WORKERS
    %% ========================================================

    for ip = 1:length(numWorkers)

        p = numWorkers(ip);

        fprintf('\n----------------------------------------\n');
        fprintf('NUMBER OF WORKERS = %d\n',p);
        fprintf('----------------------------------------\n');

        % Close any previously existing pool
        pool = gcp('nocreate');

        if ~isempty(pool)
            delete(pool);
        end

        % Create thread-based pool for parallel runs
        if p > 1

            try
                parpool("Threads",p);
            catch ME
                error(['Cannot create a thread pool. ', ...
                       'Check Parallel Computing Toolbox. ', ...
                       'MATLAB error: ',ME.message]);
            end

            useParallel = true;

        else
            % Serial baseline
            useParallel = false;
        end

        %% ----------------------------------------------------
        % CG WARM-UP
        %% ----------------------------------------------------

        [~,~,~,~] = CG_solver(A,b,tol,maxIter,useParallel);

        %% CG EXECUTION

        cgTimes = zeros(numRuns,1);

        for run = 1:numRuns

            [xCG,iterCG,tCG,resCG] = ...
                CG_solver(A,b,tol,maxIter,useParallel);

            cgTimes(run) = tCG;

        end

        timeCG(ip) = median(cgTimes);

        fprintf('\nCG RESULTS\n');
        fprintf('Iterations = %d\n',iterCG);
        fprintf('Time       = %.8f seconds\n',timeCG(ip));
        fprintf('Residual   = %.6e\n',resCG);

        %% ----------------------------------------------------
        % JACOBI WARM-UP
        %% ----------------------------------------------------

        [~,~,~,~] = Jacobi_solver( ...
            A,b,D,tol,maxIter,useParallel);

        %% JACOBI EXECUTION

        jacobiTimes = zeros(numRuns,1);

        for run = 1:numRuns

            [xJac,iterJac,tJac,resJac] = ...
                Jacobi_solver(A,b,D,tol,maxIter,useParallel);

            jacobiTimes(run) = tJac;

        end

        timeJacobi(ip) = median(jacobiTimes);

        fprintf('\nJACOBI RESULTS\n');
        fprintf('Iterations = %d\n',iterJac);
        fprintf('Time       = %.8f seconds\n',timeJacobi(ip));
        fprintf('Residual   = %.6e\n',resJac);

        %% CLOSE POOL

        pool = gcp('nocreate');

        if ~isempty(pool)
            delete(pool);
        end

    end

    %% ========================================================
    % SPEEDUP CALCULATION
    %% ========================================================

    speedupCG = timeCG(1)./timeCG;

    speedupJacobi = timeJacobi(1)./timeJacobi;

    efficiencyCG = speedupCG./numWorkers(:)*100;

    efficiencyJacobi = speedupJacobi./numWorkers(:)*100;

    %% ========================================================
    % DISPLAY RESULTS
    %% ========================================================

    fprintf('\n\n========================================\n');
    fprintf('FINAL RESULTS FOR MATRIX ORDER %d\n',N);
    fprintf('========================================\n');

    T = table(numWorkers(:),timeCG,speedupCG, ...
        efficiencyCG,timeJacobi,speedupJacobi, ...
        efficiencyJacobi, ...
        'VariableNames', ...
        {'Workers','CG_Time','CG_Speedup','CG_Efficiency', ...
         'Jacobi_Time','Jacobi_Speedup','Jacobi_Efficiency'});

    disp(T);

    %% SAVE RESULTS

    for ip = 1:length(numWorkers)

        results = [results;
            N,numWorkers(ip), ...
            timeCG(ip),speedupCG(ip),efficiencyCG(ip), ...
            timeJacobi(ip),speedupJacobi(ip), ...
            efficiencyJacobi(ip)];

    end

    %% ========================================================
    % CG SPEEDUP PLOT
    %% ========================================================

    figure;

    plot(numWorkers,speedupCG,'-o', ...
        'LineWidth',1.8,'MarkerSize',7);

    hold on;

    plot(numWorkers,numWorkers,'--k','LineWidth',1.2);

    grid on;

    xlabel('Number of Workers');
    ylabel('Speedup');

    title(sprintf('CG Speedup - Matrix Order %d',N));

    legend('Actual Speedup','Ideal Speedup', ...
        'Location','northwest');

    xticks(numWorkers);

    %% ========================================================
    % JACOBI SPEEDUP PLOT
    %% ========================================================

    figure;

    plot(numWorkers,speedupJacobi,'-s', ...
        'LineWidth',1.8,'MarkerSize',7);

    hold on;

    plot(numWorkers,numWorkers,'--k','LineWidth',1.2);

    grid on;

    xlabel('Number of Workers');
    ylabel('Speedup');

    title(sprintf('Jacobi Speedup - Matrix Order %d',N));

    legend('Actual Speedup','Ideal Speedup', ...
        'Location','northwest');

    xticks(numWorkers);

end

%% ============================================================
% SAVE ALL RESULTS
%% ============================================================

FinalResults = array2table(results, ...
    'VariableNames', ...
    {'MatrixOrder','Workers', ...
     'CG_Time','CG_Speedup','CG_Efficiency', ...
     'Jacobi_Time','Jacobi_Speedup','Jacobi_Efficiency'});

writetable(FinalResults,'OpenMP_MATLAB_Results.csv');

fprintf('\n\n========================================\n');
fprintf('ALL COMPUTATIONS COMPLETED\n');
fprintf('Results saved in OpenMP_MATLAB_Results.csv\n');
fprintf('========================================\n');

disp(FinalResults);


%% ============================================================
% FUNCTION 1: CONJUGATE GRADIENT
%% ============================================================

function [x,iter,elapsed,res] = ...
    CG_solver(A,b,tol,maxIter,useParallel)

    N = length(b);

    x = zeros(N,1);

    r = b-A*x;

    p = r;

    normb = max(norm(b),1e-30);

    rr = r'*r;

    res = norm(r)/normb;

    if res <= tol
        iter = 0;
        elapsed = 0;
        return;
    end

    tic;

    for iter = 1:maxIter

        % Matrix-vector multiplication
        Ap = parallelMatVec(A,p,useParallel);

        denominator = p'*Ap;

        if denominator <= 1e-30
            warning('CG breakdown: Matrix may not be SPD.');
            break;
        end

        alpha = rr/denominator;

        % Update solution and residual
        x = x + alpha*p;

        rnew = r-alpha*Ap;

        rrnew = rnew'*rnew;

        res = sqrt(rrnew)/normb;

        if res <= tol
            r = rnew;
            break;
        end

        beta = rrnew/rr;

        p = rnew+beta*p;

        r = rnew;

        rr = rrnew;

    end

    elapsed = toc;

    % Calculate actual final residual
    rFinal = b-parallelMatVec(A,x,useParallel);

    res = norm(rFinal)/normb;

end


%% ============================================================
% FUNCTION 2: JACOBI ITERATION
%% ============================================================

function [x,iter,elapsed,res] = ...
    Jacobi_solver(A,b,D,tol,maxIter,useParallel)

    N = length(b);

    x = zeros(N,1);

    xnew = zeros(N,1);

    normb = max(norm(b),1e-30);

    res = norm(b)/normb;

    if res <= tol
        iter = 0;
        elapsed = 0;
        return;
    end

    tic;

    for iter = 1:maxIter

        % Jacobi update
        % Every row is independent

        if useParallel

            parfor i = 1:N

                row = A(i,:);

                sigma = row*x-D(i)*x(i);

                xnew(i) = (b(i)-sigma)/D(i);

            end

        else

            for i = 1:N

                row = A(i,:);

                sigma = row*x-D(i)*x(i);

                xnew(i) = (b(i)-sigma)/D(i);

            end

        end

        % Update solution
        x = xnew;

        % Calculate residual
        r = b-parallelMatVec(A,x,useParallel);

        res = norm(r)/normb;

        if res <= tol
            break;
        end

    end

    elapsed = toc;

end


%% ============================================================
% FUNCTION 3: PARALLEL MATRIX-VECTOR MULTIPLICATION
%% ============================================================

function y = parallelMatVec(A,x,useParallel)

    N = length(x);

    y = zeros(N,1);

    if useParallel

        parfor i = 1:N

            y(i) = A(i,:)*x;

        end

    else

        for i = 1:N

            y(i) = A(i,:)*x;

        end

    end

end
