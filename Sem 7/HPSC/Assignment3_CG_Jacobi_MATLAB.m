%% ASSIGNMENT 3 (PART ii): CG AND JACOBI USING MATLAB PARALLEL COMPUTING
% Solves the Assignment-2 structured symmetric systems:
%       A*x = b
% for matrix orders 100, 1000 and 10000, using 1, 2, 4 and 8 workers.
%
% IMPORTANT:
% MATLAB parfor / Parallel Computing Toolbox is used here. This is MATLAB
% parallel computation, NOT an OpenMP C implementation. If the assignment
% strictly requires OpenMP, MATLAB results should be described as a MATLAB
% parallel implementation rather than OpenMP.
%
% INPUT FILES expected in the Current Folder:
%   Kmat_100.txt       Fvec_100.txt
%   Kmat_1000.txt      Fvec_1000.txt
%   Kmat_10000.txt     Fvec_10000.txt
%
% Matrix-file interpretation follows the Assignment-2 MATLAB script:
%   data = load('Kmat_1000.txt');
%   N = sqrt(length(data));
%   A = reshape(data,N,N);
% MATLAB reshape is column-major. This is correct for a symmetric matrix.
%
% Run this entire file. Results are saved in Assignment3_results.csv.

clc; clear; close all;

%% USER SETTINGS
matrixOrders = [100, 1000, 10000];
workerCounts = [1, 2, 4, 8];

tol = 1e-10;
maxIter = 200000;

matrixPattern = 'Kmat_%d.txt';
rhsPattern    = 'Fvec_%d.txt';

results = table('Size',[0 8], ...
    'VariableTypes',{'double','string','double','double','double','double','double','string'}, ...
    'VariableNames',{'N','Solver','Workers','Iterations','Time_s', ...
                     'Residual','Speedup','Status'});

%% RUN EACH MATRIX SIZE
for s = 1:numel(matrixOrders)
    Nexpected = matrixOrders(s);

    matrixFile = sprintf(matrixPattern,Nexpected);
    rhsFile    = sprintf(rhsPattern,Nexpected);

    fprintf('\n====================================================\n');
    fprintf(' MATRIX ORDER: %d\n',Nexpected);
    fprintf('====================================================\n');

    if ~isfile(matrixFile) || ~isfile(rhsFile)
        fprintf('Skipping N=%d: missing %s or %s\n', ...
            Nexpected,matrixFile,rhsFile);
        continue;
    end

    % Load exactly as in the Assignment-2 MATLAB program.
    rawA = load(matrixFile);
    rawA = rawA(:);

    N = round(sqrt(numel(rawA)));
    if N ~= Nexpected || N*N ~= numel(rawA)
        error(['Matrix file %s has %d numeric entries; this is not ' ...
               'a %d-by-%d dense-entry file.'], ...
               matrixFile,numel(rawA),Nexpected,Nexpected);
    end

    A = reshape(rawA,N,N);
    clear rawA;

    b = load(rhsFile);
    b = b(:);

    if numel(b) ~= N
        error('RHS file %s has %d entries; expected %d.', ...
              rhsFile,numel(b),N);
    end

    % Convert to sparse storage after loading. This is important for FEM
    % matrices, which are normally sparse.
    A = sparse(A);

    % Basic matrix checks.
    symmetryError = norm(A-A','fro') / max(norm(A,'fro'),eps);
    d = full(diag(A));

    fprintf('Loaded %d-by-%d matrix; nnz(A) = %d\n',N,N,nnz(A));
    fprintf('Relative symmetry error = %.3e\n',symmetryError);

    if any(abs(d) < eps)
        error('Jacobi cannot be applied: at least one diagonal entry is zero.');
    end

    % CG requires a symmetric positive-definite matrix.
    % This check is diagnostic; chol may be costly for very large matrices.
    if symmetryError > 1e-10
        warning('A is not symmetric to the requested tolerance. CG may fail.');
    end

    % Store solver timings for speedup calculation.
    timeCG = nan(size(workerCounts));
    timeJ  = nan(size(workerCounts));

    for w = 1:numel(workerCounts)
        nw = workerCounts(w);

        % Start/reconfigure the local parallel pool BEFORE timing.
        pool = gcp('nocreate');
        if ~isempty(pool) && pool.NumWorkers ~= nw
            delete(pool);
            pool = [];
        end
        if isempty(pool)
            parpool('local',nw);
        end

        fprintf('\n--- %d MATLAB worker(s) ---\n',nw);

        % CG
        [xCG,itCG,tCG,resCG,flagCG] = cg_parallel(A,b,tol,maxIter);
        timeCG(w) = tCG;

        fprintf('CG     : iter = %6d | time = %.6f s | relres = %.3e | %s\n', ...
            itCG,tCG,resCG,flagCG);

        results = [results; {N,"CG",nw,itCG,tCG,resCG,NaN,string(flagCG)}]; %#ok<AGROW>

        % Jacobi
        [xJ,itJ,tJ,resJ,flagJ] = jacobi_parallel(A,b,tol,maxIter);
        timeJ(w) = tJ;

        fprintf('Jacobi : iter = %6d | time = %.6f s | relres = %.3e | %s\n', ...
            itJ,tJ,resJ,flagJ);

        results = [results; {N,"Jacobi",nw,itJ,tJ,resJ,NaN,string(flagJ)}]; %#ok<AGROW>

        % Keep solutions available for optional comparison.
        if nw == 1
            xCG_1 = xCG; %#ok<NASGU>
            xJ_1  = xJ;  %#ok<NASGU>
        end
    end

    % Calculate speedup relative to one worker for each solver.
    idxCG = find(results.N == N & results.Solver == "CG");
    idxJ  = find(results.N == N & results.Solver == "Jacobi");

    if ~isempty(idxCG) && isfinite(timeCG(1)) && timeCG(1) > 0
        results.Speedup(idxCG) = timeCG(1) ./ timeCG(:);
    end
    if ~isempty(idxJ) && isfinite(timeJ(1)) && timeJ(1) > 0
        results.Speedup(idxJ) = timeJ(1) ./ timeJ(:);
    end
end

%% DISPLAY AND SAVE
fprintf('\n\n================ FINAL RESULTS ================\n');
disp(results);
writetable(results,'Assignment3_results.csv');
fprintf('Saved results to Assignment3_results.csv\n');

% Plot speedup for each matrix order and solver.
validRows = isfinite(results.Speedup);
if any(validRows)
    figure('Name','Assignment 3 Speedup','Color','w');
    tiledlayout(1,2,'TileSpacing','compact');

    nexttile;
    hold on; grid on;
    for N = matrixOrders
        ids = results.N == N & results.Solver == "CG";
        if any(ids)
            [wp,ord] = sort(results.Workers(ids));
            sp = results.Speedup(ids);
            sp = sp(ord);
            plot(wp,sp,'-o','LineWidth',1.5, ...
                'DisplayName',sprintf('N=%d',N));
        end
    end
    plot(workerCounts,workerCounts,'k--','DisplayName','Ideal');
    xlabel('Number of workers'); ylabel('Speedup');
    title('Conjugate Gradient'); legend('Location','best');

    nexttile;
    hold on; grid on;
    for N = matrixOrders
        ids = results.N == N & results.Solver == "Jacobi";
        if any(ids)
            [wp,ord] = sort(results.Workers(ids));
            sp = results.Speedup(ids);
            sp = sp(ord);
            plot(wp,sp,'-o','LineWidth',1.5, ...
                'DisplayName',sprintf('N=%d',N));
        end
    end
    plot(workerCounts,workerCounts,'k--','DisplayName','Ideal');
    xlabel('Number of workers'); ylabel('Speedup');
    title('Jacobi'); legend('Location','best');
end

%% ======================== LOCAL FUNCTIONS =============================

function [x,iter,elapsed,relres,status] = cg_parallel(A,b,tol,maxIter)
% Conjugate Gradient with parallel row-wise sparse matrix-vector product.
% Vector dot products use MATLAB's sum operation.

    n = length(b);
    x = zeros(n,1);
    r = b - parallel_spmv(A,x);
    p = r;

    rr = sum(r.*r);
    normb = max(norm(b),eps);
    relres = sqrt(rr)/normb;

    status = 'maximum iterations reached';
    iter = 0;

    tStart = tic;
    if relres <= tol
        elapsed = toc(tStart);
        status = 'converged';
        return;
    end

    for k = 1:maxIter
        Ap = parallel_spmv(A,p);
        denom = sum(p.*Ap);

        if ~isfinite(denom) || denom <= 0
            status = 'CG breakdown: p''*A*p <= 0';
            iter = k-1;
            break;
        end

        alpha = rr/denom;
        x = x + alpha*p;
        r = r - alpha*Ap;

        rrNew = sum(r.*r);
        relres = sqrt(rrNew)/normb;
        iter = k;

        if relres <= tol
            status = 'converged';
            rr = rrNew;
            break;
        end

        beta = rrNew/rr;
        p = r + beta*p;
        rr = rrNew;
    end
    elapsed = toc(tStart);

    % Report the true residual, not only the recursively updated residual.
    rTrue = b - A*x;
    relres = norm(rTrue)/normb;
end

function [x,iter,elapsed,relres,status] = jacobi_parallel(A,b,tol,maxIter)
% Jacobi iteration with separate old/new vectors and parallel row updates.

    n = length(b);
    d = full(diag(A));
    R = A - spdiags(d,0,n,n);

    x = zeros(n,1);
    normb = max(norm(b),eps);
    relres = norm(b-A*x)/normb;
    status = 'maximum iterations reached';
    iter = 0;

    tStart = tic;
    if relres <= tol
        elapsed = toc(tStart);
        status = 'converged';
        return;
    end

    for k = 1:maxIter
        xNew = zeros(n,1);

        % Each worker computes distinct entries of xNew.
        parfor i = 1:n
            xNew(i) = (b(i) - R(i,:)*x)/d(i);
        end

        x = xNew;
        iter = k;

        % True residual for a consistent stopping test.
        relres = norm(b-A*x)/normb;
        if relres <= tol
            status = 'converged';
            break;
        end
    end
    elapsed = toc(tStart);
end

function y = parallel_spmv(A,x)
% Parallel sparse matrix-vector multiplication, distributed by rows.

    n = size(A,1);
    y = zeros(n,1);
    parfor i = 1:n
        y(i) = A(i,:)*x;
    end
end
