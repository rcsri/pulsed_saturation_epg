function [x_out, resnorm_out, resid, J] = run_lsqnonlin(F, x0_red, lb_red, ub_red, pars)
% F: objective function to be minimized
% x0_red: initial guess, reduced parameter set 
% lb_red: lower bound, reduced parameter set 
% ub_red: upper bound, reduced parameter set 
% pars: additional parameters 

% parallel processing
if verLessThan('matlab', '9')
    options = optimoptions('lsqnonlin');
else
    if isfield(pars, 'optimizer_use_parallel')
        optimizer_use_parallel = pars.optimizer_use_parallel;
    else
        optimizer_use_parallel = true;
    end

    options = optimoptions('lsqnonlin', 'UseParallel', optimizer_use_parallel);
end

% max function evaluations and iterations
if verLessThan('matlab', '9')
    options.MaxFunEvals = pars.max_iter;
    options.MaxIter = pars.max_func_eval; 
else
    options.MaxIterations = pars.max_iter;
    options.MaxFunctionEvaluations = pars.max_func_eval; 
end

% finite difference type
options.Display  = 'iter';
if verLessThan('matlab', '9')
    options.FinDiffType = 'central';
else
    options.FiniteDifferenceType = 'central';
end
disp('lsqnonlin options:')
disp(options); 

% successive minimizations
x = cell(1,pars.n_iter);
resnorm = cell(1,pars.n_iter);
[x{1}, resnorm{1}, resid, ~,~,~,J] = lsqnonlin( F , x0_red, lb_red, ub_red, options); % [x_scaled{1}, resnorm{1}, resid, exitflag,output,lambda,J] 
for iter = 2:pars.n_iter    
    [x{iter}, resnorm{iter}, resid,~,~,~,J] = lsqnonlin( F , x{iter-1}, lb_red, ub_red, options);
end

% output parameter and residuals for each minimization 
for iter = 1:pars.n_iter
    fprintf('%s, resnorm = %7.7f\n',num2str(x{iter}),resnorm{iter});
end

% output result in final iteration
x_out       = x{end}; 
resnorm_out = resnorm{end}; 
