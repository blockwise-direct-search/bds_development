function [xopt, fopt, exitflag, output] = bds_norma(fun, x0, options)
%BDS solves unconstrained optimization problems without using derivatives by blockwise direct search methods.
%
%   BDS supports in MATLAB R2017b or later.
%
%   XOPT = BDS(FUN, X0) returns an approximate minimizer XOPT of the function FUN, starting the
%   calculations at X0. FUN must accept a vector input X and return a scalar.
%
%   XOPT = BDS(FUN, X0, OPTIONS) performs the computations with the options in OPTIONS. OPTIONS should be a
%   structure with the following fields.
%
%   Algorithm                           Algorithm to use. It can be "cbds" (cyclic blockwise direct search),
%                                       "pbds" (randomly permuted blockwise direct search), "rbds" (randomized
%                                       blockwise direct search), "ds" (the classical direct search without blocks).
%                                       "pads" (parallel blockwise direct search).
%                                       Default: "cbds".
%   nb                                  Number of blocks. A positive integer. Default: n if Algorithm is "cbds", "pbds",
%                                       or "rbds", 1 if Algorithm is "ds".
%   MaxFunctionEvaluations              Maximum of function evaluations. A positive integer. See also MaxFunctionEvaluations_factor.
%   MaxFunctionEvaluations_dim_factor   Factor to define the maximum number of function evaluations as a multiple
%                                       of the dimension of the problem. A positive integer. See also MaxFunctionEvaluations.
%                                       The maximum of function evaluations is min(MaxFunctionEvaluations, MaxFunctionEvaluations_factor*n) if the user
%                                       specify both MaxFunctionEvaluations and MaxFunctionEvaluations_factor; it is MaxFunctionEvaluations if the user only specifies
%                                       MaxFunctionEvaluations; it is MaxFunctionEvaluations_factor*n if the user only specifies MaxFunctionEvaluations_factor; it is
%                                       min(get_default_constant("MaxFunctionEvaluations"), get_default_constant("MaxFunctionEvaluations_factor")*n) if
%                                       the user specifies neither MaxFunctionEvaluations nor MaxFunctionEvaluations_factor.
%   direction_set                       A matrix whose columns will be used to define the polling directions.
%                                       If options does not contain direction_set, then the polling directions will be
%                                       {e_1, -e_1, ..., e_n, -e_n}. Otherwise, direction_set should be a matrix of n
%                                       rows, and the polling directions will be {d_1, -d_1, ..., d_m, -d_m}, where d_i
%                                       is the i-th column of direction_set, and m is the number of columns of direction_set.
%                                       If necessary, we will first extend direction_set by adding some columns to make
%                                       sure that rank(direction_set) = n, so that the polling directions make a
%                                       positive spanning set. See get_direction_set.m for details.
%   is_noisy                            A flag deciding whether the problem is noisy or not. Default: false.
%   expand                              Expanding factor of step size. A real number no less than 1.
%                                       It depends on the dimension of the problem and the algorithm and whether the problem is noisy.
%   shrink                              Shrinking factor of step size. A positive number less than 1.
%                                       It depends on the dimension of the problem and the algorithm and whether the problem is noisy.
%   forcing_function                    The forcing function used for deciding whether the step achieves a sufficient
%                                       decrease. A function handle. Default: @(alpha) alpha^2. See also reduction_factor.
%   reduction_factor                    Factors multiplied to the forcing function when deciding whether the step achieves
%                                       a sufficient decrease. A 3-dimentional vector such that
%                                       reduction_factor(1) <= reduction_factor(2) <= reduction_factor(3),
%                                       reduction_factor(1) >= 0, and reduction_factor(2) > 0.
%                                       reduction_factor(0) is used for deciding whether to update the base point;
%                                       reduction_factor(1) is used for deciding whether to shrink the step size;
%                                       reduction_factor(2) is used for deciding whether to expand the step size.
%                                       Default: [0, eps, eps]. See also forcing_function.
%   StepTolerance                       Lower bound of the step size. If the step size is smaller than StepTolerance,
%                                       then the algorithm terminates. A (small) positive number. Default: 1e-10.
%   ftarget                             Target of the function value. If the function value is smaller than or equal to
%                                       ftarget, then the algorithm terminates. A real number. Default: -Inf.
%   polling_inner                       Polling strategy in each block. It can be "complete" or "opportunistic".
%                                       Default: "opportunistic".
%   cycling_inner                       Cycling strategy employed within each block. It is used only when polling_inner
%                                       is "opportunistic". It can be 0, 1, 2, 3, 4. See cycling.m for details.
%                                       Default: 3.
%   permuting_period                    It is only used in PBDS, which shuffles the blocks every permuting_period
%                                       iterations. A positive integer. Default: 1.
%   replacement_delay                   It is only used for RBDS. Suppose that replacement_delay is r. If block i
%                                       is selected at iteration k, then it will not be selected at iterations
%                                       k+1, ..., k+r. An integer between 0 and nb-1. Default: 0.
%   seed                                The seed for permuting blocks in PBDS or randomly choosing one block in RBDS.
%                                       It is only for reproducibility in experiments. A positive integer.
%   output_xhist                        Whether to output the history of points visited. Default: false.
%   output_alpha_hist                   Whether to output the history of step sizes. Default: false.
%   output_block_hist                   Whether to output the history of blocks visited. Default: false.
%
%   [XOPT, FOPT] = BDS(...) returns an approximate minimizer XOPT and its function value FOPT.
%
%   [XOPT, FOPT, EXITFLAG] = BDS(...) also returns an EXITFLAG that indicates the exit
%   condition. The possible values of EXITFLAG are 0, 1, 2, 3.
%
%   0    The StepTolerance of the step size is reached.
%   1    The target of the objective function is reached.
%   2    The maximum number of function evaluations is reached.
%   3    The maximum number of iterations is reached.
%
%   [XOPT, FOPT, EXITFLAG, OUTPUT] = BDS(...) returns a
%   structure OUTPUT with the following fields.
%
%   fhist        History of function values.
%   xhist        History of points visited (if output_xhist is true).
%   alpha_hist   History of step size for every iteration (if alpha_hist is true).
%   blocks_hist  History of blocks visited (if block_hist is true).
%   funcCount    The number of function evaluations.
%   message      The information of EXITFLAG.
%
%   ***********************************************************************
%   Authors:    Haitian LI (hai-tian.li@connect.polyu.hk)
%               and Zaikun ZHANG (zaikun.zhang@polyu.edu.hk)
%               Department of Applied Mathematics,
%               The Hong Kong Polytechnic University
%   ***********************************************************************
%   All rights reserved.
%

% Set options to an empty structure if it is not provided.
if nargin < 3
    options = struct();
end

% Transpose x0 if it is a row.
x0_is_row = isrow(x0);
x0 = double(x0(:));

% Set the default value of debug_flag. If options do not contain debug_flag,
% then debug_flag is false.
if isfield(options, "debug_flag")
    debug_flag = options.debug_flag;
else
    debug_flag = false;
end

% Check the inputs of the user when debug_flag is true.
if debug_flag
    verify_preconditions(fun, x0, options);
end

% If FUN is a string, then convert it to a function handle.
if ischarstr(fun)
    fun = str2func(fun);
end
% Redefine fun to accept columns if x0 is a row, as we use columns internally.
if x0_is_row
    fun = @(x)fun(x');
end

% Get the dimension of the problem.
n = length(x0);
% Set the default Algorithm of BDS, which is "cbds".
if ~isfield(options, "Algorithm")
    options.Algorithm = get_default_constant("Algorithm");
end

% Get the direction set.
D = get_direction_set(n, options);

% Get the number of directions.
m = size(D, 2);

% Get the number of blocks.
if isfield(options, "nb")
    % The number of directions should be greater or equal to the number of blocks.
    nb = min(m, options.nb);
elseif strcmpi(options.Algorithm, "cbds") || strcmpi(options.Algorithm, "pbds") ...
        || strcmpi(options.Algorithm, "rbds") || strcmpi(options.Algorithm, "pads") ...
        || strcmpi(options.Algorithm, "sCBDS")
    % Default value is set as n, which is good for canonical with 2n directions. For
    % other situations, other value may be good.
    nb = n;
elseif strcmpi(options.Algorithm, "ds")
    nb = 1;
end

% Set indices of blocks as 1:nb.
all_block_indices = 1:nb;

% Set the default value of noisy.
if ~isfield(options, "is_noisy")
    options.noisy = get_default_constant("is_noisy");
end

% Set the value of expand and shrink according to the dimension of the problem
% and whether the problem is noisy or not, also according to the Algorithm.
% n == 1 is treated as a special case, and we consider the Algorithm to be "ds".
if strcmpi(options.Algorithm, "ds") || n == 1
    if numel(x0) <= 5
        expand = get_default_constant("ds_expand_small");
        shrink = get_default_constant("ds_shrink_small");
    else
        % Judge whether the problem is noisy or not.
        if isfield(options, "is_noisy") && options.is_noisy
            expand = get_default_constant("ds_expand_big_noisy");
            shrink = get_default_constant("ds_shrink_big_noisy");
        else
            expand = get_default_constant("ds_expand_big");
            shrink = get_default_constant("ds_shrink_big");
        end
    end
else
    if numel(x0) <= 5
        expand = get_default_constant("expand_small");
        shrink = get_default_constant("shrink_small");
    else
        % Judge whether the problem is noisy or not.
        if isfield(options, "is_noisy") && options.is_noisy
            expand = get_default_constant("expand_big_noisy");
            shrink = get_default_constant("shrink_big_noisy");
        else
            expand = get_default_constant("expand_big");
            shrink = get_default_constant("shrink_big");
        end
    end
end

% Set the value of expand if options contains expand.
if isfield(options, "expand")
    expand = options.expand;
end

% Set the value of shrink if options contains shrink.
if isfield(options, "shrink")
    shrink = options.shrink;
end

% Set the maximum number of function evaluations. If the options do not contain MaxFunctionEvaluations,
% it is set to MaxFunctionEvaluations_dim_factor*n, where n is the dimension of the problem.
if isfield(options, "MaxFunctionEvaluations")
    MaxFunctionEvaluations = options.MaxFunctionEvaluations;
else
    MaxFunctionEvaluations = get_default_constant("MaxFunctionEvaluations_dim_factor")*n;
end

% Each iteration will at least use one function evaluation. We will perform at most MaxFunctionEvaluations iterations.
% We set maxit in the following way to avoid that maxit is reached but MaxFunctionEvaluations is not exhausted when the
% algorithm terminates.
maxit = MaxFunctionEvaluations;

% Set the value of reduction factor.
if isfield(options, "reduction_factor")
    reduction_factor = options.reduction_factor;
else
    reduction_factor = get_default_constant("reduction_factor");
end

% Set the forcing function, which should be the function handle.
if isfield(options, "forcing_function")
    forcing_function = options.forcing_function;
else
    forcing_function = get_default_constant("forcing_function");
end
% If the forcing function is a string, then convert it to a function handle.
if isfield(options, "forcing_function_type")
    switch options.forcing_function_type
        case "quadratic"
            forcing_function = @(x)x.^2;
        case "cubic"
            forcing_function = @(x)x.^3;
    end
end

% Set the value of StepTolerance. The algorithm will terminate if the stepsize is less than
% the StepTolerance.
if isfield(options, "StepTolerance")
    alpha_tol = options.StepTolerance;
else
    alpha_tol = get_default_constant("StepTolerance");
end

% Set the target of the objective function.
if isfield(options, "ftarget")
    ftarget = options.ftarget;
else
    ftarget = get_default_constant("ftarget");
end

% Set the value of polling_inner. This is the polling strategy employed within one block.
if ~isfield(options, "polling_inner")
    options.polling_inner = get_default_constant("polling_inner");
end

% Set the value of cycling_inner, which represents the cycling strategy inside each block.
if isfield(options, "cycling_inner")
    cycling_inner = options.cycling_inner;
else
    cycling_inner = get_default_constant("cycling_inner");
end

% Set the value of permuting_period, which permutes the blocks every permuting_period iterations.
if strcmpi(options.Algorithm, "pbds")
    if isfield(options, "permuting_period")
        permuting_period = options.permuting_period;
    else
        permuting_period = get_default_constant("permuting_period");
    end
end

% Set replacement_delay and num_selected_blocks. This is done only when
% Algorithm is "rbds", which randomly selects num_selected_blocks blocks in each
% iteration. If replacement_delay is r, then the block that is selected in the
% current iteration will not be selected in the next r iterations. Note that
% replacement_delay cannot exceed floor(nb/num_selected_blocks)-1.
% While a larger replacement_delay can potentially improve performance, 
% we set the default value to 0 to maintain the simplicity and consistency of the algorithm.
if strcmpi(options.Algorithm, "rbds")
    if isfield(options, "num_selected_blocks")
        num_selected_blocks = min(options.num_selected_blocks, nb);
    else
        num_selected_blocks = 1;
    end

    if isfield(options, "replacement_delay")
        replacement_delay = min(options.replacement_delay, floor(nb/num_selected_blocks)-1);
    else
        replacement_delay = 0;
    end

    % fprintf("bds_norma.m: expand = %f, shrink = %f, replacement_delay = %d\n", expand, shrink, replacement_delay);

end

% Initialize the step sizes and alpha_hist, which is the history of step sizes.
if isfield(options, "output_alpha_hist")
    output_alpha_hist = options.output_alpha_hist;
else
    output_alpha_hist = get_default_constant("output_alpha_hist");
end
% If alpha_hist exceeds the maximum of memory size limit, then we will not output alpha_hist.
if output_alpha_hist
    try
        alpha_hist = NaN(nb, maxit);
    catch
        output_alpha_hist = false;
        warning("The size of alpha_hist exceeds the maximum of memory size limit.")
    end
end

if isfield(options, "alpha_init")
    if length(options.alpha_init) == 1
        alpha_all = options.alpha_init*ones(nb, 1);
    elseif length(options.alpha_init) == nb
        alpha_all = options.alpha_init;
    else
        error("The length of alpha_init should be equal to nb or equal to 1.");
    end
    % Try alpha_all = 0.5 * max(abs(x0), 1) in the canonical case.
elseif isfield(options, "alpha_init_scaling") && options.alpha_init_scaling
    %alpha_all = 0.1 * ones(nb, 1);
    %alpha_all(x0 ~= 0) = 0.1 * abs(x0(x0 ~= 0));
    alpha_all = 0.5 * max(abs(x0), ones(nb, 1));
    %alpha_all = 0.1 * max(1e-3, abs(x0));
else
    alpha_all = ones(nb, 1);
end

% Determine the indices of directions in each block.
direction_set_indices = divide_direction_set(m, nb);

% Initialize the history of function values.
fhist = NaN(1, MaxFunctionEvaluations);

% Initialize the history of points visited.
if isfield(options, "output_xhist")
    output_xhist = options.output_xhist;
else
    output_xhist = get_default_constant("output_xhist");
end
% If xhist exceeds the maximum of memory size limit, then we will not output xhist.
if output_xhist
    try
        xhist = NaN(n, MaxFunctionEvaluations);
    catch
        output_xhist = false;
        warning("xhist will be not included in the output due to the limit of memory.");
    end
end

if isfield(options, "output_block_hist")
    output_block_hist = options.output_block_hist;
else
    output_block_hist = get_default_constant("output_block_hist");
end
% Initialize the history of blocks visited.
block_hist = NaN(1, MaxFunctionEvaluations);

% To avoid that the users bring some randomized strings.
if ~isfield(options, "seed")
    options.seed = get_default_constant("seed");
end
random_stream = RandStream("mt19937ar", "Seed", options.seed);

% Initialize the exitflag where the maximum number of iterations is reached.
exitflag = get_exitflag("MAXIT_REACHED");
xbase = x0;
[fbase, fbase_real] = eval_fun(fun, xbase);
% Set the number of function evaluations.
nf = 1;
if output_xhist
    xhist(:, nf) = xbase;
end
fhist(nf) = fbase_real;
xopt = xbase;
fopt = fbase;
terminate = false;

if nf >= MaxFunctionEvaluations || fbase_real <= ftarget
    % Either MaxFunctionEvaluations has been reached at the very first function evaluation
    % or FTARGET has been reached at the very first function evaluation.
    % In this case, no further computation should be entertained, and hence,
    % no iteration should be run.
    maxit = 0;
end
if fbase_real <= ftarget
    exitflag = get_exitflag( "FTARGET_REACHED");
elseif nf >= MaxFunctionEvaluations
    exitflag = get_exitflag("MAXFUN_REACHED");
end

% fopt_all(i) records the best function values encountered in the i-th block after one iteration,
% and xopt_all(:, i) is the corresponding value of x.
fopt_all = NaN(1, nb);
xopt_all = NaN(n, nb);

for iter = 1:maxit

    if strcmpi(options.Algorithm, "ds") || strcmpi(options.Algorithm, "cbds")
        % If the Algorithm is "ds", "cbds", then we will visit all blocks in order.
        % When the Algorithm is "ds", note that num_blocks = 1 and block_indices = [1],
        % a vector of length 1.
        block_indices = all_block_indices;
    end

    if strcmpi(options.Algorithm, "pads")
        % If the Algorithm is "pads", then we will visit all blocks in order.
        % However, we will permute the blocks every iteration.
        block_indices = random_stream.randperm(nb, nb);
    end

    % Permute the blocks every permuting_period iterations if the Algorithm is "pbds".
    % Why iter-1? Since we will permute block_indices at the initial stage.
    if strcmpi(options.Algorithm, "pbds") && mod(iter - 1, permuting_period) == 0
        % Ensure that permuting_period is properly defined when the selected Algorithm is "pbds".
        % A random permutation of 1:nb is generated using randperm(nb, nb), as the latest version of 
        % bds.m employs randperm(nb, nb) instead of randperm(nb).
        block_indices = random_stream.randperm(nb, nb);
    end
    
    % Get the block that is going to be visited if the Algorithm is "rbds".
    if strcmpi(options.Algorithm, "rbds")
        % Get the blocks that are going to be visited in this iteration when the Algorithm is "rbds".
        % These blocks should not have been visited in the previous replacement_delay iterations.
        % Note that block_indices is a vector of length num_selected_blocks.
        unavailable_block_indices = unique(block_hist(max(1, (iter-replacement_delay) * num_selected_blocks) : (iter-1) * num_selected_blocks), 'stable');
        available_block_indices = setdiff(all_block_indices, unavailable_block_indices);
        % Select num_selected_blocks blocks randomly from the available blocks. randsample
        % may be more efficient than randperm when num_selected_blocks is much smaller
        % than num_blocks. However, randsample is only available in Statistics and Machine
        % Learning Toolbox.
        block_indices = available_block_indices(random_stream.randperm(length(available_block_indices), num_selected_blocks));
    end

    if strcmpi(options.Algorithm, "sCBDS")
        block_indices = [1:nb (nb-1):-1:2];
    end

    for i = 1:length(block_indices)

        % If block_indices is 1 3 2, then block_indices(2) = 3, which is the real block that we are
        % going to visit.
        i_real = block_indices(i);

        % Record the number of blocks visited.
        num_visited = sum(~isnan(block_hist));
        % Record the block that is going to be visited.
        block_hist(num_visited+1) = i_real;

        % Get indices of directions in the i-th block.
        direction_indices = direction_set_indices{i_real};

        suboptions.MaxFunctionEvaluations = MaxFunctionEvaluations - nf;
        suboptions.cycling_inner = cycling_inner;
        suboptions.reduction_factor = reduction_factor;
        suboptions.forcing_function = forcing_function;
        suboptions.ftarget = ftarget;
        suboptions.polling_inner = options.polling_inner;

        [sub_xopt, sub_fopt, sub_exitflag, sub_output] = inner_direct_search(fun, xbase,...
            fbase, D(:, direction_indices), direction_indices,...
            alpha_all(i_real), suboptions);

        % Record the step size when the algorithm executes inner_direct_search above.
        if output_alpha_hist
            alpha_hist(:, iter) = alpha_all;
        end

        % Store the history of the evaluations by inner_direct_search,
        % and accumulate the number of function evaluations.
        fhist((nf+1):(nf+sub_output.nf)) = sub_output.fhist;
        % Accumulate the points visited by inner_direct_search.
        if output_xhist
            xhist(:, (nf+1):(nf+sub_output.nf)) = sub_output.xhist;
        end
        nf = nf+sub_output.nf;

        % Record the best function value and point encountered in the i_real-th block.
        fopt_all(i_real) = sub_fopt;
        xopt_all(:, i_real) = sub_xopt;

        % Retrieve the indices of the i_real-th block in the direction set.
        direction_set_indices{i_real} = sub_output.direction_indices;

        % Whether to update xbase and fbase. xbase serves as the "base point" for the computation in the next block,
        % meaning that reduction will be calculated with respect to xbase, as shown above.
        % Note that their update requires a sufficient decrease if reduction_factor(1) > 0.
        update_base = (reduction_factor(1) <= 0 && sub_fopt < fbase) ...
                    || (sub_fopt + reduction_factor(1) * forcing_function(alpha_all(i_real)) < fbase);

        % Update the step size alpha_all according to the reduction achieved.
        if sub_fopt + reduction_factor(3) * forcing_function(alpha_all(i_real)) < fbase
            alpha_all(i_real) = expand * alpha_all(i_real);
        elseif sub_fopt + reduction_factor(2) * forcing_function(alpha_all(i_real)) >= fbase
            % alpha_all(i_real) = max(shrink * alpha_all(i_real), alpha_threshold);
            alpha_all(i_real) = shrink * alpha_all(i_real);
        end

        % If the scheme is not "parallel", then we will update xbase and fbase after finishing the
        % direct search in the i_real-th block. For "parallel", we will update xbase and fbase after
        % one iteration of the outer loop.
        if ~strcmpi(options.Algorithm, "pads")
            if update_base
                xbase = sub_xopt;
                fbase = sub_fopt;
            end
        end

        % If sub_output.terminate is true, then inner_direct_search returns
        % boolean value of terminate because either the maximum number of function
        % evaluations or the target of the objective function value is reached.
        % In both cases, the exitflag is set by inner_direct_search.
        if sub_output.terminate
            terminate = true;
            exitflag = sub_exitflag;
            break;
        end

        % Terminate the computations if the largest component of step size is below a
        % given StepTolerance.
        if max(alpha_all) < alpha_tol
            terminate = true;
            exitflag = get_exitflag("SMALL_ALPHA");
            break;
        end
    end

    % Update xopt and fopt. Note that we do this only if the iteration encounters a strictly better point.
    % Make sure that fopt is always the minimum of fhist after the moment we update fopt.
    % The determination between fopt_all and fopt is to avoid the case that fopt_all is
    % bigger than fopt due to the update of xbase and fbase.
    [~, index] = min(fopt_all, [], "omitnan");
    if fopt_all(index) < fopt
        fopt = fopt_all(index);
        xopt = xopt_all(:, index);
    end

    % If the algorithm is "pads", then we only update xbase and fbase before the calculation or after, which
    % implies that in one iteration, each block will use the same xbase and fbase.
    if strcmpi(options.Algorithm, "pads")
        if (reduction_factor(1) <= 0 && fopt < fbase) || fopt + reduction_factor(1) * forcing_function(min(alpha_all)) < fbase
            xbase = xopt;
            fbase = fopt;
        end
    end

    % Check whether one of SMALL_ALPHA, MAXFUN_REACHED, and FTARGET_REACHED is reached.
    if terminate
        break;
    end

end

% Truncate HISTORY into a vector of nf length.
output.funcCount = nf;
output.fhist = fhist(1:nf);

if output_xhist
    output.xhist = xhist(:, 1:nf);
end

if output_alpha_hist
    output.alpha_hist = alpha_hist(1:min(iter, maxit));
end

num_blocks_visited = sum(~isnan(block_hist));
if output_block_hist
    output.blocks_hist = block_hist(1:num_blocks_visited);
end

switch exitflag
    case {get_exitflag("SMALL_ALPHA")}
        output.message = "The StepTolerance of the step size is reached.";
    case {get_exitflag("MAXFUN_REACHED")}
        output.message = "The maximum number of function evaluations is reached.";
    case {get_exitflag("FTARGET_REACHED")}
        output.message = "The target of the objective function is reached.";
    case {get_exitflag("MAXIT_REACHED")}
        output.message = "The maximum number of iterations is reached.";
    otherwise
        output.message = "Unknown exitflag";
end

% Transpose xopt if x0 is a row.
if x0_is_row
    xopt = xopt';
end

% verify_postconditions is to detect whether the output is in the right form when debug_flag is true.
if debug_flag
    verify_postconditions(fun, xopt, fopt, exitflag, output);
end
