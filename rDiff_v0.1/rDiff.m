
function []=rDiff()
% rDiff()
%

%%% Add paths  %%%
fprintf('Set the paths\n')
CFG.paths = set_rDiff_paths();


%%% Read configuration file %%%
fprintf('Load configuration\n')
CFG = configure_rDiff(CFG);
CFG = process_configure_rDiff(CFG);

%%% Get read counts %%%

%load the gene structure
load(CFG.genes_path , 'genes');

tic
% mask the regions which overlap with other genes
fprintf('Compute regions common to multiple genes\n')
[genes]=detect_overlapping_regions(genes);
toc;tic
%Precompute testing regions
fprintf('Compute alternative regions\n')
[genes]=compute_testing_region(CFG,genes);
toc;tic
%Get the gene expression
if CFG.estimate_gene_expression
    fprintf('Measure gene expression\n')
    get_read_counts(CFG,genes);
end
toc;tic
%%% Estimate variance function %%%
if CFG.perform_nonparametric
    variance_function_nonparametric_1=[];
    variance_function_nonparametric_2=[];
    [variance_function_nonparametric_1, variance_function_nonparametric_2]=estimate_variance_nonparametric(CFG,genes);
end
toc;tic

if CFG.perform_parametric
    variance_function_parametric_1=[];
    variance_function_parametric_2=[];
    [variance_function_parametric_1, variance_function_parametric_2]=estimate_variance_parametric(CFG,genes); 
end
toc;tic
%%% Perform tests &  Write output %%%
toc;tic
%Run the prametric tests
if or(CFG.perform_parametric,CFG.perform_poisson)
    perform_parametric_tests(CFG,genes,variance_function_parametric_1, variance_function_parametric_2)
end
toc;tic
%Run the nonparametric tests
if or(CFG.perform_nonparametric,CFG.perform_mmd)
    perform_nonparametric_tests(CFG,genes,variance_function_nonparametric_1, variance_function_nonparametric_2)
end
toc
return



