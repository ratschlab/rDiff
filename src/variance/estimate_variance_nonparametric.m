function [VARIANCE1, VARIANCE2]=estimate_variance_nonparametric(CFG,genes)

fprintf('Estimating variance function for rDiff.nonparametric\n')

VARIANCE1=[];
VARIANCE2=[];

%Getthe gene expression
fprintf('Loading gene expression\n')

if isempty(CFG.Counts_gene_expression)
    EXPR_TAB_FILENAME=[CFG.out_base 'Gene_expression.tab'];
else
    EXPR_TAB_FILENAME=CFG.Counts_gene_expression;
end
try
    Gene_expression=importdata(EXPR_TAB_FILENAME,'\t',1);
catch
    error(['Could not open: ' EXPR_TAB_FILENAME])
end
%C=importdata('../out/release_test/Gene_expression.tab','\t',1);
  

%Get the counts
fprintf('Loading alternative region counts\n')
if isempty(CFG.Counts_rDiff_nonparametric)
    IN_FILENAME=[CFG.out_base 'Nonparametric_region_counts.mat'];
    load(IN_FILENAME,'Counts_rDiff_nonparametric')
else
    IN_FILENAME=[CFG.out_base CFG.Counts_rDiff_nonparametric];
    load(IN_FILENAMEc,'Counts_rDiff_nonparametric') 
end



%Iterate over the functions to be generated
%compute means and variances

if CFG.compute_variance_function_1
    fprintf('estimating variance function for sample 1\n')
    %Get the samples to use for the for estimation the variance function
    SAMPLE_IX=find(CFG.SAMPLES==1);
    VARIANCE1=estimate_variance_helper(CFG,SAMPLE_IX,Counts_rDiff_nonparametric,Gene_expression.data);
    if not(isempty(CFG.save_variance_function_1))
        VARIANCE_FUNTION_OUTPATH=[CFG.out_base CFG.save_variance_function_1];
        save(VARIANCE_FUNTION_OUTPATH,'VARIANCE1');
    else
        %Use previously estimated variance function
        if not(isempty(CFG.variance_function_1))
            try
                VARIANCE_FUNTION_INPATH=[CFG.out_base CFG.variance_function_1];
                VARIANCE1=load(VARIANCE_FUNTION_INPATH);
            catch
                error(['Could not load variance function for sample 1 from: ' VARIANCE_FUNTION_INPATH])
            end
        end
    end
end

if CFG.compute_variance_function_2
    fprintf('estimating variance function for sample 2\n')
   %Get the samples to use for the for estimation the variance function
    SAMPLE_IX=find(CFG.SAMPLES==2);
    VARIANCE2=estimate_variance_helper(CFG,SAMPLE_IX,Counts_rDiff_nonparametric,Gene_expression.data);
    if not(isempty(CFG.save_variance_function_2))
        VARIANCE_FUNTION_OUTPATH=[CFG.out_base CFG.save_variance_function_2];
        save(VARIANCE_FUNTION_OUTPATH,'VARIANCE2');
    else
        %Use previously estimated variance function
        if not(isempty(CFG.variance_function_2))
            try
                VARIANCE_FUNTION_INPATH=[CFG.out_base CFG.variance_function_2];
                VARIANCE2=load(VARIANCE_FUNTION_INPATH);
            catch
                error(['Could not load variance function for sample 2 from: ' VARIANCE_FUNTION_INPATH])
            end
        end
    end
end



