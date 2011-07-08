function CFG = configure_difftest(CFG)
% configure_rquant(CFG)
%
% -- input --
% CFG.organism
% CFG.exp
% CFG.read_len
% CFG.gene_source
% CFG.read_maps_select
% CFG.fname
% CFG.paths






%%%%% directories from which to load read data and genes %%%%%


% This is the directory in which the results and temporary results
% are saved.
CFG.out_base = '../out/analysis_artificial_variance_03_06_2011/';

%hits is the directory where the Bam-Files are
CFG.real_data_base = '/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/data_sim/simulation_2011-04-12/';

%These are the Bam-filenames
CFG.BAM_FILES={'genes_expr_1_1_weak_bias01.bam','genes_expr_1_2_weak_bias01.bam','genes_expr_2_1_weak_bias01.bam','genes_expr_2_2_weak_bias01.bam',};

BAM_FILES=CFG.BAM_FILES;
IDX1=CFG.IDX1;
IDX2=CFG.IDX2;
CFG.expr1_bam = fullfile(CFG.real_data_base,BAM_FILES{IDX1});
CFG.expr2_bam = fullfile(CFG.real_data_base,BAM_FILES{IDX2});

% this is the path where the genes-structure is saved
CFG.genes_path='/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/data_sim/simulation_2011-04-12/genes.mat';










%%%%% paramaters
CFG.nb_of_chunks=400;
CFG.SEQUENCED_LENGTH=75;	
		
%%%%% rproc settings %%%%%
CFG.use_rproc = 1; % 1: cluster submission or 0: locally
if CFG.use_rproc,
  CFG.rproc_num_jobs              = 400;
  CFG.rproc_memreq                = 4000;
  CFG.rproc_par.priority          = 105;
  CFG.rproc_par.resubmit          = 3;
  CFG.rproc_par.mem_req_resubmit  = [8000 12000 30000];
  CFG.rproc_par.time_req_resubmit = [60*60 100*60 90*60];
  CFG.rproc_par.express           = 0;
  CFG.rproc_par.immediately_bg    = 0;
  CFG.rproc_par.immediately       = 0;
  CFG.rproc_par.arch              = 64;
  CFG.rproc_par.identifier        = '';
  CFG.rproc_par.verbosity         = 0;
  CFG.rproc_time                  = 5*60; % mins
end
		
