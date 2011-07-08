function []=get_read_counts(CFG,genes)



addpath('/fml/ag-raetsch/share/software/matlab_tools/rproc/')
addpath('/fml/ag-raetsch/share/svn/tools/utils');
addpath('~/svn/tools/utils'); 
addpath('~/svn/tools/genomes');
addpath ~/svn/projects/RNASeq/difftest/variance
%%%% paths


 
JOB_INFO = rproc_empty();
JB_NR=1;
for RUN=1:size(CFG.BAM_FILES,2)
  CFG.paths = set_difftest_paths()
  %%%% configuration
  CFG.IDX1=RUN;
  CFG.IDX2=RUN;
  
  CFG = configure_difftest(CFG);
    
  %%%% paths
  
  
  nb_of_chunks=CFG.nb_of_chunks;
  
  idx=[(1:size(genes,2))',ceil((1:size(genes,2))*nb_of_chunks/size(genes,2))']; 
  
  % submit jobs to cluster
  for i = 1:nb_of_chunks
    PAR.genes = genes(idx(idx(:,2)==i,1));      
    %set parameters
    CFG.IDX=i;
    CFG.rproc_memreq = 5000;
    CFG.JB_NR=i; 
    CFG.rproc_par.mem_req_resubmit = [5000 10000 32000];
    %CFG.rproc_par.force_octave =1;
    CFG.rproc_par.identifier = sprintf('ar.f%i-',i);  
    CFG.outfile_prefix=[CFG.BAM_FILES{RUN} '_' num2str(i) '_of_' num2str(nb_of_chunks)];
    fprintf(1, 'Submitting job %i to cluster\n',i);
    PAR.CFG=CFG;
    JOB_INFO(JB_NR) = rproc('get_gene_variance_caller_single_alt2', PAR,CFG.rproc_memreq, CFG.rproc_par, CFG.rproc_time); %pause(1);
    JB_NR=JB_NR+1;
  end
  %results are written
  
end
[JOB_INFO num_crashed] = rproc_wait(JOB_INFO, 60, 1, -1);


