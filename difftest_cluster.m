function PROBS=difftest_cluster()
% difftest_cluster(CFG)
%
% -- input --
% CFG: configuration struct


addpath('/fml/ag-raetsch/share/software/matlab_tools/rproc/')
addpath('/fml/ag-raetsch/share/svn/tools/utils');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/variance/final');
addpath /fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/experimental/
addpath /fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/variance/nbin/

%%%% paths

JB_NR=1;
TASK_NR=1;


REPLICATES=[1,2];

NR_OF_TASK=1;
 
CFG.IDX1=1;
CFG.IDX2=1;
CFG.paths = set_difftest_paths();
%%%% configuration
CFG = configure_difftest(CFG);
%%%% load genes
load(CFG.genes_path , 'genes');

[genes]=mask_dubl(genes,0);


JOB_INFO = rproc_empty();
RESULTS=cell(1,NR_OF_TASK);
%Prepare Variance functions
get_read_counts(CFG,genes);
JB_NR=1


CFG.IDX1=1;
CFG.IDX2=2;
CFG.paths = set_difftest_paths();
CFG.SUFF_NAME='NonAlt_rep';
%%%% configuration
CFG = configure_difftest(CFG);
%%%% load genes
nb_of_chunks=CFG.nb_of_chunks;
idx=[(1:size(genes,2))',ceil((1:size(genes,2))*nb_of_chunks/size(genes,2))']; 

FRAC_OF_OUTLRS=0;
COUNTER=1;
REPLICATES=[1 1];
REGION_IDEN=1;

BIN_LEN=50;
VARIANCE=get_gene_variance_biol_bin2(CFG,genes,REPLICATES,FRAC_OF_OUTLRS,BIN_LEN,REGION_IDEN); 
VARIANCES={VARIANCE};

CFG.VARIANCES=VARIANCES;
for i = 1:nb_of_chunks
  PAR.genes_idx = idx(idx(:,2)==i,1);
  %set parameters
  CFG.rproc_memreq = 10000;
  CFG.rproc_par.mem_req_resubmit = [16000 24000 32000];
  CFG.rproc_par.identifier = sprintf('and%i-',JB_NR);
  fprintf(1, 'Submitting job %i to cluster\n',JB_NR);
  CFG.JB_NR=JB_NR;
  PAR.CFG=CFG;
  JOB_INFO(JB_NR) = rproc('diff_expr_caller', PAR, CFG.rproc_memreq,CFG.rproc_par, CFG.rproc_time); %pause(1);
  JB_NR=JB_NR+1;
end

[JOB_INFO num_crashed] = rproc_wait(JOB_INFO, 60, 1, -1);
if num_crashed>0, pause(330); end; % some delay to wait until results are written
JB_NR=1;  


ERRORS_NR=[];
JB_NR=0;
idx=[(1:size(genes,2))',ceil((1:size(genes,2))*CFG.nb_of_chunks/ size(genes,2))'];
NR_OF_TASK=1;
RESULTS=cell(1,NR_OF_TASK) 
for i=1:NR_OF_TASK
  TEMP=cell(size(genes));
  RESULTS{i}=TEMP;
  for j=1:CFG.nb_of_chunks
    JB_NR=JB_NR+1;
    OUT_FILENAME=[CFG.out_base  CFG.SUFF_NAME '_' int2str(JB_NR)   '.mat'];
    try
      load(OUT_FILENAME)
      IDX=idx(idx(:,2)==j,1);
      for k=1:length(IDX)
	RESULTS{i}{IDX(k)}=STAT{k};
      end
    catch
      ERRORS_NR=[ERRORS_NR; JB_NR];
      IDX=idx(idx(:,2)==j,1);
      for k=1:length(IDX)
	RESULTS{i}{IDX(k)}={};
      end      
    end
    clear STAT;
  end
end


NR_OF_TASK=1;
S=size(genes,2);
NR_OF_TESTS=length(RESULTS{1}{1});
%Format Data



P_TABLE=ones(S,3);

for j=1:S
  NR_OF_TESTS=length(RESULTS{1}{j});
  T=RESULTS{1}{j};
  for i=3:5
    if not(isempty(T))
      if (isempty(T{i}))
	P_TABLE(j,i-2)=1;
      else	
	P_TABLE(j,i-2)=T{i};
      end
    end
  end
end
P_TABLE(P_TABLE>1)=1;


OUT_FILENAME=[CFG.out_base  'P_Values.tab'];

fid=fopen(OUT_FILENAME,'w');
fprintf(fid,'gene\tmmd p-value\tNB p-value\tpoisson p-value\n')

for i=1:S
  fprintf(fid,'%s\t%f\t%f\t%f\n',genes(i).name,P_TABLE(i,1),P_TABLE(i,2),P_TABLE(i,3));
end
fclose(fid);

fprintf('Saved the results under: %s\n',OUT_FILENAME)