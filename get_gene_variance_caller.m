function []=get_gene_variance_caller(PAR)

samtools='/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/samtools/';
CFG = PAR.CFG;

genes = PAR.genes;
%BAM_FILES=CFG.BAM_FILES;
JB_NR=CFG.JB_NR;

%load ../POS.mat

clear PAR;

%%%% paths
%addpath(CFG.paths);
addpath('~/svn/projects/mGene_core/experiments/thaliana/');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/experiments/');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/read_utils/');
out_base=CFG.out_base;
%sim_data_base=CFG.sim_data_base;
real_data_base=CFG.real_data_base;

expr2_bam=CFG.expr2_bam;
expr1_bam=CFG.expr1_bam;

%genes=PAR.genes;
STAT=cell(size(genes));
do_save = 0;
out_base = '/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/out';
save_dir=fullfile(out_base,'cache');


for i=1:size(genes,2)
  RESULT=cell(1,7);
  gene = genes(i);
  
  RESULT{4}=JB_NR;
  RESULT{1}=gene.name;
  load_only = false;
  if isempty(gene.exons)
    PV=1;
    RESULT{2}=[inf,inf];
    RESULT{3}=[inf,inf];
    RESULT{5}=inf;
    STAT{i}=RESULT;
    continue;
  end

  gene.name 
  %  if or(strcmp(expr1_bam,'/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/reads/Col_0_nss_R1.bam'),strcmp(expr1_bam,'/fml/ag-raetsch/nobackup/projects/sequencing_runs/A_thaliana_magic/reads/Col_0_nss_R2.bam'))
    [reads1,reads2] = get_reads_dual2(expr1_bam,expr2_bam,gene,samtools,false,10);

    %else
    %[reads1,reads2] = get_reads_duals_2_mod_both(expr1_bam,expr2_bam,gene,samtools,false,10);
    %end
 
  
   
  
  
  RESULT{2}=[size(reads1,1),size(reads2,1)]
  reads1=reads1(sum(reads1,2)>30,:);
  reads2=reads2(sum(reads2,2)>30,:);
  
  
  NR_OF_TRANS=size(gene.transcripts,2);
  [SPLICINGEVENTS,SEQUENCE,EXONSEQUENCE]=splicingsequence(gene);
  EXON_IDX=sum(EXONSEQUENCE,1)==NR_OF_TRANS;
  
  %Get starts of the reads
  [TEMP,START1]=max(reads1,[],2);
  [TEMP,START2]=max(reads2,[],2);
  S1=size(reads1,1);
  new_reads1=sparse((1:S1)',START1,ones(S1,1),S1,size(reads1,2),S1);
  S2=size(reads2,1);
  new_reads2=sparse((1:S2)',START2,ones(S2,1),S2,size(reads2,2),S2);
   
  
  %Count the reads that start in the non-diferential exons
  IDX_VECT=zeros(1,size(reads2,2));
  for j=find(EXON_IDX)
    IDX_VECT(SPLICINGEVENTS(j):(SPLICINGEVENTS(j+1)-1))=1;
  end
  
  S1=new_reads1(:,IDX_VECT>0);
  S2=new_reads2(:,IDX_VECT>0);

  
  
  POSITIONS=sum(IDX_VECT);
  
  COUNT1=new_reads1*IDX_VECT';
  COUNT2=new_reads2*IDX_VECT';
  RESULT{3}=[sum(COUNT1),sum(COUNT2,1)]/POSITIONS
  RESULT{5}=POSITIONS
  RESULT{6}=sum(S1,1);
  RESULT{7}=sum(S2,1);
  
  [size(reads1,1),size(reads2,1)]
  % old and weighted poisson new ,weighted regions reads and
  % unexplained reads

  
  
  clear reads1;
  clear reads2;
  STAT{i}=RESULT;
  
end;

OUT_FILENAME=[CFG.out_base 'count_temp_' int2str(JB_NR)   '.mat'];
    save(OUT_FILENAME,'STAT')
