function STAT = diff_expr_caller(PAR)
% genes = opt_transcripts_caller(PAR)
%
% -- input --
% PAR contains
%     CFG: configuration struct
%     genes: struct defining genes with start, stops, exons etc.
%     profile_weights: weights of profile functions
%     intron_dists: distances to closest intron
%
% -- output --
% genes: struct with additional fields of eg. estimated transcript weights
addpath ../tests
addpath ../variance/final
addpath ../experimental

samtools='/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/samtools/';
CFG = PAR.CFG;
SUFF_NAME=CFG.SUFF_NAME; 
load(CFG.genes_path , 'genes');
[genes]=mask_dubl(genes,0);

VARIANCES=CFG.VARIANCES;
genes = genes(PAR.genes_idx);



BAM_FILES=CFG.BAM_FILES;
JB_NR=CFG.JB_NR;


clear PAR;

%%%% paths
addpath(CFG.paths);
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/variance/nbin');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/experiments/');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/tests/');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/read_utils/');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/tests/sequence_tools/');
addpath(genpath('/fml/ag-raetsch/home/drewe/svn/tools/chronux/'))
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/tests/functions/');

out_base=CFG.out_base;
real_data_base=CFG.real_data_base;

expr2_bam=CFG.expr2_bam;
expr1_bam=CFG.expr1_bam;

CFG.SEQUENCED_LENGTH=80;
SEQUENCED_LENGTH=80;
%make map: gene->file

IDX1=CFG.IDX1;
IDX2=CFG.IDX2;

BAM_FILES{IDX1}
BAM_FILES{IDX2}


%genes=PAR.genes;
STAT=cell(size(genes));
do_save = 0;
out_base = '/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/out';
save_dir=fullfile(out_base,'cache');

for i=1:size(genes,2)
  %try
  RESULT=cell(1,7)
  gene = genes(i);
  
  RESULT{7}=JB_NR;
  RESULT{1}=gene.name;
  load_only = false;
   
  if or(isempty(gene.exons),gene.stop-gene.start<=SEQUENCED_LENGTH)
     continue;
  end
  
  PV1=1;
  for l=1:size(gene.transcripts,2)
    EXONS=gene.exons{l}
    if sum(EXONS(:,2)-EXONS(:,1))<=SEQUENCED_LENGTH	  
	  continue;
    end
  end
  if PV1==2
    continue;
  end
  gene.name 
    
   
  
  %adjust start and stop positions on the negative strand
  if strcmp(gene.strand,'-') 
    gene.start=gene.start+1;
    gene.stop=gene.stop+1;
    [reads1,reads2] = get_reads_dual2(expr1_bam,expr2_bam,gene,samtools,false,10);
    [reads1,FLAG]=remove_reads_from_other_genes(reads1,gene);
    [reads2,FLAG]=remove_reads_from_other_genes(reads2,gene);
    gene.start=gene.start-1;
    gene.stop=gene.stop-1;
  else
    [reads1,reads2] = get_reads_dual2(expr1_bam,expr2_bam,gene,samtools,false,10);
    [reads1,FLAG]=remove_reads_from_other_genes(reads1,gene);
    [reads2,FLAG]=remove_reads_from_other_genes(reads2,gene);
  end
  
  reads1=reads1(sum(reads1,2)>70,:);
  reads2=reads2(sum(reads2,2)>70,:);
  %% TEST with read cliping
  

  
  
  
  
  CLIP=0;
  if 1==0
  if (size(reads1,1)>0)
    for i=1:size(reads1,1)
      reads1(i,find(reads1(i,:),2,'first'))=0;
    end
    for i=1:size(reads1,1)
      reads1(i,find(reads1(i,:),2,'last'))=0;
    end
  end
  if (size(reads2,1)>0)
    for i=1:size(reads2,1)
      reads2(i,find(reads2(i,:),2,'first'))=0;
    end    
    for i=1:size(reads2,1)
      reads2(i,find(reads2(i,:),2,'last'))=0;
    end
  end
  CLIP=4;
  end
  
  
  %% END TEST
  
  [size(reads1,1),size(reads2,1)]
  
  FLAG=0;
  for k=1:size(gene.transcripts,2)
    if sum(gene.exons{k}(:,2)-gene.exons{k}(:,1))<75
      FLAG=1;
    end
  end
  if FLAG==1
    RESULT{2}={'gene shorter than 75 bp'};
    STAT{i}=RESULT;
    continue
  end
  
  
  COUNTER=3;
  V_COUNTER=1
  if (size(reads1,1)+size(reads2,1)>0)
    if(min(size(reads1,1),size(reads2,1))>0)
      PV=1;
      INFO='';
      [PV,INFO] =diff_mmd_variance(reads1,reads2,gene,CFG.VARIANCES{V_COUNTER},CFG.VARIANCES{V_COUNTER});
      RESULT{COUNTER}=PV;
      COUNTER=COUNTER+1;
    else
      PV=1;
      RESULT{COUNTER}=PV;
      COUNTER=COUNTER+1;
    end
  else
    PV=1;
    RESULT{COUNTER}=PV;
    COUNTER=COUNTER+1;
  end
      
  if and(size(reads1,1)+size(reads2,1)>0,size(gene.transcripts,2)>1)
    if(min(size(reads1,1),size(reads2,1))>0)   
      [P_VALUE, INFO]= diff_nbin7({reads1},{reads2},gene,SEQUENCED_LENGTH,CFG.VARIANCES{V_COUNTER},CFG.VARIANCES{V_COUNTER});
      RESULT{COUNTER}=P_VALUE;
      COUNTER=COUNTER+1;
    else
      PV=1;
      RESULT{COUNTER}=PV;
      COUNTER=COUNTER+1;
    end
  else
    PV=1;
    RESULT{COUNTER}=PV;
    COUNTER=COUNTER+1;
  end
  
  
  
  [SPLICINGEVENTS,SEQUENCE,EXONSEQUENCE]=splicingsequence(gene);
  [UNIQUE_NEW_EXONS,GRAPHNODES,ORDER_OF_GRAPHNODE,EIRS_IN_SEQ]=transform_single_end_reads(SPLICINGEVENTS,SEQUENCE,EXONSEQUENCE,CFG.SEQUENCED_LENGTH-CLIP);
  [NEW_READS1,UNEXPLAINED_READS1,UNEXPLAINED_INDEX1]= convert_reads_to_region_indicators(reads1,UNIQUE_NEW_EXONS,GRAPHNODES,ORDER_OF_GRAPHNODE,EIRS_IN_SEQ,gene);
  [NEW_READS2,UNEXPLAINED_READS2,UNEXPLAINED_INDEX2]=convert_reads_to_region_indicators(reads2,UNIQUE_NEW_EXONS,GRAPHNODES,ORDER_OF_GRAPHNODE,EIRS_IN_SEQ,gene);
  
  if and(size(reads1,1)+size(reads2,1)>0,size(gene.transcripts,2)>1)
    if(min(size(reads1,1),size(reads2,1))>0)
      [PV,INFO] =  diff_poisson_bonf_3_unequal_segment(NEW_READS1,NEW_READS2,gene,CFG.SEQUENCED_LENGTH-CLIP);
      RESULT{COUNTER}=PV;
      COUNTER=COUNTER+1;
    else
      PV=1;
      RESULT{COUNTER}=PV;
      COUNTER=COUNTER+1;
    end
  else
    PV=1;
    RESULT{COUNTER}=PV;
    COUNTER=COUNTER+1;
  end
  
  
  STAT{i}=RESULT;
  
end;
%keyboard
OUT_FILENAME=[CFG.out_base  SUFF_NAME '_' int2str(JB_NR)   '.mat'];
save(OUT_FILENAME,'STAT')
