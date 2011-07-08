function []=get_gene_variance_caller_single(PAR)

samtools='/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/samtools/';
CFG = PAR.CFG;

genes = PAR.genes;
%BAM_FILES=CFG.BAM_FILES;
JB_NR=CFG.JB_NR;

%load ../POS.mat

clear PAR;

%%%% paths
%addpath(CFG.paths);
addpath('/fml/ag-raetsch/home/drewe/svn/projects/mGene_core/experiments/thaliana/');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/experiments/');
addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/read_utils/');
out_base=CFG.out_base;
%sim_data_base=CFG.sim_data_base;
real_data_base=CFG.real_data_base;

expr1_bam=CFG.expr1_bam;

%genes=PAR.genes;
STAT=cell(size(genes));
do_save = 0;
out_base = '/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/out';
save_dir=fullfile(out_base,'cache');
addpath ../../experimental


for i=1:size(genes,2)
  RESULT=cell(1,7);
  gene = genes(i);
  %set default return
  RESULT{1}=gene.name;
  RESULT{2}=JB_NR;
  RESULT{3}=[inf,inf];
  RESULT{4}=[inf,inf];
  RESULT{5}=[inf,inf];
  RESULT{6}=0;
  RESULT{7}=0;
  
  if isempty(gene.exons)
    STAT{i}=RESULT;
    continue;
  end

  
  CHR_PREFIX='';
  %[reads1] = get_reads_single(expr1_bam,gene,samtools,false,10,CHR_PREFIX);
  if strcmp(gene.strand,'-') 
    gene.start=gene.start+1;
    gene.stop=gene.stop+1;
    [reads1] = get_reads_single(expr1_bam,gene,samtools,false,10,CHR_PREFIX);
    [reads1,FLAG]=remove_reads_from_other_genes(reads1,gene);    
    gene.start=gene.start-1;
    gene.stop=gene.stop-1;
  else
    [reads1] = get_reads_single(expr1_bam,gene,samtools,false,10,CHR_PREFIX);
    [reads1,FLAG]=remove_reads_from_other_genes(reads1,gene);  
  end
  reads1=reads1(sum(reads1,2)>40,:);
 
  RESULT{3}=[size(reads1,1)]
  

  NR_OF_TRANS=size(gene.transcripts,2);
  [SPLICINGEVENTS,SEQUENCE,EXONSEQUENCE]=splicingsequence(gene);
  EXON_IDX=sum(EXONSEQUENCE,1)<NR_OF_TRANS;

 
  %Determine non-unique regions
  NON_UNIQUE_IDX_VECT_STARTS=zeros(1,size(reads1,2));%This vector is for the readstarts
  if not(isempty(gene.non_unique_regions))
    NON_UNIQUE_IDX_VECT=zeros(1,size(reads1,2));    
    REGIONS=gene.non_unique_regions;
    for j=1:size(gene.non_unique_regions,1)
      CURR_REGIONS=REGIONS(j,:);
      CURR_REGIONS=max([CURR_REGIONS;[gene.start,gene.start]],[],1);
      CURR_REGIONS=min([CURR_REGIONS;[gene.stop,gene.stop]],[],1);
      CURR_REGIONS=CURR_REGIONS-gene.start+1;
      NON_UNIQUE_IDX_VECT(CURR_REGIONS(1):CURR_REGIONS(2))=1;
      if (CURR_REGIONS(2)-CURR_REGIONS(1)+1)>CFG.SEQUENCED_LENGTH
	NON_UNIQUE_IDX_VECT_STARTS(CURR_REGIONS(1):(CURR_REGIONS(2)-CFG.SEQUENCED_LENGTH+1))=1;
      end
    end
    %find reads which could belong to other genes
    READ_IDX=sum(reads1(:,NON_UNIQUE_IDX_VECT>0),2)==sum(reads1,2);    
    COVERED_REGION_LENGTH=sum(sum(reads1(READ_IDX,:),1)>0);
    COVERED_REGION_IDX=(sum(reads1(READ_IDX,:),1)>0);
    reads1=reads1(not(READ_IDX),:);
  else 
    COVERED_REGION_LENGTH=0;
  end
  
  %Get starts of the reads
  [TEMP,START1]=max(reads1,[],2); 
  
  S1=size(reads1,1);
  new_reads1=sparse((1:S1)',START1,ones(S1,1),S1,size(reads1,2),S1);
   
  
  %Count the reads that start in the non-diferential exons
   
  
  %Determine nonspliced regions
  EXON_IDX2=sum(EXONSEQUENCE,1);
  IDX_VECT=zeros(1,size(reads1,2)); 
  OVERHEAD=CFG.SEQUENCED_LENGTH-1;
  STATE=0;  %This means, is alt spliced
  for j=length(EXON_IDX2):(-1):1
    if EXON_IDX2(j)==NR_OF_TRANS
      if SPLICINGEVENTS(j+1)-SPLICINGEVENTS(j)<OVERHEAD
	OVERHEAD=OVERHEAD-(SPLICINGEVENTS(j+1)-SPLICINGEVENTS(j)+1);
      else
	IDX_VECT(SPLICINGEVENTS(j):(SPLICINGEVENTS(j+1)-1-OVERHEAD))=1;
	OVERHEAD=0;
      end
    else
      if EXON_IDX2(j)>0
	OVERHEAD=CFG.SEQUENCED_LENGTH-1;
      end
    end    
  end


  S1=new_reads1(sum(new_reads1(:,IDX_VECT>0),2)>0,:);
  
  if or(isempty(S1),sum(IDX_VECT>0)-COVERED_REGION_LENGTH<1)
    RESULT{4}=[];
    RESULT{5}=0;
    RESULT{6}=[];
    continue;
  end
  
  
  IDX_VECT=and(IDX_VECT,not(NON_UNIQUE_IDX_VECT_STARTS));
  POSITIONS=max(sum(IDX_VECT>0),1);
  
  COUNT1=S1*(IDX_VECT'>0);
  RESULT{4}=[sum(COUNT1)]/POSITIONS;
  RESULT{5}=POSITIONS;
  RESULT{6}=sum(S1(:,IDX_VECT>0),1);
  
  [size(reads1,1)];
  
  clear reads1;
  STAT{i}=RESULT;
  
end;

OUT_FILENAME=[CFG.out_base CFG.outfile_prefix '.mat'];
    save(OUT_FILENAME,'STAT')
