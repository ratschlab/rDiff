function [VARIANCE]=get_gene_variance_biol_bin(CFG,genes,REPLICATES,FRAC_OF_OUTLRS,BIN_LEN,REGION)



addpath('/fml/ag-raetsch/share/software/matlab_tools/rproc/')
addpath('/fml/ag-raetsch/share/svn/tools/utils');
addpath('~/svn/tools/utils'); 
addpath('~/svn/tools/genomes');

%%%% paths

JOB_INFO = rproc_empty();
 
CFG.paths = set_difftest_paths()
%%%% configuration
CFG = configure_difftest(CFG);

nb_of_chunks=CFG.nb_of_chunks;
idx=[(1:size(genes,2))',ceil((1:size(genes,2))*nb_of_chunks/size(genes,2))'];

NR_OF_REPLICATES=sum(REPLICATES==1);
POS_OF_REPLICATES=find(REPLICATES);

ERROR_NR=[];
PROBS=cell(NR_OF_REPLICATES,size(genes,2));
for r=1:NR_OF_REPLICATES
  for i= 1:nb_of_chunks
    CFG.outfile_prefix=[CFG.BAM_FILES{POS_OF_REPLICATES(r)} '_' num2str(i) '_of_' num2str(nb_of_chunks) '_alt'];
    FILENAME=[CFG.out_base CFG.outfile_prefix '.mat'];
    try
      load(FILENAME,'STAT');
      temp_idx=idx(idx(:,2)==i,1);
      for j=1:size(temp_idx)  
	PROBS{r,temp_idx(j)} = STAT{j};
      end
    catch
      ERROR_NR=[ERROR_NR;i];
    end
  end
end

genes_count=cell(NR_OF_REPLICATES,size(genes,2));
COUNTS=zeros(NR_OF_REPLICATES,size(genes,2));
if REGION==1
  for r=1:NR_OF_REPLICATES
    for i=1:size(genes,2)
      T=PROBS{r,i};
      if isempty(T)
	COUNTS(r,i)=0;
	genes_count{r,i}=[];   
      else
	if not(isempty(T{4}))
	  COUNTS(r,i)=T{4}(1);
	  genes_count{r,i}=T{6};
	end
      end
    end
  end
else
  for r=1:NR_OF_REPLICATES
    for i=1:size(genes,2)
      T=PROBS{r,i};
      if isempty(T)
	COUNTS(r,i)=0;
	genes_count{r,i}=[];   
      else
	if not(isempty(T{7}))
	  COUNTS(r,i)=T{7}(1);
	  genes_count{r,i}=T{9};
	end
      end
    end
  end
end
COUNTS(isnan(COUNTS))=0;
COUNTS(isinf(COUNTS))=0;
COUNTS=COUNTS';

% HUBER ANSATZ

addpath(genpath('/fml/ag-raetsch/home/drewe/svn/tools/chronux/'))
addpath ../../tests

%CONDITION=COUNTS;

CONDITION=sum(COUNTS,2);

%Estimate Library size and create variance model

[S]=estimate_lib_size(COUNTS);

%LENGTH=0;
%for i=1:size(genes,2)
%  if and(not(isempty(genes_count{1,i})),not(isempty(genes_count{2,i})))
%    LENGTH=length(genes_count{1,i})+LENGTH;
%  end
%end


FRAC_OF_OUTLRS=0.0;
EXPR=COUNTS;
[VAL,POS]=sort(abs(log(EXPR(:,1)./EXPR(:,2))));
OUTL_IDX=floor(length(EXPR(:,1))*(1-FRAC_OF_OUTLRS));
  
%BIN_LEN=1;
LENGTH=0;
for i=POS(1:OUTL_IDX)'
  %genes_count{1,i}=genes_count{1,i}(3:(end-3));
  %genes_count{2,i}=genes_count{2,i}(3:(end-3));
end


for i=POS(1:OUTL_IDX)'
  if and(not(isempty(genes_count{1,i})),not(isempty(genes_count{2,i})))
    %g=[genes_count{1,i};genes_count{2,i}];
    %NON_ZERO=sum(g,[],1)>0;
    %genes_count{1,i}=genes_count{1,i}(NON_ZERO);
    %genes_count{2,i}=genes_count{2,i}(NON_ZERO);
        
    BINS=floor(length(genes_count{1,i})/BIN_LEN);
    T1=[];
    T2=[];
    if BIN_LEN>1
      for j=1:BINS
	T1=[T1 sum(genes_count{1,i}(((j-1)*BIN_LEN+1):(j*BIN_LEN)))];
	T2=[T2 sum(genes_count{2,i}(((j-1)*BIN_LEN+1):(j*BIN_LEN)))];
      end
      genes_count{1,i}=T1;
      genes_count{2,i}=T2;
      LENGTH=BINS+LENGTH;
    else
      LENGTH=length(genes_count{1,i})+LENGTH;
    end
  end
end
LENGTH
MEAN=zeros(1,LENGTH);
VAR=zeros(1,LENGTH);

SS=[];
COUNTER=0;
%for i=1:size(genes,2)
for i=POS(1:OUTL_IDX)'
  if and(not(isempty(genes_count{1,i})),not(isempty(genes_count{2,i})))
    if sum(genes_count{1,i})+sum(genes_count{2,i})>10
      LENGTH=length(genes_count{1,i});
      TS(1)=sum(genes_count{1,i});
      TS(2)=sum(genes_count{2,i});
      S(1)=2*TS(1)/(sum(TS));
      S(2)=2*TS(2)/(sum(TS));
      SS=[SS;S];S=[1,1];
      TEMP_MEAN=mean([genes_count{1,i}/S(1);genes_count{2,i}/S(2)]);
      TEMP_VAR=var([genes_count{1,i}/S(1);genes_count{2,i}/S(2)]);
      MEAN((COUNTER+1):(COUNTER+LENGTH))=TEMP_MEAN;
      VAR((COUNTER+1):(COUNTER+LENGTH))=TEMP_VAR;
      %if and(min(TEMP_MEAN)<50,max(TEMP_VAR)>4000),break,end
      COUNTER=COUNTER+LENGTH;
    end
  end
end

ZERO_IDX=MEAN>0;
MEAN=MEAN(ZERO_IDX);
VAR=VAR(ZERO_IDX);

%Estimate Library size and create variance model

SAMPLEBINS=100;
SAMPLES_PER_BIN=500;
[MEAN,POS]=sort(MEAN);
MAX_VAL=max(MEAN);
MIN_VAL=min(MEAN);

BOUNDS=exp(linspace(log(MIN_VAL),log(MAX_VAL),SAMPLEBINS+1));
SAMPLE=[]
for i=1:SAMPLEBINS
  NR_IN_BIN=length((find(MEAN>=BOUNDS(i),1,'first'):find(MEAN<BOUNDS(i+1),1,'last')));
  if NR_IN_BIN==0
    continue;
  elseif NR_IN_BIN<=SAMPLES_PER_BIN
    SAMPLE=[SAMPLE find(MEAN>=BOUNDS(i),1,'first'):find(MEAN<BOUNDS(i+1),1,'last')];
  else
    CHUNK_SAMPLE=find(MEAN>=BOUNDS(i),1,'first'):find(MEAN<BOUNDS(i+1),1,'last');
    TEMP_SAMPLE=randperm(NR_IN_BIN);
    SAMPLE=[SAMPLE CHUNK_SAMPLE(TEMP_SAMPLE(1:SAMPLES_PER_BIN))];
  end
end

  
VAR=VAR(POS);



%MAXP=length(MEAN);
%SAMPLE=randperm(length(MEAN));
%SAMPLE=SAMPLE(1:min(length(MEAN),100000));

[VARIANCE]=estimate_variance(ones(size(MEAN(SAMPLE)')),VAR(SAMPLE)',1,MEAN(SAMPLE)');


%clear MEAN
%clear VAR
%clear ZERO_IDX
%clear SAMPLE
%clear genes_count


pp=0.1:0.1:1000;

figure;
plot(MEAN(SAMPLE),VAR(SAMPLE),'.')
hold on
plot(pp',predict_variance(1,1,pp',VARIANCE),'r')
plot(pp,pp,'g')
xlim([0,500])
ylim([0,1000])

RR=[]; for i=1:1000,TIDX=and(MEAN>i,MEAN<i+1);RR=[RR mean(VAR(TIDX))];end 
plot(1:1000,RR,'k')
