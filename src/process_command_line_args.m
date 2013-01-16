function [CFG]=process_command_line_args(CFG,ARGS)

%This function adds the command line arguments to the CFG config variable

%Parse the ARGS 
    
%Check which interpreter is used
    if size(ver('Octave'),1)
        INTERPR = 1;
    else
        INTERPR = 0;
    end
    if INTERPR
        ARGS = strsplit(ARGS,':');
    else
        ARGS=regexp(ARGS,':','split');
    end
    
    
    %Process Bamfiles
    if INTERPR
        BAMS1 = strsplit(ARGS{2},',');
        BAMS2 = strsplit(ARGS{3},',');
    else
        BAMS1=regexp(ARGS{2},',','split');
        BAMS2=regexp(ARGS{3},',','split');
    end
   
    CFG.BAM_FILES={BAMS1{:},BAMS2{:}};
     
    %Name of the experiment. Use the FILENAMES if the entries are empty.
    CFG.NAMES=CFG.BAM_FILES;
    
    
    % Give the directory where the input-files are
    CFG.data_dir = [ARGS{7} '/'];
    
    % Indicate to which sample the bam-files belong
    CFG.SAMPLES=[repmat(1,1,size(BAMS1,2)),repmat(2,1,size(BAMS2,2))];

    % Location of the gene structure
    CFG.genes_path=ARGS{1};
    
    % Output directory
    
    CFG.out_base =  [ARGS{6} '/'];
    
    % Output directory for temporary files
    CFG.out_base_temp =   [CFG.out_base '/temp/'];
    mkdir(CFG.out_base_temp);

    if strcmp(ARGS{5},'poisson')
        CFG.perform_poisson=1;
    end
    if strcmp(ARGS{5},'param')
        CFG.perform_parametric=1;
    end
    if strcmp(ARGS{5},'nonparam')
        CFG.perform_nonparametric=1;
    end
return
