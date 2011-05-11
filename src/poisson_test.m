function [p_values]=poisson_test(genes,bam1,bam2)
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2011 Philipp Drewe
%   Copyright (C) 2009-2011 Max Planck Society
%


min_read_length = 0;

p_values=ones(size(genes,2),1);

for j=1:size(genes,2)
  fprintf(1, 'gene %i: %s\n', j, genes(j).name);
  gene=genes(j);
  eidx = [gene.start,gene.stop];
  size(genes(j).transcripts,2)>1;
 
    if size(genes(j).transcripts,2)>1
      [mask1, read_intron_list] = get_reads(bam1, gene.chr,gene.start,gene.stop, '0');
      if isempty(mask1)
	reads1=zeros(0,gene.stop-gene.start+1);
      else
	reads1=sparse(mask1(1,:)',mask1(2,:)',ones(size(mask1,2),1),max(mask1(1,:)),gene.stop-gene.start+1);
      end
      
      [mask2, read_intron_list] = get_reads(bam2, gene.chr,gene.start,gene.stop, '0');
      if isempty(mask2)
	reads2=zeros(0,gene.stop-gene.start+1);
      else
	reads2=sparse(mask2(1,:)',mask2(2,:)',ones(size(mask2,2),1),max(mask2(1,:)),gene.stop-gene.start+1);
      end

    N1 = size(reads1,1);
    N2 = size(reads2,1);
  
  
    if (N1<N2) 
        RP = randperm(N2); 
	RP = RP(1:N1); 
	reads2 = reads2(RP,:); 
    elseif (N2<N1) 
        RP = randperm(N1); 
	RP = RP(1:N2); 
	reads1 = reads1(RP,:); 
    end; 
    
  
    [P_VALUE, STRUCT]= diff_poisson_bonf_2(reads1,reads2,gene);
    p_values(j)=P_VALUE;
    else
      p_values(j)=1; 
    end
end
