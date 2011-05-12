function rdiff(anno_fname, track1, track2, out_fname, test_meth)
% RDIFF   Performs differential expression testing from RNA-Seq measurements.
%   rdiff(anno_fname, track1, track2, out_fname, test_meth)
%
%   -- input --
%   anno_fname: name of file containing genes
%   track1:     name of BAM file 1
%   track2:     name of BAM file 2
%   out_fname:  name of result file with p-values
%   test_meth:  test method ('poisson' or 'poisson_exp 'or 'mmd')
%
%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 3 of the License, or
%   (at your option) any later version.
%
%   Written (W) 2009-2011 Regina Bohnert
%   Copyright (C) 2009-2011 Max Planck Society
%


% rDiff paths
global RDIFF_PATH RDIFF_SRC_PATH

% interpreter paths
global INTERPRETER MATLAB_BIN_PATH OCTAVE_BIN_PATH

% SAMTools path
global SAMTOOLS_DIR

addpath(sprintf('%s/mex', RDIFF_PATH));
addpath(sprintf('%s/tools', RDIFF_PATH));
addpath(sprintf('%s', RDIFF_SRC_PATH));

samtools_dir = sprintf('%s/', SAMTOOLS_DIR);

if ~exist(sprintf('%s.bai', track1), 'file')
  command = sprintf('%s./samtools index %s', samtools_dir, track1);
  [s m] = unix(command);
  if ~exist(sprintf('%s.bai', track1), 'file')
    sprintf('\nbai file for %s could not be created\n', track1); 
  end
end

if ~exist(sprintf('%s.bai', track2), 'file')
  command = sprintf('%s./samtools index %s', samtools_dir, track2);
  [s m] = unix(command);
  if ~exist(sprintf('%s.bai', track2), 'file')
    sprintf('\nbai file for %s could not be created\n', track2); 
  end
end

[ret, timedate] = unix('date');
timedate(timedate==sprintf('\n')) = [];
fprintf(1, '\n*** rDiff started with %s test %s *** \n\n', test_meth, timedate);

%%%% load genes %%%%%
load(anno_fname, 'genes');

%%%% rDiff %%%%%
fprintf(1, 'testing %i genes for differential expression...\n', length(genes));
p_values = eval(sprintf('%s_test(genes, track1, track2)', test_meth)); 
fprintf(1, 'done.\n');

%%%%% write to txt file %%%%%
write_rdiff_result(genes, p_values, out_fname, test_meth);

[ret, timedate] = unix('date');
timedate(timedate==sprintf('\n')) = [];
fprintf(1, '\n*** rDiff finished with %s test %s *** \n\n', test_meth, timedate);