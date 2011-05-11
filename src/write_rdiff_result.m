function write_rdiff_result(genes, p_values, out_fname, test_meth) 
% WRITE_RDIFF_RESULT   Writes rDiff results to file.
%
%   write_rdiff_result(genes, out_fname, source, mapped_reads, read_len) 
%
%   -- input --
%   genes:     struct defining genes with start, stops, exons etc.
%   p_values:  p-values of test
%   out_fname: name of result file with p-values
%   test_meth: test method ('poisson' or 'mmd')
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


if exist(out_fname, 'file')
  fprintf('replacing file %s...\n', out_fname);
  [fd msg] = fopen(out_fname, 'w+');
  disp(msg);
else
  fprintf('creating file %s...\n', out_fname);
  [fd msg] = fopen(out_fname, 'w+');
  disp(msg);
end

fprintf(fd, '# gene name\ttesting method\tp-value\n');
for g = 1:length(genes),
  fprintf(fd, '%s\t%s\t%.4f\n', genes(g).name, test_meth, p_values(g));
end

fclose(fd);
