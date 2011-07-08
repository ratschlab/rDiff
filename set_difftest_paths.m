function [difftest_paths] = set_difftest_paths()
% [difftest_paths] = set_difftest_paths()



difftest_paths ='.:./../samtools/:./../read_utils:./../variance:./../tests:./../mmd:./../nnmf:./../visualisation:/fml/ag-raetsch/share/svn/tools/genomes/:/fml/ag-raetsch/share/svn/tools/utils/';
;

addpath(difftest_paths);