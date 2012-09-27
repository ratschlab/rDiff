function [rDiff_paths] = set_rDiff_paths()
% [difftest_paths] = set_rDiff_paths()

%Initialize rDiff_paths
rDiff_paths='';

%Add paths to rDiff_paths
rDiff_paths=[rDiff_paths ':' './tests/'];
rDiff_paths=[rDiff_paths ':' './tools/'];
rDiff_paths=[rDiff_paths ':' './tools/read_utils/'];
rDiff_paths=[rDiff_paths ':' './variance/'];
rDiff_paths=[rDiff_paths ':' genpath('./chronux/')];
rDiff_paths=[rDiff_paths ':' '/fml/ag-raetsch/home/drewe/svn/tools/utils/'];
rDiff_paths=[rDiff_paths ':' '/fml/ag-raetsch/home/drewe/svn/tools/rproc/'];
addpath(rDiff_paths);