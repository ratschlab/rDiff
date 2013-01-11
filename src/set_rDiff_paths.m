function [rDiff_paths] = set_rDiff_paths()
% [difftest_paths] = set_rDiff_paths()

%Initialize rDiff_paths
rDiff_paths='';

%get the path to the source files
BASEPATH_SRC=getenv('RDIFF_SRC_PATH')


%Add paths to rDiff_paths
rDiff_paths=[rDiff_paths ':' [BASEPATH_SRC '/tests/']];
rDiff_paths=[rDiff_paths ':' [BASEPATH_SRC '/tools/']];
rDiff_paths=[rDiff_paths ':' [BASEPATH_SRC '/tools/read_utils/']];
rDiff_paths=[rDiff_paths ':' [BASEPATH_SRC '/variance/']];
rDiff_paths=[rDiff_paths ':' genpath([BASEPATH_SRC '/locfit/'])];
addpath(rDiff_paths);
