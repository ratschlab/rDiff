function [rDiff_paths] = set_rDiff_paths()
% [difftest_paths] = set_rDiff_paths()

%Initialize rDiff_paths
rDiff_paths='';


%ADDPATH('/FML/AG-RAETSCH/SHARE/SOFTWARE/MATLAB_TOOLS/RPROC/')
%ADDPATH('/FML/AG-RAETSCH/SHARE/SVN/TOOLS/UTILS');
%ADDPATH('/FML/AG-RAETSCH/HOME/DREWE/SVN/PROJECTS/RNASEQ/DIFFTEST/VARIANCE/FINAL');
%ADDPATH ../../EXPERIMENTAL
%addpath('/fml/ag-raetsch/share/software/matlab_tools/rproc/')
%addpath('/fml/ag-raetsch/share/svn/tools/utils');
%addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/variance/final');
%addpath ../../experimental
%addpath /fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/variance/nbin/
%addpath tools/read_utils/
%addpath('/fml/ag-raetsch/share/software/matlab_tools/rproc/')
%addpath('/fml/ag-raetsch/share/svn/tools/utils');
%addpath('~/svn/tools/utils'); 
%addpath('~/svn/tools/genomes');
%addpath ~/svn/projects/RNASeq/difftest/variance
%addpath('/fml/ag-raetsch/home/drewe/svn/projects/mGene_core/experiments/thaliana/');
%addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/experiments/');
%addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/read_utils/');
%addpath('/fml/ag-raetsch/home/drewe/svn/projects/RNASeq/difftest/tests/sequence_tools/');


%addpath('/fml/ag-raetsch/share/software/matlab_tools/rproc/')
%addpath('/fml/ag-raetsch/share/svn/tools/utils');
%addpath('~/svn/tools/utils'); 
%addpath('~/svn/tools/genomes'); 

%addpath(genpath('/fml/ag-raetsch/home/drewe/svn/tools/chronux/'))
%addpath ../../tests


%Add paths to rDiff_paths
rDiff_paths=[rDiff_paths ':' './tests/'];
rDiff_paths=[rDiff_paths ':' './tools/'];
rDiff_paths=[rDiff_paths ':' './tools/read_utils/'];
rDiff_paths=[rDiff_paths ':' './variance/'];
rDiff_paths=[rDiff_paths ':' genpath('./locfit/')];
rDiff_paths=[rDiff_paths ':' '/fml/ag-raetsch/home/drewe/svn/tools/utils/'];
rDiff_paths=[rDiff_paths ':' '/fml/ag-raetsch/home/drewe/svn/tools/rproc/'];
addpath(rDiff_paths);