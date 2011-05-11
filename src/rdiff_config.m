function rdiff_config
% RDIFF_CONFIG   Sets a few global variables with system dependent paths.
%
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Written (W) 2009-2011 Regina Bohnert, Gunnar Raetsch
% Copyright (C) 2009-2011 Max Planck Society
%


% rDiff paths
global RDIFF_PATH RDIFF_SRC_PATH

% interpreter paths
global INTERPRETER MATLAB_BIN_PATH OCTAVE_BIN_PATH

% SAMTools path
global SAMTOOLS_DIR

% configuration (adapt to the user's configuration)
RDIFF_PATH = getenv('RDIFF_PATH');
RDIFF_SRC_PATH = getenv('RDIFF_SRC_PATH');
INTERPRETER = getenv('INTERPRETER');
MATLAB_BIN_PATH = getenv('MATLAB_BIN_PATH');
OCTAVE_BIN_PATH = getenv('OCTAVE_BIN_PATH');
SAMTOOLS_DIR = getenv('SAMTOOLS_DIR');

% switch off a few expected warnings
addpath(sprintf('%s/tools', RDIFF_PATH));
engine = determine_engine();
if isequal(engine, 'octave'),
  warning('off', 'Octave:precedence-change');
  warning('off', 'Octave:function-name-clash');
  warning('off', '');
  warning('off', 'Octave:num-to-str');
  warning('off', 'Octave:function-name-clash');
  warning('off', 'Octave:divide-by-zero');
  warning('off', 'Octave:future-time-stamp');
  warning('off', 'Octave:assign-as-truth-value');
else
  warning('off', 'MATLAB:typeaheadBufferOverflow');
end

% make sure no process stops with a debug prompt
global g_ignore_keyboard
g_ignore_keyboard = 1;
