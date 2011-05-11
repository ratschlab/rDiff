function rdiff_gendata(organism, fasta_fname, gff3_fname, alignment1_fname, alignment2_fname, info_fname)
% rdiff_gendata(organism, fasta_fname, gff3_fname, alignment1_fname, alignment2_fname, info_fname)

RDIFF_GALAXY_DIR='/home/galaxy/svn/projects/RNASeq_galaxy/rdiff.web/' ;

unix(sprintf('cp %s/examples/%s.fasta %s', RDIFF_GALAXY_DIR, organism, fasta_fname)) ;
unix(sprintf('cp %s/examples/%s.gff3 %s', RDIFF_GALAXY_DIR, organism, gff3_fname)) ;
unix(sprintf('cp %s/examples/%s-SRX001872.sam %s', RDIFF_GALAXY_DIR, organism, alignment1_fname)) ;
unix(sprintf('cp %s/examples/%s-SRX001875.sam %s', RDIFF_GALAXY_DIR, organism, alignment2_fname)) ;
unix(sprintf('cp %s/examples/%s.info %s', RDIFF_GALAXY_DIR, organism, info_fname)) ;

return

function ret=keyboard_allowed()
% ret=keyboard_allowed()
%
% returns 1, if a keyboard command would be allowed

global g_ignore_keyboard
global THIS_IS_A_RPROC_PROCESS  

if isequal(g_ignore_keyboard, 1),
  ret=0 ;
  return ;
end ;

if isequal(THIS_IS_A_RPROC_PROCESS, 1),
  ret=0 ;
  return ;
end ;

ret=1 ;


return
