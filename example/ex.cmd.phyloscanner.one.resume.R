#	setup output and working directories
work.dir			<- getwd()
out.dir				<- getwd()
#	link to phyloscanner.py
prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner.py'
#   prefix that identifies phyloscanner output
prefix.infiles      <- '/work/or105/phyloscanner/ptyr1_'
#	define phyloscanner arguments
pty.args			<- list(	prog.pty=prog.pty, 
								prog.mafft=NA, 
								prog.raxml=NA, 
								data.dir=NA, 
								work.dir=work.dir, 
								out.dir=out.dir, 
								alignments.file=system.file(package="phyloscan", "HIV1_compendium_AD_B_CPX_v2.fasta"),
								alignments.root='REF_CPX_AF460972', 
								alignments.pairwise.to='REF_B_K03455',
								window.automatic= '', 
								merge.threshold=0, 
								min.read.count=1, 
								quality.trim.ends=23, 
								min.internal.quality=23, 
								merge.paired.reads=TRUE, 
								no.trees=FALSE, 
								dont.check.duplicates=FALSE,
								num.bootstraps=1,
								all.bootstrap.trees=TRUE,
								strip.max.len=350, 
								min.ureads.individual=NA, 
								win=c(800,9400,25,250), 
								keep.overhangs=FALSE, 
								duplicated.raw.threshold=3,
								duplicated.ratio.threshold=1/200,				
								select=NA)
#	generate bash script
cmd                 <- phsc.cmd.phyloscanner.one.resume(prefix.infiles, pty.args)
#	print bash script to screen
cat(cmd)
#   this can be started eg through system on a UNIX machine, or written to file and then started manually (recommended)
#system(cmd)
#file.out <- 'ZZZ'
#cat(cmd, file=file.out)
