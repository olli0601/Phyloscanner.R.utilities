#	setup input, output and working directories
pty.data.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
work.dir			<- getwd()
out.dir				<- getwd()
#	link to phyloscanner.py
prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner.py'
#   input file with reference to bams and refs
file.bam            <- 'XXX'
file.ref            <- 'YYY'
#	define phyloscanner arguments
pty.args			<- list(	prog.pty=prog.pty,
                                prog.mafft='mafft',
                                prog.raxml='"raxmlHPC-AVX -m GTRCAT"',
                                data.dir=pty.data.dir,
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
                                min.ureads.individual=NA,
                                win=c(800,2400,250,250),
                                keep.overhangs=FALSE,
                                duplicated.raw.threshold=3,
                                duplicated.ratio.threshold=1/200,
                                select=NA)
#	generate bash script
cmd                 <- phsc.cmd.phyloscanner.one(file.bam, file.ref, pty.args)
#	print bash script to screen
cat(cmd)
#   this can be started eg through system on a UNIX machine, or written to file and then started manually (recommended)
#system(cmd)
#file.out <- 'ZZZ'
#cat(cmd, file=file.out)
