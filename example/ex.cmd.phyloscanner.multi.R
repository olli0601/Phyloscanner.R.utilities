# 	load data.table of predefined phyloscanner runs
data(runs_rakai_couples)
#	setup input, output and working directories
pty.data.dir		<- '/Users/Oliver/duke/2016_PANGEAphylotypes/data'
work.dir			<- getwd()
out.dir				<- getwd()
#	link to phyloscanner.py
prog.pty			<- '/work/or105/libs/phylotypes/phyloscanner.py'
#	produce bash scripts for the first two phyloscanner runs only
pty.select			<- 1:2
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
                                select=pty.select)
#	generate bash scripts
pty.c				<- phsc.cmd.phyloscanner.multi(runs_rakai_couples, pty.args)
#	print first bash script to screen
pty.c[1,cat(CMD)]
