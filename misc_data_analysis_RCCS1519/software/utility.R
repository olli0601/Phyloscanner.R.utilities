
#' @export
#' @title Generate bash commands for multiple phyloscanner runs
#' @param pty.runs Data.table of individual assignments to phyloscanner runs, with columns 'PTY_RUN' (run id), 'SAMPLE_ID' (ID of individuals that are assigned to that run). Optional columns: 'RENAME_ID' (new ID for each bam file in phyloscanner output).
#' @param pty.args List of phyloscanner input variables. See examples.
#' @return Data.table with columns 'PTY_RUN' (run id) and 'CMD' (bash commands for that run). 
#' @description This function generates bash commands for multiple phyloscanner runs, that can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.multi.R  
phsc.cmd.phyloscanner.multi <- function(pty.runs, pty.args) 		
{
  stopifnot( any(colnames(pty.runs)=='PTY_RUN') )
  stopifnot( any(colnames(pty.runs)=='SAMPLE_ID') )
  stopifnot( any(colnames(pty.runs)=='BAM') )
  stopifnot( any(colnames(pty.runs)=='REF') )	
  ptyd	<- copy(pty.runs)		
  if(!any(is.na(pty.args[['select']])))
    ptyd<- subset(ptyd, PTY_RUN%in%pty.args[['select']])			
  #	if background alignment is specified in pty.runs, use it
  if(any(colnames(ptyd)=='BACKGROUND_ID'))
  {
    tmp	<- unique(data.table(BACKGROUND_ID=ptyd[, BACKGROUND_ID], FILE=file.path(pty.args[['data.dir']], ptyd[, BACKGROUND_ID])))
    #	check if files exist
    tmp	<- tmp[, list(EXISTS=file.exists(FILE)), by=c('BACKGROUND_ID','FILE')]
    if(any(tmp[, !EXISTS]))
      warning('\nCould not find location of BACKGROUND files for all runs in pty.runs, n=', tmp[, length(which(!EXISTS))],'\nRuns with missing background alignments are ignored. Please check.')
    ptyd <- merge(ptyd, subset(tmp, EXISTS, c(BACKGROUND_ID, FILE)), by='BACKGROUND_ID')
    set(ptyd, NULL, 'BACKGROUND_ID', NULL)
    setnames(ptyd, 'FILE', 'BACKGROUND_ID')
  }
  setkey(ptyd, PTY_RUN)
  #	add legacy arguments
  pty.args['all.bootstrap.trees']	<- FALSE
  pty.args['num.bootstraps']	<- 0		
  #	write pty.run files and get pty command lines
  pty.c		<- ptyd[, {
    #	PTY_RUN<- z <- 1; BAM<- subset(ptyd, PTY_RUN==z)[, BAM]; REF<- subset(ptyd, PTY_RUN==z)[, REF]
    #	SAMPLE_ID<- subset(ptyd, PTY_RUN==z)[, SAMPLE_ID]; RENAME_ID<- subset(ptyd, PTY_RUN==z)[, RENAME_ID]; BACKGROUND_ID<- subset(ptyd, PTY_RUN==z)[, BACKGROUND_ID]
    file.input		<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_input.csv',sep=''))
    tmp				<- cbind(BAM[!is.na(BAM)&!is.na(REF)], REF[!is.na(BAM)&!is.na(REF)])
    if(exists('RENAME_ID'))
      tmp			<- cbind(tmp, RENAME_ID[!is.na(BAM)&!is.na(REF)])
    write.table(tmp, file=file.input, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')
    file.patient	<- file.path(pty.args[['work.dir']], paste('ptyr',PTY_RUN,'_patients.txt',sep=''))
    tmp				<- paste0(SAMPLE_ID[!is.na(BAM)&!is.na(REF)],'.bam')
    if(exists('RENAME_ID'))
      tmp			<- unique(RENAME_ID[!is.na(BAM)&!is.na(REF)])
    if(exists('UNIT_ID'))
      tmp			<- unique(UNIT_ID[!is.na(BAM)&!is.na(REF)])				
    cat( paste(tmp,collapse='\n'), file= file.patient	)
    if(exists('BACKGROUND_ID'))
    {
      if(length(unique(BACKGROUND_ID))!=1)
        stop('\nrun',PTY_RUN,' Expected one background alignment file, found: ',paste(unique(BACKGROUND_ID), collapse=''),'. Please check.')
      pty.args$alignments.file<- unique(BACKGROUND_ID)
    }				
    cmd			<- phsc.cmd.phyloscanner.one(	pty.args, 
                                        file.input, 
                                        file.patient)
    #cmd			<- paste(cmd, pty.cmd.evaluate.fasta(pty.args[['out.dir']], strip.max.len=pty.args[['strip.max.len']], select=paste('^ptyr',PTY_RUN,'_In',sep=''), min.ureads.individual=pty.args[['min.ureads.individual']]), sep='')
    #cat(cmd)
    list(CMD= cmd)				
  },by='PTY_RUN']
  pty.c
}	

#' @export
#' @title Generate bash command for a single phyloscanner run
#' @param pty.args List of phyloscanner input variables. See examples.
#' @param file.input File name of the file that contains the list of bam files, reference files, and potentially aliases
#' @param file.patient File name of the file that contains the list of unique individuals/units to which the bam files correspond. Multiple bam files are allowed per individual/unit. 
#' @return Character string of phyloscanner commands.
#' @description This function generates bash commands for a single phyloscanner run, that can be called via 'system' in R, or written to file to run on a UNIX system.
#' @example example/ex.cmd.phyloscanner.one.R    
phsc.cmd.phyloscanner.one<- function(pty.args, file.input, file.patient)
{	
  stopifnot(is.character(file.input),is.character(file.patient))
  #	copy input variables into namespace	 
  attach(pty.args)	
  #	sense checks
  stopifnot(is.logical(all.bootstrap.trees))	
  #	define window coordinates
  window.coord			<- integer(0)
  if(!nchar(pty.args[['window.automatic']]))
  {
    stopifnot(length(pty.args[['win']])==4)
    window.coord		<- seq(pty.args[['win']][1], pty.args[['win']][2]-max(pty.args[['win']][3:4]),pty.args[['win']][3])
    window.coord		<- as.vector(rbind( window.coord,window.coord-1+pty.args[['win']][4] ))											
  }	
  #	
  if(!nchar(window.automatic))	stopifnot( is.numeric(window.coord), !length(window.coord)%%2)
  if(nchar(window.automatic))		stopifnot( !length(window.coord) )
  merge.paired.reads			<- ifelse(!is.na(merge.paired.reads) & merge.paired.reads, '--merge-paired-reads', NA_character_)
  keep.overhangs				<- ifelse(!is.na(keep.overhangs) & keep.overhangs, '--keep-overhangs', NA_character_)
  no.trees					<- ifelse(!is.na(no.trees) & no.trees, '--no-trees', NA_character_)
  num.bootstraps				<- ifelse(is.na(no.trees) & !is.na(num.bootstraps) & is.numeric(num.bootstraps) & num.bootstraps>1, as.integer(num.bootstraps), NA_integer_)
  dont.check.duplicates		<- ifelse(!is.na(dont.check.duplicates) & dont.check.duplicates, '--dont-check-duplicates', NA_character_)
  dont.check.recombination	<- ifelse(!is.na(dont.check.recombination) & dont.check.recombination==FALSE, '--check-recombination', NA_character_)	
  #	create local tmp dir
  cmd		<- paste("CWD=$(pwd)\n",sep='\n')
  cmd		<- paste(cmd,"echo $CWD\n",sep='')
  tmpdir	<- paste('pty','_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')	
  tmpdir	<- paste("$CWD/",tmpdir,sep='')
  cmd		<- paste(cmd,'mkdir -p "',tmpdir,'"\n',sep='')
  #	copy files to local tmp dir
  cmd		<- paste(cmd,'cp "',file.input,'" "',tmpdir,'"\n',sep='')	
  cmd		<- paste(cmd,'cp "',file.patient,'" "',tmpdir,'"\n',sep='')
  #	cd to tmp dir
  cmd		<- paste(cmd, 'cd "',tmpdir,'"\n', sep='')	
  cmd		<- paste(cmd, prog.pty,' "',basename(file.input),'" ',sep='')	
  cmd		<- paste(cmd, '--merging-threshold-a', merge.threshold,'--min-read-count',min.read.count,'--quality-trim-ends', quality.trim.ends, '--min-internal-quality',min.internal.quality,'--keep-output-together')
  if(!is.na(merge.paired.reads))
    cmd	<- paste(cmd,' ',merge.paired.reads,sep='')	
  if(!is.na(dont.check.duplicates))
    cmd	<- paste(cmd,' ',dont.check.duplicates,sep='')
  if(!is.na(dont.check.recombination))
    cmd	<- paste(cmd,' ',dont.check.recombination,sep='')	
  if(!is.na(alignments.pairwise.to))		
    cmd	<- paste(cmd,' --pairwise-align-to ',alignments.pairwise.to,sep='')
  if(!is.na(default.coord)){cmd <- paste(cmd, " --excision-ref ", alignments.pairwise.to,sep='')}
  if(nchar(window.automatic))
    cmd	<- paste(cmd,' --auto-window-params ', window.automatic,sep='')
  if(!nchar(window.automatic))
    cmd	<- paste(cmd,' --windows ', paste(as.character(window.coord), collapse=','),sep='')
  if(!is.na(alignments.file))
    cmd	<- paste(cmd,' --alignment-of-other-refs "',alignments.file,'"',sep='')	
  if(!is.na(no.trees))		
    cmd	<- paste(cmd, no.trees)
  if(!is.na(num.bootstraps))
    cmd	<- paste(cmd, ' --num-bootstraps ', num.bootstraps, sep='')
  if(!is.na(discard.improper.pairs) & discard.improper.pairs==TRUE)
    cmd	<- paste(cmd, ' --discard-improper-pairs ', sep = '')
  if(!is.na(keep.overhangs))
    cmd	<- paste(cmd, keep.overhangs)
  cmd		<- paste(cmd, '--x-mafft',prog.mafft)
  if(is.na(no.trees))
    cmd	<- paste(cmd, '--x-raxml',prog.raxml)	
  if(!is.na(default.coord)){
    if(default.coord){	  
      cmd	<- paste0(cmd, " --excision-coords \'823,824,825,892,893,894,907,908,909,1012,1013,1014,1156,1157,1158,1384,1385,1386,1444,1445,1446,1930,1931,1932,1957,1958,1959,2014,2015,2016,2023,2024,2025,2080,2081,2082,2134,2135,2136,2191,2192,2193,2280,2281,2282,2283,2284,2285,2298,2299,2300,2310,2311,2312,2316,2317,2318,2319,2320,2321,2322,2323,2324,2340,2341,2342,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2373,2374,2375,2379,2380,2381,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2400,2401,2402,2409,2410,2411,2412,2413,2414,2415,2416,2417,2424,2425,2426,2430,2431,2432,2436,2437,2438,2439,2440,2441,2442,2443,2444,2457,2458,2459,2460,2461,2462,2463,2464,2465,2469,2470,2471,2472,2473,2474,2478,2479,2480,2481,2482,2483,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2514,2515,2516,2517,2518,2519,2520,2521,2522,2526,2527,2528,2529,2530,2531,2535,2536,2537,2670,2671,2672,2679,2680,2681,2703,2704,2705,2709,2710,2711,2733,2734,2735,2742,2743,2744,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2769,2770,2771,2772,2773,2774,2778,2779,2780,2811,2812,2813,2814,2815,2816,2817,2818,2819,2823,2824,2825,2841,2842,2843,2847,2848,2849,2850,2851,2852,2856,2857,2858,2865,2866,2867,2871,2872,2873,2892,2893,2894,2895,2896,2897,2901,2902,2903,2904,2905,2906,2952,2953,2954,2961,2962,2963,3000,3001,3002,3015,3016,3017,3018,3019,3020,3030,3031,3032,3042,3043,3044,3084,3085,3086,3090,3091,3092,3099,3100,3101,3111,3112,3113,3117,3118,3119,3135,3136,3137,3171,3172,3173,3177,3178,3179,3180,3181,3182,3189,3190,3191,3192,3193,3194,3204,3205,3206,3210,3211,3212,3222,3223,3224,3228,3229,3230,3237,3238,3239,3246,3247,3248,3249,3250,3251,3255,3256,3257,3261,3262,3263,3396,3397,3398,3501,3502,3503,3546,3547,3548,3705,3706,3707,4425,4426,4427,4449,4450,4451,4503,4504,4505,4518,4519,4520,4590,4591,4592,4641,4642,4643,4647,4648,4649,4656,4657,4658,4668,4669,4670,4671,4672,4673,4692,4693,4694,4722,4723,4724,4782,4783,4784,4974,4975,4976,5016,5017,5018,5067,5068,5069,", paste(seq(6615,6811,by=1),collapse = ','),  ",", paste(seq(7110,7636,by=1),collapse = ',') ,",7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7873,7874,7875,7876,7877,7881,7882,7883,7884,7885,7886,", paste(seq(9400,9719,by=1),collapse=','),"\' ")
    }else{
      cmd	<- paste(cmd, " --excision-coords \'823,824,825,892,893,894,907,908,909,1012,1013,1014,1156,1157,1158,1384,1385,1386,1444,1445,1446,1930,1931,1932,1957,1958,1959,2014,2015,2016,2023,2024,2025,2080,2081,2082,2134,2135,2136,2191,2192,2193,2280,2281,2282,2283,2284,2285,2298,2299,2300,2310,2311,2312,2316,2317,2318,2319,2320,2321,2322,2323,2324,2340,2341,2342,2346,2347,2348,2349,2350,2351,2352,2353,2354,2355,2356,2357,2358,2359,2360,2373,2374,2375,2379,2380,2381,2385,2386,2387,2388,2389,2390,2391,2392,2393,2394,2395,2396,2400,2401,2402,2409,2410,2411,2412,2413,2414,2415,2416,2417,2424,2425,2426,2430,2431,2432,2436,2437,2438,2439,2440,2441,2442,2443,2444,2457,2458,2459,2460,2461,2462,2463,2464,2465,2469,2470,2471,2472,2473,2474,2478,2479,2480,2481,2482,2483,2496,2497,2498,2499,2500,2501,2502,2503,2504,2505,2506,2507,2514,2515,2516,2517,2518,2519,2520,2521,2522,2526,2527,2528,2529,2530,2531,2535,2536,2537,2670,2671,2672,2679,2680,2681,2703,2704,2705,2709,2710,2711,2733,2734,2735,2742,2743,2744,2748,2749,2750,2751,2752,2753,2754,2755,2756,2757,2758,2759,2769,2770,2771,2772,2773,2774,2778,2779,2780,2811,2812,2813,2814,2815,2816,2817,2818,2819,2823,2824,2825,2841,2842,2843,2847,2848,2849,2850,2851,2852,2856,2857,2858,2865,2866,2867,2871,2872,2873,2892,2893,2894,2895,2896,2897,2901,2902,2903,2904,2905,2906,2952,2953,2954,2961,2962,2963,3000,3001,3002,3015,3016,3017,3018,3019,3020,3030,3031,3032,3042,3043,3044,3084,3085,3086,3090,3091,3092,3099,3100,3101,3111,3112,3113,3117,3118,3119,3135,3136,3137,3171,3172,3173,3177,3178,3179,3180,3181,3182,3189,3190,3191,3192,3193,3194,3204,3205,3206,3210,3211,3212,3222,3223,3224,3228,3229,3230,3237,3238,3239,3246,3247,3248,3249,3250,3251,3255,3256,3257,3261,3262,3263,3396,3397,3398,3501,3502,3503,3546,3547,3548,3705,3706,3707,4425,4426,4427,4449,4450,4451,4503,4504,4505,4518,4519,4520,4590,4591,4592,4641,4642,4643,4647,4648,4649,4656,4657,4658,4668,4669,4670,4671,4672,4673,4692,4693,4694,4722,4723,4724,4782,4783,4784,4974,4975,4976,5016,5017,5018,5067,5068,5069,7863,7864,7865,7866,7867,7868,7869,7870,7871,7872,7873,7874,7875,7876,7877,7881,7882,7883,7884,7885,7886\' ",sep='')
    }
  }
  cmd		<- paste(cmd, '\n')
  run.id	<- gsub('_input.csv','',basename(file.input))
  #	process RAxML files
  if(is.na(no.trees) & (is.na(num.bootstraps) | (!is.na(num.bootstraps) & all.bootstrap.trees)))
    cmd	<- paste(cmd, 'for file in RAxML_bestTree*.tree; do\n\tmv "$file" "${file//RAxML_bestTree\\./',run.id,'_}"\ndone\n',sep='')
  if(is.na(no.trees) & !is.na(num.bootstraps) & !all.bootstrap.trees)
    cmd	<- paste(cmd, 'for file in RAxML_bipartitions.MLtreeWbootstraps*.tree; do\n\tmv "$file" "${file//RAxML_bipartitions.MLtreeWbootstraps/',run.id,'_}"\ndone\n',sep='')	
  #cmd	<- paste(cmd, "for file in AlignedReads*.fasta; do\n\tsed 's/<unknown description>//' \"$file\" > \"$file\".sed\n\tmv \"$file\".sed \"$file\"\ndone\n",sep='')		
  if(!is.na(alignments.file) & !is.na(keep.overhangs))
  {
    cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tcat "$file" | awk \'{if (substr($0,1,4) == ">REF") censor=1; else if (substr($0,1,1) == ">") censor=0; if (censor==0) print $0}\' > NoRef$file\ndone\n', sep='')		
    cmd	<- paste(cmd, 'for file in NoRefAlignedReads*.fasta; do\n\t',phsc.cmd.mafft.add(alignments.file,'"$file"','Ref"$file"', options='--keeplength --memsave --parttree --retree 1'),'\ndone\n',sep='')		
    cmd	<- paste(cmd, 'for file in RefNoRefAlignedReads*.fasta; do\n\t','mv "$file" "${file//RefNoRefAlignedReads/',run.id,'_}"\ndone\n',sep='')		
  }
  if(realignment==TRUE){
    cmd	<- paste0(cmd, 'echo Performing realignment...\n')
    cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\t', gsub('"','',mafft.opt) , '"$file" > "${file//.fasta/_v2.fasta}" \n done \n', sep='')
  }
  if(is.na(alignments.file) || is.na(keep.overhangs))
  {
    cmd	<- paste(cmd, 'for file in AlignedReads*.fasta; do\n\tmv "$file" "${file//AlignedReads/',run.id,'_}"\ndone\n',sep='')	
  }	
  #	move Duplicate Read Counts - only for backward compatibility
  #cmd	<- paste(cmd, 'for file in DuplicateReadCountsProcessed_*.csv; do\n\tmv "$file" "${file//DuplicateReadCountsProcessed_/',run.id,'_DuplicateReadCounts_}"\ndone',sep='')
  #	run phyloscanner tools and compress output
  if(is.na(no.trees))
    cmd	<- paste(cmd, phsc.cmd.process.phyloscanner.output.in.directory(tmpdir, file.patient, pty.args), sep='\n')
  #
  out.dir2<- out.dir	
  if(!is.na(no.trees))
  {
    out.dir2<- file.path(out.dir,paste0(run.id,'_trees'))
    cmd		<- paste(cmd,'\nmkdir -p ',out.dir2)		
  }
  cmd		<- paste(cmd, '\ncp ',run.id,'* "',out.dir2,'"\n',sep='')	
  #	zip up everything else
  tmp		<- ''
  if(length(window.coord)==2)
    tmp	<- window.coord[1]
  if(is.null(mem.save) || is.na(mem.save) || mem.save==0)
  {
    cmd		<- paste(cmd, 'for file in *; do\n\tzip -ur9XTjq ',paste(run.id,'_otherstuff',tmp,'.zip',sep=''),' "$file"\ndone\n',sep='')
    cmd		<- paste(cmd, 'cp ',paste(run.id,'_otherstuff',tmp,'.zip',sep=''),' "',out.dir2,'"\n',sep='')		
  }
  #	clean up
  cmd		<- paste(cmd,'cd $CWD\nrm -r "',tmpdir,'"\n',sep='')
  detach(pty.args)
  cmd
}


#' @export
#' @title Produce a single iqtree shell command. 
#' @return	Character string
cmd.iqtree<- function(infile.fasta, outfile=infile.fasta, pr=PR, pr.args='-m GTRCAT --HKY85 -p 42')
{		
  cmd<- paste("#######################################################
  # start: IQTREE
  #######################################################\n",sep='')
  cmd<- paste(cmd,"CWD=$(pwd)\n",sep='')
  cmd<- paste(cmd,"echo $CWD\n",sep='')	
  tmpdir.prefix	<- paste('rx_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"),sep='')
  tmpdir			<- paste("$CWD/",tmpdir.prefix,sep='')
  tmp.in			<- basename(infile.fasta)
  tmp.out			<- basename(outfile)	
  cmd<- paste(cmd,"mkdir -p ",tmpdir,'\n',sep='')
  cmd<- paste(cmd,'cp "',infile.fasta,'" ',file.path(tmpdir,tmp.in),'\n', sep='')	
  cmd<- paste(cmd,'cd "',tmpdir,'"\n', sep='')	
  cmd<- paste(cmd, pr,' ',pr.args,' -s ', tmp.in, ' -pre ',tmp.out,'\n', sep='') 
  cmd<- paste(cmd, "rm ", tmp.in,'\n',sep='')	
  cmd	<- paste(cmd, 'cp ',paste0(basename(outfile),'.iqtree'),' "',dirname(outfile),'"\n',sep='')
  cmd	<- paste(cmd, 'cp ',paste0(basename(outfile),'.treefile'),' "',dirname(outfile),'"\n',sep='')
  cmd<- paste(cmd, 'for file in *; do\n\tzip -ur9XTjq ',basename(outfile),'.zip "$file"\ndone\n',sep='')	
  cmd<- paste(cmd, 'cp ',basename(outfile),'.zip "',dirname(outfile),'"\n',sep='')
  cmd<- paste(cmd,'cd $CWD\n', sep='')
  cmd<- paste(cmd, "rm -r ", tmpdir,'\n',sep='')
  cmd<- paste(cmd, "#######################################################
  # end: IQTREE
  #######################################################\n",sep='')
  cmd
}


#' @export
cmd.hpcwrapper.cx1.ic.ac.uk<- function(hpc.select=1, hpc.walltime=24, hpc.mem=HPC.MEM, hpc.nproc=1, hpc.q=NA, hpc.load=HPC.CX1.IMPERIAL.LOAD, hpc.array=NA)
{
  wrap<- "#!/bin/sh"
  tmp	<- paste("#PBS -l walltime=",hpc.walltime,":59:59,pcput=",hpc.walltime,":45:00",sep='')
  wrap<- paste(wrap, tmp, sep='\n')		
  tmp	<- paste("#PBS -l select=",hpc.select,":ncpus=",hpc.nproc,":mem=",hpc.mem,sep='')
  wrap<- paste(wrap, tmp, sep='\n')
  wrap<- paste(wrap, "#PBS -j oe", sep='\n')
  if(!is.na(hpc.array))
    wrap<- paste(wrap, "\n#PBS -J 1-", hpc.array, sep='')
  if(!is.na(hpc.q))
    wrap<- paste(wrap, "\n#PBS -q",hpc.q, sep='')
  wrap<- paste(wrap, hpc.load, sep='\n')
  wrap
}


#' @export
#' @title Generate bash commands to extract queue status information.
#' @param jobid ID pointing to the job we want to query for.
#' @param var_name Name of the variable where we want to store the bash output.
#' @return Character string containing the line of code generated
#' @description This function produces code to extract the number of subjobs that are either queued, running or finished according to the qstat command. Additionaly, the output number can be stored in a bash variable.
#' @example  a <- qstat.info(option='f', jobid='124', var_name='Q'); cat(a)
qstat.info <- function(jobid=NA, var_name=NA, option='Running')
{       
        #
        # MATCH OPTIONS:
        options_vec <- c('Running', 'Queued', 'Finished')
        names(options_vec) <- c('running', 'queued', 'finished')
        option <- pmatch(tolower(option), names(options_vec))
        if(is.na(option)) return(NA)
        option <- unname(options_vec[option])

        # BUILD COMMAND:
        cmd <- 'qstat'
        if(!is.na(jobid)){cmd <- paste0(cmd, '| grep "', jobid,'"')}
        cmd <- paste0(cmd, ' | grep -Eo "',option,':[0-9]+" | grep -Eo [0-9]+')

        # Store variable if wanted:
        if(!is.na(var_name))
        {
                cmd <- paste0(var_name, '=$(' ,cmd, ')')
        }

        cmd <- paste0(cmd, '\n')
        return(cmd)
}
