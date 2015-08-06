# GetBaseCounts 1.4.0

#### Calculate the read depth and allele depth in the a BAM file for all the sites in a given VCF file


Usage:

[REQUIRED ARGUMENTS]

	--fasta                 <string>                        Input reference sequence file
	--bam                   <string>                        Input bam file
	--vcf                   <string>                        Input vcf file, it needs to be sorted within each chromosome to optimize the running speed, gzipped vcf file is supported
	--output                <string>                        Output file

[OPTIONAL ARGUMENTS]

	--thread                <int>                           Number of thread. Default 1
	--sort_output                                           Sort output file by genomic position, this option requires addtional memory
	--compress_output                                       Compress the output and write gzipped file directly
	--maq                   <int>                           Mapping quality threshold. Default 15
	--baq                   <int>                           Base quality threshold, Default 20
	--cov                   <int>                           Minimum coverage applied to BASEQ_depth. Default 0
	--filter_duplicate      [0, 1]                          Whether to filter reads that are marked as duplicate. 0=off, 1=on. Default 1
	--filter_improper_pair  [0, 1]                          Whether to filter reads that are marked as improperly paired. 0=off, 1=on. Default 1
	--filter_qc_failed      [0, 1]                          Whether to filter reads that are marked as failed quality control. 0=off, 1=on. Default 0
	--filter_indel          [0, 1]                          Whether to filter reads that contain indels. 0=off, 1=on. Default 0
	--filter_non_primary    [0, 1]                          Whether to filter reads that are marked as non primary alignment. Default 1
	--suppress_warning      <int>                           Only print a limit number of warnings for each type. Default 3
	--help                                                  Print command line usage


[ADVANCED ARGUMENTS, CHANGING THESE ARGUMENTS MAY SIGNIFICANTLY AFFECT MEMORY USAGE AND RUNNING TIME. USE WITH CAUTION]

	--max_block_size        <int>                           The maximum number of vcf entries that can be processed at once per thread. Default 10000
	--max_block_dist        <int>                           The longest spanning region (bp) of vcf chunks that can be processed at once per thread. Default 1000000
