	//Parameters for the coalescence simulation program : simcoal.exe
	4 samples to simulate :
	//Population effective sizes (number of genes)
	NPOP0
	NPOP1
	NPOP2
	NPOP3
	//Haploid Samples sizes and samples age
	34
	34
	34
	34
	//Growth rates: negative growth implies population expansion
	0
	0
	0
	0
	//Number of migration matrices : 0 implies no migration between demes
	2
	//Migration matrix 0
	0 0 0 M1
	0 0 0 0
	0 0 0 0
	M1 0 0 0
	//Migration matrix 1
	0 0 0 0
	0 0 0 0
	0 0 0 0
	0 0 0 0
	//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
	3 historical event
	TDIV1 2 0 1 1 0 0
	TDIV2 1 3 1 1 0 0
	TDIV3 0 3 1 1 0 1
	//Number of independent loci [chromosome]
	1 0
	//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
	1
	//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
	FREQ 1 0 2.9e-9 OUTEXP
