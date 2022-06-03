//Parameters for the coalescence simulation program : fsimcoal2.exe
5 samples to simulate :
//Population effective sizes (number of genes)
NPOP0
NPOP1
NPOP2
NPOP3
NPOPG
//Samples sizes and samples age
8
8
24
14
0
//Growth rates  : negative growth implies population expansion
0
0
0
0
0
//Number of migration matrices : 0 implies no migration between demes
0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
6 historical event
TADM1 2 4 0.13 1 0 0
TDIV1 2 3 1 RES1 0 0
TADM2 3 0 0.03 1 0 0
TDIVG 4 3 1 RESG 0 0
TDIV2 1 3 1 RES2 0 0
TDIV3 0 3 1 RES3 0 0
//Number of independent loci [chromosome]
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ  1   0   1.49e-8 OUTEXP
