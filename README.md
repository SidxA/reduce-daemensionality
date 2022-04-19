# plan

- the point is to localize information
- for that we need to split the time series into intervals
- interval selection might be based on brute split, driving period, merging after splitting, 

## ideas

- equal splitting lengthes might make the later comparison of local modes easier
- can we create a structural splitting, by the split and merge with shannon entropy?


- the localized multiscale ssa is able to detect a frequency shift aka devilsstaircase in the enso data

- there is the possibility of 'sifting' the signal to get rid of riding waves and to create local mean zero
- there is the possibility of smooth orthogonal projection (bell)
- there is the possibility of merging short intervals in order to avoid redundancy

# computation structure

## 1 mode extraction

- fit the data into a shape that is most convenient for the extraction approaches
- PCA in MultiVariateStats.jl uses general matrix input
- DiffusionMap in ManifoldLearning.jl uses general matrix input
- output storage: log everything in individual run folders!
	
- question of window length: lets make a trial with the PCA, then check the SSA
	

## 2 differentiating properties

## 3 coupling model

# done

- logging program <create_log_structure.jl> with init <init_logging()> returns string directory variable for logging appendix
- testfile are accessed by <<DataFrame(CSV.File(path, header=1, delim=","))>>
