# structure for the cluster: todo

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
