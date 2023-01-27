# dimensionality reduction methods on fluxnet series and further processing

## file contents

- toolbox.jl basic functionality
- struct_blockdata.jl object structures for dimensionality reduction
- iterate_blockdata.jl parallelized dimensionality reduction computations
- phase_comparison.jl further processing and plotting (very messy)


## functionality

- read measurement series of length **N**
- center measurement series, with observation average and variance
- perform delay embedding with fixed delay embedding parameter W to create data matrix
- centralize data matrix
- chose number **k** of reduced dimension from **k<P, P = N - W +1**
- SSA computes **k** (left) singular vectors of data matrix: modes
- NLSA first samples the diffusion distance distribution for kernel scale parameter **e** computation
- NLSA computes diffusion kernel of data matrix and performs kd-tree to create diffusion distance matrix
- NLSA computes k eigenvectors of diffusion distance matrix: modes
- estimate mode amplitude by variance coverage of original data matrix
- create reconstructed time series from modes
- hilbert transform quaternizes individual mode to create instantanious phase in analytic signal representation: protophase
- half protophase zero count gives period length to estimate mode frequency
- showcase characteristics based on the variance and frequency of the individual modes

## purpose

- dimensionality reduction methods can do an additive decomposition of time series
- this decomposition is datadriven and orthogonal, which poses the question wether it can separate different timescales in time series
- since modes are quasiperiodic individual timescales can be estimated by frequencies
- dimensionality reduction suffers from artifacts linked to orthogonality: variance compression & degeneracy
- how do linear (SSA, keeps global metric) and nonlinear (NLSA,keeps local metric - diffusion distance) differ in attributing timescales?
- how do these artifacts play out?
- how does the embedding length parameter influence the attributed timescales?
