# dimensionality reduction methods on fluxnet series and further processing

## file contents

- toolbox.jl basic functionality
- struct_blockdata.jl object structures for dimensionality reduction
- iterate_blockdata.jl parallelized dimensionality reduction computations
- phase_comparison.jl further processing and plotting (very messy)


## functionality

- read measurement series of length N
- center measurement series, with observation average and variance
- perform delay embedding with fixed delay embedding parameter W to create data matrix
- centralize data matrix
- chose number k of reduced dimension from k<P, P = N - W +1
- SSA computes k (left) singular vectors of data matrix: modes
- NLSA first samples the diffusion distance distribution for kernel scale parameter computation
- NLSA computes diffusion kernel of data matrix and performs kd-tree to create diffusion distance matrix
- NLSA computes k eigenvectors of diffusion distance matrix: modes
- estimate mode amplitude by variance coverage of original data matrix
- hilbert transform quaternizes individual modes to create instantanious phase in analytic signal representation: protophase
- half protophase zero count gives period length to estimate mode frequency

- showcase characteristics based on the variance and frequency of the individual modes
