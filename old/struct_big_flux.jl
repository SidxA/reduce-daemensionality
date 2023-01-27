"""
the big fluxnet run with boundary
signal -> boundary -> centralizer ( -> seasons -> centralizer)

one struct per kernel

fit SSA and diff in after one another
save every single one individually

the function only needs the signal and the saving information

so the packaging provides one signal per package
need location and variable indexing system that suits the analysis
yes the indices are a good way

test how long one package takes
"""

using Distributed
@everywhere begin


    include("/net/home/lschulz/scripts/toolbox.jl")
    include("/net/home/lschulz/scripts/lSSA_deseasonalization.jl")

    
#choice of W: 16 * 365.25 = 5844


#suppose that we read in a single series that has been already preprocessed

#the idea is again to use time series of identical length to fit in the new series in the existing structure to save ram
#fixed N,W,k
#fixed alpha,t
"""
should be one struct initialized with N W k
"""
mutable struct matrixholder
    #requirements for function that puts in data
    N           ::Int64
    W           ::Int64
    k           ::Int64
    signal      ::Vector{Float32}
    emb         ::Matrix{Float32}

    #requirements for function that finds out epsilon
    eps         ::Int64

    #requirements for dimensionality reduction
    lambda      ::Vector{Float32}
    EOF         ::Matrix{Float32}
    PC         ::Matrix{Float32}
    RC         ::Matrix{Float32}

    #the init function
    function matrixholder(N,W,k)
        signal      = Vector{Float32}(undef,N)
        emb         = Matrix{Float32}(undef,N-W+1,W)
        eps         = 0
        EOF         = Matrix{Float32}(undef,W,k)
        lambda      = Vector{Float32}(undef,k)
        PC         = Matrix{Float32}(undef,N-W+1,k)
        RC         = Matrix{Float32}(undef,N,k)
        return new(
            N,
            W,
            k,
            signal,
            emb,
            eps,
            lambda,
            EOF,
            PC,
            RC
        )
    end
end
"""
#the data put
"""
#include("/net/home/lschulz/scripts/lSSA_deseasonalization.jl")
#we need different data put in functions
#with and without de-seasonalization by lSSA by Miguel      N from 12052 to 11895
#SY SN                                                      methods = Array(1:2), include lSSA_de with signal fix
    #function that works on object                          hardcoded N = 11895
    function put_data_by_season(d::matrixholder,raw_signal::Vector{Float32},season::Int64)

        if season == 1
            signal = centralizer(deseasonalize_lSSA(centralizer(boundary(raw_signal,365)))[1:d.N])
        elseif season == 2
            signal = centralizer(deseasonalize_gSSA(centralizer(boundary(raw_signal,365)))[1:d.N])
        elseif season == 3
            signal = centralizer(boundary(raw_signal,365)[1:d.N])
        end
        d.signal = signal
        d.emb = centralized_embed_lag(d.signal,d.W)'
    end


function centralized_embed_lag(data::Vector{Float32},W::Int64)
    #centralize
    function centralizer(data)
        m = mean(data)
        cv = std(data)
        return  (data.-m)./cv
    end
    Y = []
    for i=1:length(data)-W+1
        Y = append!(Y,centralizer(data[i:i+W-1]))
    end
    return reshape(float.(Y),W,length(data)-W+1)::Matrix{Float32}
end

"""
#the epsilon findout
"""
#required struct to iterate with different epsilon
mutable struct epsilon
    data_samples :: Matrix{Float32}
    eps         :: Vector{Float32}
    L           :: Vector{Float32}
    Weight      :: Matrix{Float32}
    function epsilon(W,stepnumber,sampling_size)
        min_eps = Float32(10^0)
        max_eps = Float32(10^7)

        eps             = 10 .^ range(log10(min_eps),log10(max_eps),length=stepnumber)
        data_samples    = Matrix{Float32}(undef,W,sampling_size)
        L               = Vector{Float32}(undef,length(eps))
        Weight          = Matrix{Float32}(undef,sampling_size,sampling_size)

        return new(data_samples,eps,L,Weight)
    end
end

# iteration function that returns the median
function fit_epsilon(data::Matrix{Float32},stepnumber,sampling_size,it_number,W)
    P = size(data)[2]
    object = epsilon(W,stepnumber,sampling_size)
    fit_eps_L = Vector{Float32}(undef,it_number)
    sample = rand(1:P,sampling_size,it_number)
    for t in 1:it_number
        object.data_samples = data[:,sample[:,t]]

        for (i,eps) in enumerate(object.eps)
            for i in 1:sampling_size, j in 1:sampling_size
                object.Weight[i,j] = exp(- norm(object.data_samples[:,i] - object.data_samples[:,j])^2 / eps)
            end
            object.L[i] = sum(object.Weight)
        end

        p0 = ones(3)
        model(eps,p) = p[1] .* atan.(eps .- p[2]) .+ p[3]
        p_midpoint = coef(curve_fit(model, log10.(object.eps), log10.(object.L), p0))
        a = p_midpoint[2]
        b = median(log10.(object.L[1:8]))
        one_over_e = model(b+(a-b)/exp(1),p_midpoint)
        fit_eps_L[t] = 10^(one_over_e)
    end
    return Int64(floor(median(fit_eps_L)))
end

#function that works on object                      # FIXED SAMPLING,STEPS,ITER
function put_epsilon(d::matrixholder)
    stepnumber      = 16
    sampling_size   = 128
    it_number       = 8

    d.eps = fit_epsilon(d.emb,stepnumber,sampling_size,it_number,d.W)
end
"""
#functions that fit in the model
"""
#little struct to hold all the variable types for the fit?
struct diffholder
    f::DiffMap{Float32}
    function diffholder(data::Matrix{Float32},k,t,alpha,eps)
        return new(fit(DiffMap,data,maxoutdim=k,t=t, α=alpha, ɛ=eps))
    end
end

#type-stable extraction
function gettheproj_diff(diff::diffholder)
    return Matrix(diff.f.proj')::Matrix{Float32},diff.f.λ::Vector{Float32}
end

#function that works on the object                  # FIXED t=1 alpha = 1
function put_EOF_diff(d::matrixholder)

    t = 1
    alpha = 1.0
    eps = d.eps

    diff = diffholder(d.emb,d.k,t,alpha,eps)
    d.EOF,d.lambda = gettheproj_diff(diff)

end

struct ssaholder
    f::PCA{Float32}
    function ssaholder(data::Matrix{Float32},k)
        return new(fit(PCA,data',maxoutdim=k,method=:auto))
    end
end

#type-stable extraction
function gettheproj_ssa(ssa::ssaholder)
    return Matrix(ssa.f.proj)::Matrix{Float32},(ssa.f.prinvars ./ ssa.f.prinvars[1])::Vector{Float32}
end

#function that works on the object
function put_EOF_ssa(d::matrixholder)
    ssa = ssaholder(d.emb,d.k)
    d.EOF,d.lambda = gettheproj_ssa(ssa)

end

"""
#building PC,RC   
"""                  
function calculate(d::matrixholder)
    d.PC = d.emb * d.EOF # hcat([pc!(d.signal,d.EOF[:,i],d.N-d.W-1,d.W) for i in 1:d.k]...) # 
    d.RC = hcat([reconstructor(d.PC[:,i],d.EOF[:,i],d.N,d.W) for i in 1:d.k]...)
end

"""
#save the results to jld2
"""
function save_results(d::matrixholder,method::String,outdir::String,loc_ind::Int64,var_ind::Int64,season::Int64)

    jldsave(outdir*method*"_$(loc_ind)_$(var_ind)_$(season).jld2",
    EOF = Matrix(d.EOF),PC = Matrix(d.PC),RC = Matrix(d.RC),lambda = d.lambda,eps = d.eps, signal = d.signal)
end

"""
#combine it all : packaging system
"""
# outdir loc_ind var_ind
struct package
    outdir::String
    loc_ind::Int64
    var_ind::Int64
    function package(outdir,loc_ind,var_ind)
        return new(outdir,loc_ind,var_ind)
    end
end

function doitall(d::matrixholder,wholedata::Array{Float32},p::package)
    outdir= p.outdir
    loc_ind= p.loc_ind
    var_ind= p.var_ind
    for season in 1:3
        put_data_by_season(d,wholedata[:,loc_ind,var_ind],season)
        put_epsilon(d)
        put_EOF_diff(d)
        calculate(d)
        save_results(d,"diff",outdir,loc_ind,var_ind,season)
        put_EOF_ssa(d)
        calculate(d)
        save_results(d,"ssa",outdir,loc_ind,var_ind,season)
    end
end

"""
# the big put together
"""

    
    N = 11895
    W = 5844
    k = 48

    d = matrixholder(N,W,k)

    N_raw = 11322
    spots=71
    varlist = ["TA_ERA","VPD_ERA","PA_ERA"]
    wholedata = reshape(hcat([Float32.(load("/net/scratch/lschulz/FLUXERAI/F32_$(flux).jld2")["data"]) for flux in varlist]...),(N_raw,spots,length(varlist)))
    outdir = "data/fluxfull/"

function packaging(outdir)
    P = package[]
    for i in 1:spots,j in 1:length(varlist)
        P = push!(P,package(outdir,i,j))
    end
    return P
end

function packaging2(outdir)
    L = package[]
    for loc_ind in 1:71, var_ind in 1:3,j=1:3
        name = "ssa"*"_$(loc_ind)_$(var_ind)_$j.jld2"
        if isempty(findall(x->x==name,readdir("data/fluxfull")))
            L = push!(L,package(outdir,loc_ind,var_ind))
        end
    end
    return L
end
end

Threads.@threads for p=packaging2(outdir)
    doitall(d,wholedata,p)
    end