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

function moving_window_deseason(signal::Vector{Float32},windowlength,overlap)
    #here we iterate the mid part and then the overlp part
    L = Float32[]
    full_windows = Int(floor(length(signal)/(windowlength-overlap)))-2
        #missing first longer clear part, first overlap
        w_1 = 1:windowlength
        w_2 = windowlength-overlap+1:windowlength-overlap+windowlength
        averages_1 = signal[w_1] .- mean(signal[w_1])
        averages_2 = signal[w_2] .- mean(signal[w_2])
        clear = averages_1[1:end-overlap] #no overlap on the first window with zeroth window
        overlaps = mean(hcat(averages_1[end-overlap+1:end],averages_2[1:overlap]),dims=2)[:]
        L = append!(L,clear) #windowlenth-overlap
        L = append!(L,overlaps) #overlap
    #each window starts with the clear middle part and then has the overlap part
    for window_ind = 1:full_windows-1
        w_1 = (windowlength-overlap)*window_ind+1:(windowlength-overlap)*window_ind+windowlength
        w_2 = (windowlength-overlap)*(window_ind+1)+1:(windowlength-overlap)*(window_ind+1)+windowlength
        averages_1 = signal[w_1] .- mean(signal[w_1])
        averages_2 = signal[w_2] .- mean(signal[w_2])
        clear = averages_1[1+overlap:end-overlap]
        overlaps = mean(hcat(averages_1[end-overlap+1:end],averages_2[1:overlap]),dims=2)[:]
        L = append!(L,clear)        #windowlenth-overlap
        L = append!(L,overlaps)     #overlap
    end
        #also missing last clear part, last overlap, end window
        w_1 = (windowlength-overlap)*(full_windows)+1:(windowlength-overlap)*(full_windows+1) +overlap
        w_2 = (windowlength-overlap)*(full_windows+1) +overlap +1:length(signal)
        averages_1 = signal[w_1] .- mean(signal[w_1])
        averages_2 = signal[w_2] .- mean(signal[w_2])
        clear = averages_1[1+overlap:end-overlap]
        overlaps = mean(hcat(averages_1[end-overlap+1:end],averages_2[1:overlap]),dims=2)[:]
        L = append!(L,clear)        #windowlenth-overlap
        L = append!(L,overlaps)     #overlap
        L = append!(L,averages_2)

    return L
end

#we need different data put in functions
#with and without de-seasonalization by lSSA by Miguel      N from 12052 to 11895
#SY SN                                                      methods = Array(1:2), include lSSA_de with signal fix
#function that works on object                              hardcoded N = 11895

function put_data_gSSA(d::matrixholder,signal::Vector{Float32},sampleyear::Int64)
    d.signal = centralizer(deseasonalize_gSSA(centralizer(signal),sampleyear::Int64))
    d.emb = centralized_embed_lag(d.signal,d.W)'
end

function put_data_lSSA(d::matrixholder,signal::Vector{Float32},sampleyear::Int64)
    d.signal = centralizer(deseasonalize_lSSA(centralizer(signal),sampleyear::Int64))
    d.emb = centralized_embed_lag(d.signal,d.W)'
end

function put_data_window(d::matrixholder,signal::Vector{Float32},sampleyear::Int64)
    windowlength=Int(floor(sampleyear/6))
    overlap=Int(floor(sampleyear/12))
    d.signal = centralizer(moving_window_deseason(centralizer(signal),windowlength,overlap))
    d.emb = centralized_embed_lag(d.signal,d.W)'
end

function put_data_raw(d::matrixholder,signal::Vector{Float32})
    d.signal = centralizer(signal)
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

#function that works on object
function put_epsilon(d::matrixholder)
    stepnumber      = 32
    sampling_size   = 64
    it_number       = 4

    d.eps = fit_epsilon(d.emb,stepnumber,sampling_size,it_number,d.W)
end
"""
#functions that fit in the model
"""
#little struct to hold all the variable types for the fit?
struct diffholder
    f::DiffMap{Float32}
    function diffholder(data::Matrix{Float32},k,t,alpha,eps)
        return new(fit(DiffMap,data,maxoutdim=k,t=1, α=alpha, ɛ=eps))
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
    d.lambda = Float32.(diag(1/d.W .* transpose(d.EOF) * cov(d.emb) * d.EOF))
end

"""
#save the results to jld2
"""
function save_results(d::matrixholder,name::String)
    jldsave("$name.jld2",
    EOF = Matrix(d.EOF),PC = Matrix(d.PC),RC = Matrix(d.RC),lambda = d.lambda,eps = d.eps, signal = d.signal)
end

