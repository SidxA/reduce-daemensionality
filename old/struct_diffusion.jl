
#suppose that we read in a single series that has been already preprocessed

#the idea is again to use time series of identical length to fit in the new series in the existing structure to save ram
#fixed N,W,k
#fixed alpha,t

#should be one struct initialized with N W k

    mutable struct diffusion
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
        function diffusion(N,W,k)
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

#the data put

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

    #function that works on object
    function put_data(d::diffusion,signal::Vector{Float32})
        d.signal = signal
        d.emb = centralized_embed_lag(d.signal,d.W)'
    end


#the epsilon findout

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
        for t in 1:it_number
            object.data_samples = data[:,rand(1:P,sampling_size)]

            for (i,eps) in enumerate(object.eps)
                for i in 1:sampling_size, j in 1:sampling_size
                    object.Weight[i,j] = exp(- norm(data[:,i] - data[:,j])^2 / eps)
                end
                object.L[i] = sum(object.Weight)
            end

            p0 = ones(3)
            model(eps,p) = p[1] .* atan.(eps .- p[2]) .+ p[3]
            fit_eps_L[t] = coef(curve_fit(model, object.eps, log10.(object.L), p0))[2]
        end
        return - Int64(floor(median(fit_eps_L)))
    end

    #function that works on object                      # FIXED SAMPLING,STEPS,ITER
    function put_epsilon(d::diffusion)
        stepnumber      = 16
        sampling_size   = 128
        it_number       = 8

        d.eps = fit_epsilon(d.emb,stepnumber,sampling_size,it_number,d.W)
    end

#function that fits in diffmodel

    #little struct to hold all the variable types for the fit?
    mutable struct holder
        f::DiffMap{Float32}
        function holder(data::Matrix{Float32},k,t,alpha,eps)
            return new(fit(DiffMap,data,maxoutdim=k,t=t, α=alpha, ɛ=eps))
        end
    end

    function gettheproj(diff::holder)
        return Matrix(diff.f.proj')::Matrix{Float32},diff.f.λ::Vector{Float32}
    end

    #function that works on the object                  # FIXED t=1 alpha = 1
    function put_EOF(d::diffusion)
        N = d.N
        W = d.W
        P = N - W +1
        k = d.k
        t = 1
        alpha = 1.0
        eps = d.eps

        diff = holder(d.emb,k,t,alpha,eps)
        d.EOF,d.lambda = gettheproj(diff)

    end

#building PC,RC                      
    function calculate(d::diffusion)
        d.PC = hcat([pc!(d.signal,d.EOF[:,i],d.N-d.W-1,d.W) for i in 1:d.k]...) # d.emb * d.EOF
        d.RC = hcat([reconstructor(d.PC[:,i],d.EOF[:,i],d.N,d.W) for i in 1:d.k]...)
    end

#combine it all
    function doitall(signal::Vector{Float32})
        W = Int(12*365.25)
        k = 64
        d = diffusion(length(signal),W,k)
        put_data(d,signal)
        put_epsilon(d)
        put_EOF(d)
        calculate(d)
        return d
    end
