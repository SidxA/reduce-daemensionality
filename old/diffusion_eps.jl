function epsi(data::Matrix{Float32},stepnumber,sampling_size)

    """
    function Weighting(Weights,eps)
        for i in 1:P, j in 1:P
            Weights[i,j] = exp(- norm(data[:,i] - data[:,j])^2 / eps)
        end
    end
    """
    #L(eps::Float32) = sum(Weights)

    data = data[:,rand(1:size(data)[2],sampling_size)]
    P = sampling_size


    #min_eps = Float32(1/P * sum([min(deleteat!([norm(data[:,i]-data[:,j])^2 for j=1:P],i)...) for i=1:P]))
    #max_eps = Float32(1/P * sum([max(deleteat!([norm(data[:,i]-data[:,j])^2 for j=1:P],i)...) for i=1:P]))

    min_eps = Float32(10^0)
    max_eps = Float32(10^7)

    eps_list = Vector{Float32}(undef,stepnumber)
    eps_list = 10 .^ range(log10(min_eps),log10(max_eps),length=stepnumber)

    L_list = Vector{Float32}(undef,stepnumber)
    Weight = Array{Float32}(undef,P,P)
    for (i,eps) in enumerate(eps_list)
        #Weighting(Weights,eps)
        for i in 1:P, j in 1:P
            Weight[i,j] = exp(- norm(data[:,i] - data[:,j])^2 / eps)
        end
        L_list[i] = sum(Weight)
    end

        p0 = ones(3)

        #the atan and the htan fit
        #maybe investigate this later
        #beware eps2 is negative for some reason

        # atan is much more stable lets only use it!

        model1(eps,p) = p[1] .* tanh.(eps .- p[2]) .+ p[3]
        eps1 = coef(curve_fit(model1, eps_list, log10.(L_list), p0))[2]

        model2(eps,p) = p[1] .* atan.(eps .- p[2]) .+ p[3]
        eps2 = coef(curve_fit(model2, eps_list, log10.(L_list), p0))[2]
    return Int.(floor.((eps1,eps2)))
end

