using Distributed

#addprocs(8)

@everywhere begin
    include("/net/home/lschulz/scripts/toolbox.jl")

    function parameters()
        parameterset = []
        i=1
        for epsilon = [0.0,0.1,0.5]
            for W = [20,40,100,200,400]

                p = 0
                w = 2.5 * 10^-1
                parameterset = push!(parameterset,[i,p,epsilon,w,W])
                i+=1

                p = 1
                w = 2.5 * 10^-5
                parameterset = push!(parameterset,[i,p,epsilon,w,W])
                i+=1
                p = 2
                w = 2.5 * 10^-9
                parameterset = push!(parameterset,[i,p,epsilon,w,W])
                i+=1
            end
        end
        println("total parameter combinations $i")
        return parameterset
    end

    function  pipeline(p,epsilon,w,W)
        #this runs by model
        #p = 0
        #epsilon = 0.01
        #w1 = 0.025
        #w2 = 0.025
        w1 = w2 = w
        N = 5
    
        #then performs the decomposition of the individual signals
        k = 3
        #W = 400
    
        #signal boundaries
        N_bound1 = 1
        Bool_centre1 = true
        N_bound2 = 40
        Bool_centre2 = false
    
        #result: table with signal,gSSA,gDiff protophases and table with variances
        # IS PUT TO MODEL T
        phases = Array{Float64}(undef,N,8000)
        variances = Array{Float64}(undef,(2*k+1)*4,N)
        protophases = Array{Float64}(undef,(2*k+1)*4,N,8000)
    
        #model
        model = harmonic_model(p,epsilon,w1,w2,N) #Tx2N
        #now per signal of model
        for signal_index = 1:N
            phases[signal_index,:]      = model[:,signal_index]
            signal                      = model[:,N+signal_index]
            
    
        #individually: cenralized and not, boundary and not...
            components = hcat(
                components_g(signal,W,k,N_bound1,Bool_centre1),
                components_g(signal,W,k,N_bound1,Bool_centre2),
                components_g(signal,W,k,N_bound2,Bool_centre1),
                components_g(signal,W,k,N_bound2,Bool_centre2))
    
        #protophase per component and save the variances
            for c_index in 1:size(components)[2]
                HilbertT = hilbert_transform(components[:,c_index])
                protophases[c_index,signal_index,:] = atan.(imag(HilbertT),real(HilbertT))
                variances[c_index,signal_index] = var(components[:,c_index])
            end
        end
    
        return phases,variances,protophases
    end

    function put_the_pipe(para)
        i = Int.(para[1])

        println("start $i")


        p = para[2]
        epsilon = para[3]
        w = para[4]
        W = Int.(para[5])
        phases,variances,protophases = pipeline(p,epsilon,w,W)
        npzwrite("/net/home/lschulz/logs/KW_25/runs/model/"*"$i"*"_phases",phases)
        npzwrite("/net/home/lschulz/logs/KW_25/runs/model/"*"$i"*"_protophases",protophases)
        npzwrite("/net/home/lschulz/logs/KW_25/runs/model/"*"$i"*"_variances",variances)
        return nothing
    end

    parameterset = parameters()
end


#pmap(W -> put_the_pipe(W), parameterset)

#rmprocs(8)

function crude_normalization(protophase::Vector)
    for (i,x) in enumerate(protophase)
        if x >= 2pi
            protophase[i:end] .-= 2*pi
        end
    end
    return protophase
end


function phasemodel(p,epsilon,w,phi)
    N=5
    x = 1:8000
    fL = Array{Float64}(undef,0)
    pL = Array{Float64}(undef,0)
    for n in 1:N
        wT(x,w) = w*x^p
        eps_add = (rand(8000) .- 1/2 ) * epsilon
        eps_mul = 1 .+ (rand(8000) .- 1/2 ) * epsilon
        phaseshift = rand()*2*pi
        wL = vcat([wT(x,w) for x in 1:4000],[wT(x,w) for x in 4001:8000])
        pL = append!(pL,wL.*x.+phi.+ eps_add)
        fL = append!(fL,@. sin(x*wL+phi.+ eps_add)*eps_mul)
    end
    return hcat(reshape(pL,8000,N),reshape(fL,8000,N))
end

model = hmodel(0,0.1,0.006,0)
p0 = hcat([crude_normalization(model[:,i]) for i in 1:5]...)
signal = model[:,i]
c = components_g(signal,100,3,100,true)
ssa = sum(components[:,2:3],dims=2)
diff = sum(components[:,5:6],dims=2)
