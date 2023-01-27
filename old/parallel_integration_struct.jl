using Distributed
include("/net/home/lschulz/scripts/toolbox.jl")
addprocs(50)

@everywhere include("/net/home/lschulz/scripts/toolbox.jl")
#@everywhere include("/net/home/lschulz/scripts/tools.jl")
@everywhere include("/net/home/lschulz/scripts/struct_l_diff.jl")

@everywhere function f(parameters_diffusion)

    #loop the 3 sets
    for (setnumber,setname) in enumerate([dir*"models/harmonic_test",dir*"models/period_double_ra1",dir*"models/period_double_ra2"])
        #load the set and initialize the struct
        lo = local_diffusion(npzread(setname),length(npzread(setname)),1/1,1,400,10,0.5,512,1)
        #loop over the parameters
        i=1
        for W in Array(100:100:6000)
            #set the field and put the rc in a npz
            setfield!(lo,:W,W)
            setfield!(lo,:P,length(npzread(setname))-W+1)
            setfield!(lo,:k,10)
            setfield!(lo,:t,parameters_diffusion[1])
            setfield!(lo,:e,parameters_diffusion[2])
            setfield!(lo,:a,parameters_diffusion[3])
            npzwrite(dir*"runs/"*string("l_diff_",setnumber,"_",parameters_diffusion,"_",i),loop_the_loc_w(lo))
            i+=1
        end
    end
    return nothing
end

@everywhere parameters_diffusion = vcat([(t,e,a,) for t in 1:3, e in 2 .^Array(6:10),a in Array(0.0:0.05:1.0)]...)
#t,e,a,loc_size,repl,W,loc_rate

pmap(W -> f(W), parameters_diffusion)
rmprocs(8)