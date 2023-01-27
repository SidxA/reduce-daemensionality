# so i want packages to read out th eprotophases of every single component in every single file

using Distributed

@everywhere begin

    include("/net/home/lschulz/scripts/toolbox.jl")

    datadir="/net/scratch/lschulz/ssasweep/" # put the /
    outdir = "/net/scratch/lschulz/ssa_sweep_results/"

    function readout_protophases(package,datadir,outdir) #package is a filename
        filename = package
        filenamesplit = split(filename,"_")

        method = filenamesplit[1]
        i = parse(Int64,filenamesplit[2])
        name = filenamesplit[3]
        N = parse(Int64,filenamesplit[4])
        W = parse(Int64,filenamesplit[5])
        k = parse(Int64,filenamesplit[6])
        a = parse(Float64,filenamesplit[7])
        e = parse(Float64,filenamesplit[8])
        t = parse(Float64,filenamesplit[9])

        P = N - W +1
        
        dictfile = load(datadir*filename)
    
        EOF = dictfile["EOF"]
        PC = dictfile["PC"]
        RC = hcat([cutoff(dictfile["RC"][:,i],365) for i=1:k]...)

        EOFproto =  Float32.(hcat([phases(EOF[:,i]) for i in 1:k]...))
        PCproto =   Float32.(hcat([phases(PC[i,:]) for i in 1:k]...))
        RCproto =   Float32.(hcat([phases(RC[:,i]) for i in 1:k]...))

        sigma = Float32.([norm(PC[i,:]) for i=1:k])

        EOF_T1 = Float32.(2*W ./ [iterated_hilbert(EOF[:,i],1)[1] for i=1:k])
        EOF_T4 = Float32.(2*W ./ [iterated_hilbert(EOF[:,i],4)[1] for i=1:k])
        EOF_T8 = Float32.(2*W ./ [iterated_hilbert(EOF[:,i],8)[1] for i=1:k])

        PC_T1 = Float32.(2*P ./ [iterated_hilbert(PC[i,:],1)[1] for i=1:k])
        PC_T4 = Float32.(2*P ./ [iterated_hilbert(PC[i,:],4)[1] for i=1:k])
        PC_T8 = Float32.(2*P ./ [iterated_hilbert(PC[i,:],8)[1] for i=1:k])

        RC_T1 = Float32.(2*N ./ [iterated_hilbert(RC[:,i],1)[1] for i=1:k])
        RC_T4 = Float32.(2*N ./ [iterated_hilbert(RC[:,i],4)[1] for i=1:k])
        RC_T8 = Float32.(2*N ./ [iterated_hilbert(RC[:,i],8)[1] for i=1:k])

        jldsave(outdir*"PROTO$(method)_$(i)_$(name)_$(N)_$(W)_$(k)_$(a)_$(e)_$(t)_$(rand(1000:9999)).jld2",
        EOFproto = EOFproto, PCproto = PCproto, RCproto = RCproto,
        EOF = EOF, PC = PC, RC = RC, sigma = sigma,
        EOF_T1 = EOF_T1, EOF_T4 = EOF_T4, EOF_T8 = EOF_T8,PC_T1 = PC_T1, PC_T4 = PC_T4, PC_T8 = PC_T8,
        RC_T1 = RC_T1, RC_T4 = RC_T4, RC_T8 = RC_T8)

    end

    function packs(datadir)
        fullnames = readdir(datadir)
        nameslist = []
        for i in 1:length(fullnames)
            if fullnames[i][1:5] == "ssasw"
                nameslist=push!(nameslist,i)
            end
        end
        return fullnames[nameslist]
    end
    packages = packs(datadir)
end


Threads.@threads for p=packages
    readout_protophases(p,datadir,outdir)
end

