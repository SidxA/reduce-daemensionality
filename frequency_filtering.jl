"""
code execution: extract EOF,lambda,signal from calculated files and jldsave in scratch
"""

W = Int(floor(7*365.25))

for method=["ssa","diff"], vari=1:2, spot=1:20
    read_modes_single_file(method,W,vari,spot)
end

"""
mode example
"""
using FileIO
using JLD2

for spot = 1:5
    method = "ssa"
    #method="diff"
    vari = 1
    #W = Int(floor(7*365.25))
    #name = "$(spot)_$(method)_$(W)_$(vari)"
    #outdir="/net/scratch/lschulz/fluxfullset/"
    yearsamples=365
    preproc="raw."
    filename = create_file_list(outdir,method,W,vari,preproc)[spot]

    file = load(filename)
    eof = file["EOF"]
    lambda = file["lambda"]
    pc = file["PC"]
    rc = file["RC"]
    signal = file["signal"]
    N = length(signal)

    harmonics_eof = hcat([eof[:,l] .+ 0.5*i for (i,l) in enumerate(1:6)]...)
    F,ax,s = series(1:W,harmonics_eof',solid_color="black")
    save(savedirname*"ssa_$(spot)_eof.png",F)

    harmonics_rc = hcat([rc[:,l] .+ i for (i,l) in enumerate(1:6)]...)
    F,ax,s = series(1:N,harmonics_rc',solid_color="black")
    save(savedirname*"ssa_$(spot)_rc.png",F)

    method="diff"

    filename = create_file_list(outdir,method,W,vari,preproc)[spot]
    file = load(filename)
    eof = file["EOF"]
    lambda = file["lambda"]
    pc = file["PC"]
    rc = file["RC"]
    signal = file["signal"]


    harmonics_eof = hcat([eof[:,l] .+ 0.5*i for (i,l) in enumerate(1:6)]...)
    F,ax,s = series(1:W,harmonics_eof',solid_color="black")
    save(savedirname*"nlsa_$(spot)_eof.png",F)

    harmonics_rc = hcat([rc[:,l] .+ i for (i,l) in enumerate(1:6)]...)
    F,ax,s = series(1:N,harmonics_rc',solid_color="black")
    save(savedirname*"nlsa_$(spot)_rc.png",F)


end


"""
    L = []
    for kap = 1:48
        zeros = count_sign_flip(eof[:,kap])
        if zeros%2 == 0
            #println("$kap \t $zeros \t $(zeros/14)")
            L = append!(L,kap)
        else
            #println("$kap \t $zeros")
        end
    end
"""



using DSP
using CairoMakie
savedirname = dir*"test"

function protophase(signal::Vector{Float32})
    ht = hilbert(Float64.(signal))
    return Float32.(atan.(imag(ht), real(ht)))
end

count_sign_flip(signal) = sum([(sign.(signal[i])!=sign.(signal[i+1])) ? 1 : 0 for i in 1:length(signal)-1])
protoperiod(protophase) = yearsamples / count_sign_flip(protophase) / 2



rec_protophase(RC::Matrix{Float32},Nit) = Float32.(hcat([protophase(RC[:,i],Nit) for i=1:size(RC)[2]]...))

Nit = 40
signal = mode
for i=1:Nit-1
    phi = protophase(signal)
    T = count_sign_flip(phi)#protoperiod(phi)
    println("$i \t $T")
    signal = phi
end


function plot_signal(signal)
    # Create a scatter plot of the signal
    F,ax,s = lines(signal)
    save(savedirname*".png",F)
end



