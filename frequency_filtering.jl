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
spot = 1
method="ssa"
vari = 1
W = Int(floor(7*365.25))
name = "$(spot)_$(method)_$(W)_$(vari)"


data = load("/net/scratch/lschulz/modes/"*"name")
lambda = data["lambda"]
signal = data["signal"]
modes = data["modes"]
yearsamples = 365
savedirname = dir*"test"

using DSP
using CairoMakie

function protophase(signal::Vector{Float32})
    ht = hilbert(Float64.(signal))
    return Float32.(atan.(imag(ht), real(ht)))
end

count_sign_flip(signal) = sum([(sign.(signal[i])!=sign.(signal[i+1])) ? 1 : 0 for i in 1:length(signal)-1])
protoperiod(protophase) = yearsamples / count_sign_flip(protophase) / 2



rec_protophase(RC::Matrix{Float32},Nit) = Float32.(hcat([protophase(RC[:,i],Nit) for i=1:size(RC)[2]]...))

Nit = 40
signal = Filius["EOF"][:,1]
for i=1:Nit-1
    phi = protophase(signal)
    T = protoperiod(phi)
    println("$i \t $T")
    signal = phi
end


function plot_signal(signal)
    # Create a scatter plot of the signal
    F,ax,s = lines(signal)
    save(savedirname*".png",F)
end

protophases = rec_protophase(RC,1)
protofreq = hcat([protofrequency(EOF[:,kk],yearsamples,1) for kk in 1:k]...)
