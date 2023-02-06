include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")

savedirname=dir*"signals/"
outdir="/net/scratch/lschulz/fluxfullset_1a/"
yearsamples=365
preproc="raw."
W = 365

"""
use the hard mode fuilter from 02 02 to create a larger combined figure for each spot
"""


savedirname=dir*"hard_filter/"
vari = 1
n_modes = 8
Nit = 3

using DSP
using CairoMakie

function protophase(signal::Vector{Float32})
    ht = hilbert(Float64.(signal))
    return Float32.(atan.(imag(ht), real(ht)))
end

count_sign_flip(signal) = sum([(sign.(signal[i])!=sign.(signal[i+1])) ? 1 : 0 for i in 1:length(signal)-1])


function hard_mode_filter(one,two,three)

    if one==two==three==2
        return true,2
    elseif one==two==three==3
        return true,3
    elseif one==two==three==4
        return true,4
    elseif one==two==three==6
        return true,6
    elseif one==two==three==8
        return true,8
    else
        return false,0
    end
end



for spot = [3,9,14,15,16]
    F = Figure()

    for (i,method) = enumerate(["ssa","diff"])
        a = Axis(F[i,1],title=method)
        filename = create_file_list(outdir,method,W,vari,preproc)[spot]

        file = load(filename)
        eof = file["EOF"]
        lambda = file["lambda"]

        indices = sortperm(lambda,rev=true)[1:n_modes]
        L_modes = Int64[]
        L_periods = Int64[]
        L_out = Int64[]
        for kap in 1:n_modes

            signal = eof[:,indices[kap]]
            L = Int64[]
            for i=1:Nit
                phi = protophase(signal)
                T = count_sign_flip(phi)
                L = append!(L,Int64(T))

                signal = phi
            end
            fil = hard_mode_filter(L...)
            if fil[1]
                L_modes = append!(L_modes,indices[kap])
                L_periods = append!(L_periods,L[1])
            else
                L_out = append!(L_out,indices[kap])
            end
        end

        #plot the results
        signal = file["signal"]
        rc = file["RC"]
        trend = 0 .+ sum(rc[:,L_modes],dims=2)[:]
        residual = signal .- trend
        combined = hcat(signal,trend .+3 ,residual .+ 6)
        s = series!(a,1:length(signal),combined')
        
        l = file["lambda"]
        l_sum = round(sum(l[L_modes]),digits=3)
        mode_periods = Int.(L_periods ./2)

        a2 = Axis(F[i,2],title = "var $l_sum \t periods $mode_periods")
        modes_out = 0.0 .+ eof[:,L_out]
        modes_in = 0.2 .+ eof[:,L_modes]
        s1 = series!(a2,1:W,modes_out',solid_color="black")

        s2 = series!(a2,1:W,modes_in',solid_color="green")
        #t = text!(1,1,text="$l_sum \t $L_periods")
    end
    save(savedirname*"_$(spot)_trend.png",F)
end
