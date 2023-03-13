
include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")
using CairoMakie


"""
plot the seasonality behavior...
"""
function gauss(x,p) # cauchy peak with two parameters
    x0 = p[1]
    gamma = p[2]
    sigma = p[3]
    #1/(sigma * sqrt(2*pi))
    @. return gamma * exp.(-(x-x0)^2/sigma^2/2)
end

function fit_gauss(yvalues,xvalues) # fit the cauchy peak to a vec of numbers
    onebins = xvalues
    bins = yvalues
    p0 = ones(3) .* 5
    return coef(curve_fit(gauss,bins,onebins,p0))
end


function mode_harmonicity_estimation_gauss(Eof,threshold)
    li_harmonics = Int64[]
    li_mixed = Int64[]
    li_h_freq = Float64[]
    li_m_freq = Float64[]
    li_residual=Int64[]
    W,k = size(Eof)

    for i in 1:k

        mode = Eof[:,i]

        Four = fft(mode) |> fftshift
        freqs = fftfreq(length(tw), 1.0/Ts) |> fftshift

        freqstart = findall(x->x>=1/12,freqs)[1]
        freqend = findall(x->x>=6,freqs)[1]

        
        val_domain = abs.(Four[freqstart:freqend])
        #val_domain ./= maximum(val_domain)

        freq_domain = freqs[freqstart:freqend]

        freq, value,sigma = fit_gauss(freq_domain,val_domain)
        peak = gauss(freq_domain,(freq,value,sigma))
        residual = val_domain .- peak

        freq_round = round(freq,digits=1)
        
        val = abs(1/sigma /sqrt(2*pi))
        
        #println("$i \t $freq \t $val \t $(maximum(residual))")

        #if maximal residual value other then single peak is larger then relative peak threshold value, its 
        if maximum(residual .+ 0.0)*threshold <= value &&  freq_round in Float64.(1:10)
            li_harmonics = append!(li_harmonics,i)
            li_h_freq = append!(li_h_freq,freq)
        elseif maximum(residual .+ 0.0)*threshold >= value && freq_round in Float64.(1:10)
            li_mixed = append!(li_mixed,i)
            li_m_freq = append!(li_m_freq,freq)
        elseif maximum(residual .+ 0.0)*threshold >= value
            li_residual = append!(li_residual,i)
        else
            #println("no peak")
        end

    end

    return li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual
end


# vielsalm: 3


threshold = 10
savedir=dir*"rec/"
variable_list = ["GPP","NEE","RECO","soil water content","soil temperature"]
meta = load("/net/home/lschulz/logs/KW_2_06/meta.jld2")["meta"]
outdir="/net/scratch/lschulz/fluxfullset/"

W = 2556
preproc = "raw."
kappa= 48
yearsamples=365.25
vari = 1
# Number of points
N= 5478 
# Sample period
Ts = 1 / 365
# Start time 
t0 = 0 
tmax = t0 + (N-1) * Ts
# time coordinate
t = t0:Ts:tmax
0:Ts:(t0+(W-1)*Ts)
tw = t0:Ts:(t0+(W-1)*Ts)


savedir = dir#*"poster/"


lw = 3
lw_s = 1
ms = 4
col_fundamental = "orangered"
col_harmonic = "blue"
col_signal = "grey50"
col_ssa = "darkgreen"
col_nlsa = "purple"
offset = 2.5


for spot in [1]

    threshold = 10
    vari=1

    #ssa

    method="ssa"
    Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
    file = load(Filename)
    signal = file["signal"]
    lambda = file["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file["EOF"][:,indices]
    PC = file["PC"][:,indices]
    RC = file["RC"][:,indices]
    lambda = lambda[indices]

    # ssa reconstruction: freq 1

    li_harmonics=[1,2]
    val = (round(sum(lambda[li_harmonics]),digits=3))
    le = length(li_harmonics)
    #val_m = (round(sum(lambda[li_mixed]),digits=3))
    #le_m = length(li_mixed)
    ssa_trend_first = sum(RC[:,li_harmonics],dims=2)[:]

    # ssa reconstruction: freq harmonics

    li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)
    val = (round(sum(lambda[li_harmonics]),digits=3))
    le = length(li_harmonics)
    #val_m = (round(sum(lambda[li_mixed]),digits=3))
    #le_m = length(li_mixed)
    ssa_trend_harm = sum(RC[:,li_harmonics],dims=2)[:]

    #nlsa 

    method="diff"
    Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
    file = load(Filename)
    signal = file["signal"]
    lambda = file["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file["EOF"][:,indices]
    PC = file["PC"][:,indices]
    RC = file["RC"][:,indices]
    lambda = lambda[indices]

    # nlsa reconstruction: freq 1

    li_harmonics=[1,2]
    val = (round(sum(lambda[li_harmonics]),digits=3))
    le = length(li_harmonics)
    #val_m = (round(sum(lambda[li_mixed]),digits=3))
    #le_m = length(li_mixed)
    nlsa_trend_first = sum(RC[:,li_harmonics],dims=2)[:]

    # nlsa reconstruction: freq harmonics

    li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)
    val = (round(sum(lambda[li_harmonics]),digits=3))
    le = length(li_harmonics)
    #val_m = (round(sum(lambda[li_mixed]),digits=3))
    #le_m = length(li_mixed)
    nlsa_trend_harm = sum(RC[:,li_harmonics],dims=2)[:]

    trajectorymatrix = hcat(
    ssa_trend_first .+ 0.0,
    ssa_trend_harm .+ 2.0,
    nlsa_trend_first.+ 4.0,
    nlsa_trend_harm .+ 6.0)

    """
    plotting procedure : fundamental versus harmonic
    """

    #signals

    years = 2005 .+ ((1:5478) ./ 365)
    fundamental = ssa_trend_first
    harmonic = ssa_trend_harm
    res_fund = signal .- fundamental
    res_harm = signal .- harmonic

    F = Figure(resolution=(800,500))
    ax = Axis(F[1,1],
    xticks = 2005:3:2020,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    )

    #signal background

    scatter!(ax,years,signal .+ 3.0,linewidth=lw,
    color=col_signal,marker=:x,markersize=ms)
    lines!(ax,years,signal .+ 3.0,linewidth=1,
    color=col_signal,marker=:x,markersize=ms,label="signal")
    #scatter!(ax,years,signal .+ 12.0,linewidth=lw,
    #color=col_signal,marker=:x,markersize=ms,label="signal")

    #fundamental
    lines!(ax,years,fundamental .+ 3.0,linewidth=lw,
    color=col_fundamental,linestyle=:solid,label="fundamental")
    lines!(ax,years,res_fund .+ 0.0,linewidth=lw_s,
    color=col_fundamental,linestyle=:solid)
    #harmonic
    lines!(ax,years,harmonic .+ 3.0,linewidth=lw,
    color=col_harmonic,linestyle=:solid,label="harmonic")
    lines!(ax,years,res_harm .+ 0.0,linewidth=lw_s,
    color=col_harmonic,linestyle=:solid)

    hideydecorations!(ax)
    hidespines!(ax, :t, :r, :l) # only top and right
    #axislegend(ax)
    #axislegend(ax)
    save(savedir*"$(spot)_seasons.png",F)

    #now the additional FFT

    ms = 8
    odi = 4.5 #log offset 
    normalizer(x) = x./maximum(x)
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
    freqstart = findall(x->x>=1/12,freqs)[1]
    freqend = findall(x->x>=6,freqs)[1]
    freq_domain = freqs[freqstart:freqend]

    spec_fun = (abs.(fft(fundamental) |> fftshift)[freqstart:freqend] |> normalizer ) .* 10^odi
    spec_res_fun = (abs.(fft(res_fund) |> fftshift)[freqstart:freqend] |> normalizer )
    spec_har = (abs.(fft(harmonic) |> fftshift)[freqstart:freqend] |> normalizer )  .* 10^odi
    spec_res_har = (abs.(fft(res_harm) |> fftshift)[freqstart:freqend] |> normalizer )
    spec_signal = (abs.(fft(signal) |> fftshift)[freqstart:freqend] |> normalizer )

    F = Figure(resolution=(800,500))
    ax = Axis(F[1,1],
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(7),
    yscale = log10)

    #signal
    #scatter!(ax,freq_domain,spec_signal,linewidth=lw,
    #color=col_signal,marker=:x,markersize=ms)
    scatter!(ax,freq_domain,spec_signal .*10^odi,linewidth=lw,
    color=col_signal,marker=:x,markersize=ms,label="signal")
    #lines!(ax,freq_domain,spec_signal,linewidth=lw_s,
    #color=col_signal,linestyle=:solid)
    lines!(ax,freq_domain,spec_signal .* 10^odi,linewidth=lw_s,
    color=col_signal,linestyle=:solid)

    #fundamental
    lines!(ax,freq_domain,spec_fun,linewidth=lw,
    color=col_fundamental,linestyle=:solid)
    lines!(ax,freq_domain,spec_res_fun,linewidth=lw,
    color=col_fundamental,linestyle=:solid)

    #harmonic
    lines!(ax,freq_domain,spec_har,linewidth=lw,
    color=col_harmonic,linestyle=:solid)
    lines!(ax,freq_domain,spec_res_har,linewidth=lw,
    color=col_harmonic,linestyle=:solid)

    hideydecorations!(ax)
    hidespines!(ax, :t, :r, :l) # only top and right
    #axislegend(ax)

    save(savedir*"$(spot)_seasons_spec.png",F)


    """
    plotting procedure : nlsa versus ssa
    """


    F = Figure(resolution=(800,500))
    ax = Axis(F[1,1],
    xticks = 2005:3:2020,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    )


    #signal background

    scatter!(ax,years,signal .+ offset,linewidth=lw,
    color=col_signal,marker=:x,markersize=ms)
    lines!(ax,years,signal .+ offset,linewidth=lw_s,
    color=col_signal,marker=:x,markersize=ms)
    #scatter!(ax,years,signal .+ 12.0,linewidth=lw,
    #color=col_signal,marker=:x,markersize=ms,label="signal")

    #ssa
    lines!(ax,years,ssa_trend_harm .+ offset,linewidth=lw,
    color=col_ssa,linestyle=:solid,label="ssa")
    lines!(ax,years,signal .- ssa_trend_harm .+ 0,linewidth=lw_s,
    color=col_ssa,linestyle=:solid,label="ssa")

    #nlsa
    lines!(ax,years,nlsa_trend_harm .+ offset,linewidth=lw,
    color=col_nlsa,linestyle=:solid,label="nlsa")
    lines!(ax,years,signal .- nlsa_trend_harm,linewidth=lw_s,
    color=col_nlsa,linestyle=:solid,label="nlsa")

    hideydecorations!(ax)
    hidespines!(ax, :t, :r, :l) # only top and right
    #axislegend(ax)

    save(savedir*"$(spot)_comp.png",F)

    #now the additional FFT

    ms = 8
    odi = 4.5
    normalizer(x) = x./maximum(x)
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
    freqstart = findall(x->x>=1/12,freqs)[1]
    freqend = findall(x->x>=6,freqs)[1]
    freq_domain = freqs[freqstart:freqend]

    spec_ssa = (abs.(fft(ssa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer) .* 10^odi
    spec_res_ssa = (abs.(fft(signal .- ssa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer)
    spec_nlsa = (abs.(fft(nlsa_trend_harm) |> fftshift)[freqstart:freqend]  |> normalizer) .* 10^odi
    spec_res_nlsa = (abs.(fft(signal .- nlsa_trend_harm) |> fftshift)[freqstart:freqend] |> normalizer )
    spec_signal = (abs.(fft(signal) |> fftshift)[freqstart:freqend] |> normalizer )


    F = Figure(resolution=(800,500))
    ax = Axis(F[1,1],
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(7),
    yscale = log10)

    #signal
    #scatter!(ax,freq_domain,spec_signal,linewidth=lw,
    #color=col_signal,marker=:x,markersize=ms)
    scatter!(ax,freq_domain,spec_signal .*10^odi,linewidth=lw,
    color=col_signal,marker=:x,markersize=ms,label="signal")
    #lines!(ax,freq_domain,spec_signal,linewidth=lw_s,
    #color=col_signal,linestyle=:solid)
    lines!(ax,freq_domain,spec_signal .* 10^odi,linewidth=lw_s,
    color=col_signal,linestyle=:solid)

    #ssa
    lines!(ax,freq_domain,spec_ssa,linewidth=lw,
    color=col_ssa,linestyle=:solid)
    lines!(ax,freq_domain,spec_res_ssa,linewidth=lw,
    color=col_ssa,linestyle=:solid)

    #nlsa
    lines!(ax,freq_domain,spec_nlsa,linewidth=lw,
    color=col_nlsa,linestyle=:solid)
    lines!(ax,freq_domain,spec_res_nlsa,linewidth=lw,
    color=col_nlsa,linestyle=:solid)

    hideydecorations!(ax)
    hidespines!(ax, :t, :r, :l) # only top and right
    #axislegend(ax)

    save(savedir*"$(spot)_comp_spec.png",F)


end
