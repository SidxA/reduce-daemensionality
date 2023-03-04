include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")
using CairoMakie

"""
plot the gaussian estimation of the harmonicity of the FFt spectra

the peak picture

"""


# link 5,8,9 !


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

method = "ssa"
spot = 1
savedirname = dir*"poster/"


Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
signal = file["signal"]
lambda = file["lambda"]
indices = sortperm(lambda,rev=true)
Eof = file["EOF"][:,indices]
PC = file["PC"][:,indices]
RC = file["RC"][:,indices]
lambda = lambda[indices]


peak_color = "blue"
freq_color = "orange"
limits = (0.9,5.1,0,1.1)

F = Figure(resolution=(800,500))


    ax1   = Axis(F[1,1],limits=limits)
    ax2   = Axis(F[2,1],limits=limits)
    ax3   = Axis(F[3,1],limits=limits,
    xminorticksvisible = true,
    xminorticks = 1:1/7:7)

    #first: clean peak
    mode = Eof[:,5]
    Four = fft(mode) |> fftshift
    freqs = fftfreq(length(tw), 1.0/Ts) |> fftshift
    freqstart = findall(x->x>=1/12,freqs)[1]
    freqend = findall(x->x>=6,freqs)[1]
    val_domain = abs.(Four[freqstart:freqend])
    val_domain ./= maximum(val_domain)
    freq_domain = freqs[freqstart:freqend]
    freq, value,sigma = fit_gauss(freq_domain,val_domain)
    peak = gauss(freq_domain,(freq,value,sigma)) .+ 10^-3
    residual = val_domain .- peak
    freq_round = round(freq,digits=1)
    val = abs(1/sigma /sqrt(2*pi))

    scatter!(ax1,freq_domain,val_domain,color=freq_color,marker=:x,markersize=12,label="mode")
    lines!(ax1,freq_domain,val_domain,color=freq_color,linediwdth=2)
    lines!(ax1,freq_domain,peak,color=peak_color,marker=:+,
    linestyle=:dot,linewidth=5,label="gauss")

    #second: peak + residual
    mode = Eof[:,9]
    Four = fft(mode) |> fftshift
    freqs = fftfreq(length(tw), 1.0/Ts) |> fftshift
    freqstart = findall(x->x>=1/12,freqs)[1]
    freqend = findall(x->x>=6,freqs)[1]
    val_domain = abs.(Four[freqstart:freqend])
    val_domain ./= maximum(val_domain)
    freq_domain = freqs[freqstart:freqend]
    freq, value,sigma = fit_gauss(freq_domain,val_domain)
    peak = gauss(freq_domain,(freq,value,sigma)) .+ 10^-3
    residual = val_domain .- peak
    freq_round = round(freq,digits=1)
    val = abs(1/sigma /sqrt(2*pi))

    scatter!(ax2,freq_domain,val_domain,color=freq_color,marker=:x,markersize=12,label="mode")
    lines!(ax2,freq_domain,val_domain,color=freq_color,linediwdth=2)
    lines!(ax2,freq_domain,peak,color=peak_color,marker=:+,
    linestyle=:dot,linewidth=5,label="gauss")

    #third: peak not at harmonic
    mode = Eof[:,8]
    Four = fft(mode) |> fftshift
    freqs = fftfreq(length(tw), 1.0/Ts) |> fftshift
    freqstart = findall(x->x>=1/12,freqs)[1]
    freqend = findall(x->x>=6,freqs)[1]
    val_domain = abs.(Four[freqstart:freqend])
    val_domain ./= maximum(val_domain)
    freq_domain = freqs[freqstart:freqend]
    freq, value,sigma = fit_gauss(freq_domain,val_domain)
    peak = gauss(freq_domain,(freq,value,sigma)) .+ 10^-3
    residual = val_domain .- peak
    freq_round = round(freq,digits=1)
    val = abs(1/sigma /sqrt(2*pi))

    scatter!(ax3,freq_domain,val_domain,color=freq_color,marker=:x,markersize=12,label="mode")
    lines!(ax3,freq_domain,val_domain,color=freq_color,linediwdth=2)
    lines!(ax3,freq_domain,peak,color=peak_color,marker=:+,
    linestyle=:dot,linewidth=5,label="gauss")

linkyaxes!(ax1,ax2)
linkyaxes!(ax2,ax3)



hidedecorations!(ax1)
hidedecorations!(ax2)
#hidedecorations!(ax3,label=false,ticks=false,ticklabels=false)
hideydecorations!(ax3)

vlines!(ax1, 1:8, color = :grey)
vlines!(ax2, 1:8, color = :grey)
vlines!(ax3, 1:8, color = :grey)

save(dir*"peaks.png",F)
