"""
use gaussian estimation of harmonicity peak again
in order to disinguish the results of 1 harmonic and multiples
"""

include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")
using CairoMakie

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

for spot=[1,3,5,9,16],vari=[1,2,4,5], method=["ssa","diff"]


#spot = 1
#vari = 1

#method = "diff"


println("$spot \t $vari")
F = Figure(resolution=(1200,800))

Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
signal = file["signal"]
lambda = file["lambda"]
indices = sortperm(lambda,rev=true)
Eof = file["EOF"][:,indices]
PC = file["PC"][:,indices]
RC = file["RC"][:,indices]
lambda = lambda[indices]


li_harmonics=[1,2]

val = (round(sum(lambda[li_harmonics]),digits=3))
le = length(li_harmonics)
val_m = (round(sum(lambda[li_mixed]),digits=3))
le_m = length(li_mixed)

trend = sum(RC[:,li_harmonics],dims=2)[:]


freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
freqstart = findall(x->x>=1/12,freqs)[1]
freqend = findall(x->x>=6,freqs)[1]
freq_domain = freqs[freqstart:freqend]

Four = fft(signal) |> fftshift
signal_domain = abs.(Four[freqstart:freqend]) .* 1000 .+ 10^-6
Four = fft(trend) |> fftshift
trend_domain = abs.(Four[freqstart:freqend]) .+ 10^-6
Four = fft(signal .- trend) |> fftshift
res_domain = abs.(Four[freqstart:freqend]) ./ 1000 .+ 10^-6

ax_nlsa = Axis(F[1,1:2],title="the first yearly double mode captures $val",
yticksvisible = false,yticklabelsvisible = false,xlabel="t (d)",ylabel="GPP")
ax_nlsa_f = Axis(F[1,3],title="FFT spec",yscale=log10,
yticksvisible = false,yticklabelsvisible = false,xlabel="f (/a)",ylabel="log rel power")

lines!(ax_nlsa,1:N,signal .- trend .-4 ,color="red",label="residual")
lines!(ax_nlsa,1:N,signal,color="grey",label="signal")
lines!(ax_nlsa,1:N,trend,color="blue",label="trend")
axislegend(ax_nlsa)

lines!(ax_nlsa_f,freq_domain,res_domain,color="red")
lines!(ax_nlsa_f,freq_domain,signal_domain,color="grey")
lines!(ax_nlsa_f,freq_domain,trend_domain,color="blue")

scatter!(ax_nlsa_f,freq_domain,res_domain,color="red",marker='x')
scatter!(ax_nlsa_f,freq_domain,signal_domain,color="grey",marker='+')
scatter!(ax_nlsa_f,freq_domain,trend_domain,color="blue",marker='o')


Label(F[0, :], text = "comparison single mode and harmonics | $(meta[spot,"name"]) $(meta[spot,"IGBP_class"]) $(variable_list[vari]) $(method)", fontsize = 30)


li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)



val = (round(sum(lambda[li_harmonics]),digits=3))
le = length(li_harmonics)
val_m = (round(sum(lambda[li_mixed]),digits=3))
le_m = length(li_mixed)

trend = sum(RC[:,li_harmonics],dims=2)[:]


freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
freqstart = findall(x->x>=1/12,freqs)[1]
freqend = findall(x->x>=6,freqs)[1]
freq_domain = freqs[freqstart:freqend]

Four = fft(signal) |> fftshift
signal_domain = abs.(Four[freqstart:freqend]) .* 1000 .+ 10^-6
Four = fft(trend) |> fftshift
trend_domain = abs.(Four[freqstart:freqend]) .+ 10^-6
Four = fft(signal .- trend) |> fftshift
res_domain = abs.(Four[freqstart:freqend]) ./ 1000 .+ 10^-6

ax_ssa = Axis(F[2,1:2],title="$le identified trend modes capture $val | $val_m lost by mixing of $le_m modes",
yticksvisible = false,yticklabelsvisible = false,xlabel="t (d)",ylabel="GPP")
ax_ssa_f = Axis(F[2,3],title="FFT spec",yscale=log10,
yticksvisible = false,yticklabelsvisible = false,xlabel="f (/a)",ylabel="log rel power")

lines!(ax_ssa,1:N,signal .- trend .- 4,color="red")
lines!(ax_ssa,1:N,signal,color="grey")

lines!(ax_ssa,1:N,trend,color="blue")


lines!(ax_ssa_f,freq_domain,res_domain,color="red")
lines!(ax_ssa_f,freq_domain,signal_domain,color="grey")
lines!(ax_ssa_f,freq_domain,trend_domain,color="blue")

scatter!(ax_ssa_f,freq_domain,res_domain,color="red",marker='x')
scatter!(ax_ssa_f,freq_domain,signal_domain,color="grey",marker='+')
scatter!(ax_ssa_f,freq_domain,trend_domain,color="blue",marker='o')

save(savedir*"$(spot)_$(vari)_$(method).png",F)
#save(dir*"test.png",F)

end