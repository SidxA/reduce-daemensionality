include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")
using CairoMakie


"""
algortihm for fitting the fft plots using cauchy peak
"""

function cauchy(x,p) # cauchy peak with two parameters
    x0 = p[1]
    gamma = p[2]
    @. return 1/(pi*gamma*(1+((x-x0)/gamma)^2))
end

function fit_cauchy(yvalues,xvalues) # fit the cauchy peak to a vec of numbers
    onebins = xvalues
    bins = yvalues
    p0 = ones(2)
    return coef(curve_fit(cauchy,bins,onebins,p0))
end



function mode_harmonicity_estimation_cauchy(eof,threshold)
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

        #Four ./= maximum(abs.(Four)) # apparently works better when un-normalized

        #freqend = findall(x->x>12,freqs)[1] # no faster then a freq of 12 per year

        val_domain = abs.(Four[freqstart:freqend])
        freq_domain = freqs[freqstart:freqend]

        freq, value = fit_cauchy(freq_domain,val_domain)
        peak = cauchy(freq_domain,(freq,value))
        residual = val_domain .- peak

        freq_round = round(freq,digits=2)
        val_conv = 1/value/pi
        println("$i \t $(freq_round) \t $(val_conv)")

        #if maximal residual value other then single peak is larger then relative peak threshold value, its 
        if maximum(residual .+ 0.0)*threshold <= 1/ value / pi &&  freq_round in Float64.(1:10)
            li_harmonics = append!(li_harmonics,i)
            li_h_freq = append!(li_h_freq,freq)
        elseif maximum(residual .+ 0.0)*threshold >= 1/ value / pi && freq_round in Float64.(1:10)
            li_mixed = append!(li_mixed,i)
            li_m_freq = append!(li_m_freq,freq)
        elseif maximum(residual .+ 0.0)*threshold >= 1/ value / pi
            li_residual = append!(li_residual,i)
        else
            println("no peak")
        end

    end

    return li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual
end

"""
algortihm for fitting the fft plots using gauss peak
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



function mode_harmonicity_estimation_gauss(eof,threshold)
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
            println("no peak")
        end

    end

    return li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual
end


"""
print the results for different threshold values
"""
spot = 19
outdir="/net/scratch/lschulz/fluxfullset_1a/"
#meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
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

#freqstart = Int(round(W/2,digits=0))-1
#freqend = Int(round(W/2 + W/4,digits=0))-1

tw = t0:Ts:(t0+(W-1)*Ts)
#meta = load("/net/home/lschulz/logs/KW_2_06/meta.jld2")["meta"]

outdir="/net/scratch/lschulz/fluxfullset/"
method = "ssa"

Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
lambda = file["lambda"]
indices = sortperm(lambda,rev=true)

Eof = file["EOF"][:,indices]
lambda = lambda[indices]

for threshold in 1:100
    li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation(Eof,threshold)
    println("threshold $threshold \t trend $(sum(lambda[li_harmonics])) \t mix $(sum(lambda[li_mixed]))")
end

"""
plot the results
"""
spot= 1
threshold = 10
method = "ssa"
Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
lambda = file["lambda"]
indices = sortperm(lambda,rev=true)
Eof = file["EOF"][:,indices]
lambda = lambda[indices]


li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)


method = "diff"

Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
lambda = file["lambda"]
indices = sortperm(lambda,rev=true)
Eof = file["EOF"][:,indices]
lambda = lambda[indices]

li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)


"""
it seems to be the case of a swith at around 1/10 of the value
there is substantial amplitude modulation on it, already on the easy ones

lets take threshold = 20 for starters

plot the trending selection results
"""

savedir=dir*"trends/"

for spot=1:20

    name = "$(spot)_7a"


    freqlimits = (1/12,6,-3,1.1)

    F = Figure(resolution=(800,800))

    xtic=1:6

    title1="$spot  $(meta[spot,"name"]) $(meta[spot,"IGBP_class"]) W = $(round(W/365,digits=0)) a"
    subtitle1 = "time series"

    #axis


    ax_modes_ssa = Axis(F[1,1:2],xlabel="time",ylabel="GPP",yticksvisible = false,
    yticklabelsvisible = false,title = "SSA modes",
    titlealign = :left,
    subtitlegap = 2,
    titlegap = 5,
    subtitlefont = :italic,
    subtitlelineheight = 0.9,
    titlelineheight = 0.9,)


    ax_freq_ssa = Axis(F[1,3:4],limits=freqlimits,#,yscale=log10
    xticks=xtic,xlabel="frequency",ylabel ="rel power",yticksvisible = false,
    yticklabelsvisible = false,title = "SSA spectrum",
    titlealign = :left,
    subtitlegap = 2,
    titlegap = 5,
    subtitlefont = :italic,
    subtitlelineheight = 0.9,
    titlelineheight = 0.9,)

    ax_modes_nlsa = Axis(F[2,1:2],xlabel="time",ylabel="GPP",yticksvisible = false,
    yticklabelsvisible = false,title = "NLSA modes",
    titlealign = :left,
    subtitlegap = 2,
    titlegap = 5,
    subtitlefont = :italic,
    subtitlelineheight = 0.9,
    titlelineheight = 0.9,)


    ax_freq_nlsa = Axis(F[2,3:4],limits=freqlimits,#,yscale=log10
    xticks=xtic,xlabel="frequency",ylabel ="rel power",yticksvisible = false,
    yticklabelsvisible = false,title = "NLSA spectrum",
    titlealign = :left,
    subtitlegap = 2,
    titlegap = 5,
    subtitlefont = :italic,
    subtitlelineheight = 0.9,
    titlelineheight = 0.9,)

    #signal


    signal = load(create_file_list(outdir,"ssa",W,vari,preproc)[spot])["signal"]

    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

    freqstart = findall(x->x>=1/12,freqs)[1]
    freqend = findall(x->x>=6,freqs)[1]

    freqs = freqs[freqstart:freqend]

    Four = fft(signal) |> fftshift
    Four = abs.(Four[freqstart:freqend])

    Four ./= maximum(abs.(Four))


    lines!(ax_modes_ssa,1:N,signal,color="grey")
    lines!(ax_freq_ssa,freqs,Four,color="grey")

    lines!(ax_modes_nlsa,1:N,signal,color="grey")
    lines!(ax_freq_nlsa,freqs,Four,color="grey")


    #modes ssa

    method = "ssa"
    Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
    file = load(Filename)
    lambda = file["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file["EOF"][:,indices]
    lambda = lambda[indices]

    li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)

    trend = sum(RC[:,li_harmonics],dims=2)[:]   .+ 0.0
    mix = sum(RC[:,li_mixed],dims=2)[:] .+ 0.0
    residual = signal .- trend .- mix

    series!(ax_modes_ssa,1:N,[trend.-3 mix.-6 residual.-9]')

    trendfft = fft(trend) |> fftshift
    trendfft = abs.(trendfft)[freqstart:freqend]
    trendfft ./= maximum(abs.(trendfft)) .+ 0.0

    mixfft = abs.(fft(mix) |> fftshift)
    mixfft = abs.(mixfft)[freqstart:freqend]
    mixfft ./= maximum(abs.(mixfft)) .+ 0.0

    residualfft = abs.(fft(residual) |> fftshift)
    residualfft = abs.(residualfft)[freqstart:freqend]
    residualfft ./= maximum(abs.(residualfft)) .+ 0.0

    series!(ax_freq_ssa,freqs,[trendfft.-1 mixfft.-2 residualfft.-3]')



    # modes nlsa

    method = "diff"

    Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
    file = load(Filename)
    lambda = file["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file["EOF"][:,indices]
    lambda = lambda[indices]

    li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)

    trend = sum(RC[:,li_harmonics],dims=2)[:] .+ 0.0
    mix = sum(RC[:,li_mixed],dims=2)[:] .+ 0.0
    residual = signal .- trend .- mix

    series!(ax_modes_nlsa,1:N,[trend.-3 mix.-6 residual.-9]')


    trendfft = fft(trend) |> fftshift
    trendfft = abs.(trendfft)[freqstart:freqend]
    trendfft ./= maximum(abs.(trendfft)) .+ 0.0

    mixfft = abs.(fft(mix) |> fftshift)
    mixfft = abs.(mixfft)[freqstart:freqend]
    mixfft ./= maximum(abs.(mixfft)) .+ 0.0

    residualfft = abs.(fft(residual) |> fftshift)
    residualfft = abs.(residualfft)[freqstart:freqend]
    residualfft ./= maximum(abs.(residualfft)) .+ 0.0

    series!(ax_freq_nlsa,freqs,[trendfft.-1 mixfft.-2 residualfft.-3]')


    #layout

    colgap!(F.layout,1,0)
    rowgap!(F.layout,1,0)
    #rowgap!(F.layout,2,0)

    save(savedir*name*".png",F) #save(savedir*name*".png",F)

end
