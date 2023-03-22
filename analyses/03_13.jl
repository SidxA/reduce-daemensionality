
include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")
using CairoMakie

"""
all the plots
"""

#uses these functions
function single_robustness(outdir,method,W,vari,preproc,kappa,spot,yearsamples)
    N = 5478
    k = 48
    filename = create_file_list(outdir,method,W,vari,preproc)[spot]
    lambda,protofreq,RC = extract_from_single_file(filename,yearsamples,N,k)
    lambda_l = sort(lambda,rev=true)[1:kappa]
    protofreq_l = protofreq[sortperm(lambda,rev=true)][1:kappa]
    RC_l = RC[:,sortperm(lambda,rev=true)][:,1:kappa]
    return RC_l,lambda_l,protofreq_l
end

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

function plotting_routine(spot,savedirname)


    #prerequisites
    #spot = 1
    outdir="/net/scratch/lschulz/fluxfullset/"
    #meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
    W = 2556
    preproc = "raw."
    kappa= 48
    yearsamples=365.25
    vari = 1
    L = 5478
    bins = 10 .^Array(-2:0.01:2)
    # Number of points
    N= 5478 
    # Sample period
    Ts = 1 / 365.25
    # Start time 
    t0 = 0 
    tmax = t0 + (N-1) * Ts
    # time coordinate
    t = t0:Ts:tmax
    tw = t0:Ts:(t0+(W-1)*Ts)

    color_ssa = "darkgreen"
    color_nlsa = "purple"
    color_signal = "grey50"
    lw = 3
    ms = 5

    years = ((1:5478) ./ 365) .+ 2005

    # plotting parameters

    meta = load("/net/home/lschulz/logs/KW_2_06/meta.jld2")["meta"]

    igbp = meta[spot,"IGBP_class"]
    name = meta[spot,"name"]
    variable = ["GPP"][vari]

    titlestring = "$variable in $spot $name $igbp"
    savedirstring = "$(variable)_$(spot)_"
    savedir = savedirname*savedirstring

    #load files
    Filename_ssa = create_file_list(outdir,"ssa",W,vari,preproc)[spot]
    Filename_nlsa = create_file_list(outdir,"diff",W,vari,preproc)[spot]
    file_ssa = load(Filename_ssa)
    rc_ssa = sum(file_ssa["RC"],dims=2)[:]
    file_nlsa = load(Filename_nlsa)
    rc_nlsa = sum(file_nlsa["RC"],dims=2)[:]

    signal = file_ssa["signal"]

    #FFT stuff
    ms = 8
    odi = 4.5 #log offset 
    normalizer(x) = x./maximum(x)
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
    freqstart = findall(x->x>=1/12,freqs)[1]
    freqend = findall(x->x>=6,freqs)[1]
    freq_domain = freqs[freqstart:freqend]


    spec_signal = (abs.(fft(signal) |> fftshift)[freqstart:freqend] |> normalizer )
    spec_ssa = (abs.(fft(rc_ssa) |> fftshift)[freqstart:freqend] |> normalizer )
    spec_nlsa = (abs.(fft(rc_nlsa) |> fftshift)[freqstart:freqend] |> normalizer )


    """
    first one signal rec + fft
    """

    #figure
    F = Figure(resolution=(800,800))


    ax_time = Axis(F[1,1],xticks=2005:3:2020,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel="time (a)",
    ylabel="GPP")

    ax_spec = Axis(F[2,1],yscale=log10,
    xlabel="frequency (/a)",
    ylabel="relative power",)

    scatter!(ax_time,years,signal,color=color_signal,markersize = ms,marker=:x)
    lines!(ax_time,years,signal,color=color_signal,markersize = 1,label = "signal")
    lines!(ax_time,years,rc_ssa,color=color_ssa,linewidth=lw,label="SSA")
    lines!(ax_time,years,rc_nlsa,color=color_nlsa,linewidth=lw,label="NLSA")

    lines!(ax_spec,freq_domain,spec_signal,color=color_signal,markersize = 1,label = "signal")
    lines!(ax_spec,freq_domain,spec_ssa,color=color_ssa,linewidth=lw,label="SSA")
    lines!(ax_spec,freq_domain,spec_nlsa,color=color_nlsa,linewidth=lw,label="NLSA")

    axislegend(ax_spec)

    #hideydecorations!(ax)
    #hidespines!(ax,:t,:r)

    #axislegend(ax)

    Label(F[0, 1], text = "reconstructions for "*titlestring,
    fontsize = 24)

    #trim!(F.layout)

    save(savedir*"reconstructions.png",F)

    """
    EOF + FFT
    """
    #new fft stuff for the W not the N !

    xtic=1:6
    years_modes = (1:W) ./ 365
    freqstart = Int(round(W/2,digits=0))+1
    tw = t0:Ts:(t0+(W-1)*Ts)
    spot=1

    #figure
    F = Figure(resolution=(800,800))


    ax_modes = Axis(F[1,1],yticksvisible = false,
    yticklabelsvisible = false,xticks=1:7,
    xlabel="time (a)",
    ylabel="individual modes")


    ax_freq = Axis(F[1,2],limits=(1/10, 7,-0.7,37),#,yscale=log10
    xticks=xtic,yticksvisible = false,
    yticklabelsvisible = false,
    xlabel="frequency (/a)",
    ylabel="relative power")

    #modes ssa

    lambda = file_ssa["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file_ssa["EOF"][:,indices]

    for k = 1:16

        mode = Eof[:,k] 
        Four = fft(mode) |> fftshift
        freqs = fftfreq(length(tw), 1.0/Ts) |> fftshift
        Four ./= maximum(abs.(Four)) 

        lines!(ax_modes,years_modes,mode .+ (k*0.07),
        color=color_ssa)
        lines!(ax_freq,freqs[freqstart:end],abs.(Four)[freqstart:end] .+ k,
        color=color_ssa)
    end

    #modes nlsa

    lambda = file_nlsa["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file_nlsa["EOF"][:,indices]

    for k = 1:16

        mode = Eof[:,k] 
        Four = fft(mode) |> fftshift
        freqs = fftfreq(length(tw), 1.0/Ts) |> fftshift
        Four ./= maximum(abs.(Four)) 

        lines!(ax_modes,years_modes,mode .+ ((18+k)*0.07),
        color=color_nlsa)
        lines!(ax_freq,freqs[freqstart:end],abs.(Four)[freqstart:end] .+ k .+ 18,
        color=color_nlsa)
    end

    #colgap!(F.layout,1,0)

    #hidespines!(ax_freq, :t, :r, :l)
    #hidespines!(ax_modes, :t, :r, :l)

    Label(F[0, 1:2], text = "EOF with FFT for "*titlestring,
    fontsize = 24)

    #trim!(F.layout)

    save(savedir*"modeshape.png",F)

    """
    seasons + FFT
    """

    #sort by threshold

    threshold = 10
    vari=1

    # ssa reconstruction: freq harmonics

    file = load(Filename_ssa)
    signal = file["signal"]
    lambda = file["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file["EOF"][:,indices]
    PC = file["PC"][:,indices]
    RC = file["RC"][:,indices]
    lambda = lambda[indices]

    li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)
    val = (round(sum(lambda[li_harmonics]),digits=3))
    le = length(li_harmonics)
    #val_m = (round(sum(lambda[li_mixed]),digits=3))
    #le_m = length(li_mixed)
    ssa_trend_harm = sum(RC[:,li_harmonics],dims=2)[:]

    # nlsa reconstruction: freq harmonics

    file = load(Filename)
    signal = file["signal"]
    lambda = file["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file["EOF"][:,indices]
    PC = file["PC"][:,indices]
    RC = file["RC"][:,indices]
    lambda = lambda[indices]

    li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)
    val = (round(sum(lambda[li_harmonics]),digits=3))
    le = length(li_harmonics)
    #val_m = (round(sum(lambda[li_mixed]),digits=3))
    #le_m = length(li_mixed)
    nlsa_trend_harm = sum(RC[:,li_harmonics],dims=2)[:]

    #the fundamental

    li_harmonics=[1,2]
    val = (round(sum(lambda[li_harmonics]),digits=3))
    le = length(li_harmonics)
    #val_m = (round(sum(lambda[li_mixed]),digits=3))
    #le_m = length(li_mixed)
    trend_first = sum(RC[:,li_harmonics],dims=2)[:]

    #the figure

    offset = 2.5
    lw = 3
    lw_s = 1
    ms = 4

    F = Figure(resolution=(800,800))


    ax_time = Axis(F[1,1],
    xticks = 2005:3:2020,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel = "time (a)",
    ylabel = "GPP"
    )


    #signal background

    scatter!(ax_time,years,signal .+ offset,linewidth=lw,
    color=color_signal,marker=:x,markersize=ms)
    lines!(ax_time,years,signal .+ offset,linewidth=lw_s,
    color=color_signal,marker=:x,markersize=ms)
    #scatter!(ax_time,years,signal .+ 12.0,linewidth=lw,
    #color=color_signal,marker=:x,markersize=ms,label="signal")

    #ssa
    lines!(ax_time,years,ssa_trend_harm .+ offset,linewidth=lw,
    color=color_ssa,linestyle=:solid,label="ssa")
    lines!(ax_time,years,signal .- ssa_trend_harm .+ 0,linewidth=lw_s,
    color=color_ssa,linestyle=:solid,label="ssa")

    #nlsa
    lines!(ax_time,years,nlsa_trend_harm .+ offset,linewidth=lw,
    color=color_nlsa,linestyle=:solid,label="nlsa")
    lines!(ax_time,years,signal .- nlsa_trend_harm,linewidth=lw_s,
    color=color_nlsa,linestyle=:solid,label="nlsa")

    #hideydecorations!(ax_time)
    #hidespines!(ax_time, :t, :r, :l) # only top and right
    #axislegend(ax_time)

    #save(savedir*"$(spot)_comp.png",F)

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


    ax_spec = Axis(F[2,1],
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(7),
    yscale = log10,
    xlabel = "frequency (/a)",
    ylabel = "relative power")

    #signal
    #scatter!(ax_spec,freq_domain,spec_signal,linewidth=lw,
    #color=color_signal,marker=:x,markersize=ms)
    scatter!(ax_spec,freq_domain,spec_signal .*10^odi,linewidth=lw,
    color=color_signal,marker=:x,markersize=ms,label="signal")
    #lines!(ax_spec,freq_domain,spec_signal,linewidth=lw_s,
    #color=color_signal,linestyle=:solid)
    lines!(ax_spec,freq_domain,spec_signal .* 10^odi,linewidth=lw_s,
    color=color_signal,linestyle=:solid)

    #ssa
    lines!(ax_spec,freq_domain,spec_ssa,linewidth=lw,
    color=color_ssa,linestyle=:solid)
    lines!(ax_spec,freq_domain,spec_res_ssa,linewidth=lw,
    color=color_ssa,linestyle=:solid)

    #nlsa
    lines!(ax_spec,freq_domain,spec_nlsa,linewidth=lw,
    color=color_nlsa,linestyle=:solid)
    lines!(ax_spec,freq_domain,spec_res_nlsa,linewidth=lw,
    color=color_nlsa,linestyle=:solid)

    #hideydecorations!(ax_spec)
    #hidespines!(a_spec, :t, :r, :l) # only top and right
    #axislegend(ax_spec)

    Label(F[0, 1], text = "seasons + FFT for "*titlestring,
    fontsize = 24)

    #trim!(F.layout)

    save(savedir*"seasonal.png",F)

    """
    the run on the years
    """


    year_ind =
    [
    1:365,
    366:730,
    731:1095,
    1096:1460,
    1462:1826, #89 1 offset
    1827:2191, 
    2192:2556,
    2557:2921,
    2923:3287, # 93 1 offset
    3288:3652,
    3653:4017,
    4018:4382,
    4384:4748, #97 1 offset
    4749:5113,
    5114:5478,
    ]

    years_plot = 2005 .+ (1:15)

    color_signal_grad =  cgrad([:grey20,:grey80],15,categorical=true,rev=true)
    color_ssa_grad =  cgrad([:darkgreen,:lightgreen],15,categorical=true,rev=true)
    color_nlsa_grad =  cgrad([:purple,:plum],15,categorical=true,rev=true)

    lw = 1


    fund = hcat([trend_first[i] for i in year_ind]...)

    ssa_harm = hcat([ssa_trend_harm[i] for i in year_ind]...)
    nlsa_harm = hcat([nlsa_trend_harm[i] for i in year_ind]...)


    f = Figure(resolution=(800,800),)

    ax_years = Axis(f[1, 1],
        yticks = ((1:15) ./ 2 .-1,  string.(years_plot)),
        xticks = ([365/2, 400 + 365/2, 800 + 365/2],
        ["fund","SSA harm","NLSA harm"]),
        ylabel="years"
    )


    for i in 1:15
        d = lines!(ax_years,1:365,fund[:,i].+(0.5*i),color = color_signal_grad[i],linewidth=lw)#,offset = i / 4,)
        d = lines!(ax_years,400 .+ (1:365),ssa_harm[:,i].+(0.5*i),color = color_ssa_grad[i],linewidth=lw)
        #d = lines!(800 .+ (1:365),nlsa_fund[:,i].+(0.5*i),color = colors[i])
        d = lines!(ax_years,800 .+ (1:365),nlsa_harm[:,i].+(0.5*i),color = color_nlsa_grad[i],linewidth=lw)
    end

    ax_years_2 = Axis(f[2,1],
    xticks = ([365/2, 400 + 365/2, 800 + 365/2],
        ["fund","SSA harm","NLSA harm"]),
        ylabel="GPP"
        )

    for i in 1:15
        d = lines!(ax_years_2,1:365,fund[:,i],color = color_signal_grad[i],linewidth=lw)#,offset = i / 4,)
        d = lines!(ax_years_2,400 .+ (1:365),ssa_harm[:,i],color = color_ssa_grad[i],linewidth=lw)
        #d = lines!(800 .+ (1:365),nlsa_fund[:,i].+(0.5*i),color = colors[i])
        d = lines!(ax_years_2,800 .+ (1:365),nlsa_harm[:,i],color = color_nlsa_grad[i],linewidth=lw)
    end

    Label(f[0, :], text = "changes in seasonal behavior for "*titlestring,
    fontsize = 24)

    #trim!(f.layout)

    save(savedir*"years.png",f)

end

function plotting_routine_1a(spot,savedirname,threshold = 10)


    #prerequisites
    #spot = 1
    outdir="/net/scratch/lschulz/fluxfullset_1a/"
    #meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
    W = 365
    preproc = "raw."
    kappa= 48
    yearsamples=365.25
    vari = 1
    L = 5478
    bins = 10 .^Array(-2:0.01:2)
    # Number of points
    N= 5478 
    # Sample period
    Ts = 1 / 365.25
    # Start time 
    t0 = 0 
    tmax = t0 + (N-1) * Ts
    # time coordinate
    t = t0:Ts:tmax
    tw = t0:Ts:(t0+(W-1)*Ts)

    color_ssa = "darkgreen"
    color_nlsa = "purple"
    color_signal = "grey50"
    lw = 1
    ms = 5

    years = ((1:5478) ./ 365) .+ 2005

    # plotting parameters

    meta = load("/net/home/lschulz/logs/KW_2_06/meta.jld2")["meta"]

    igbp = meta[spot,"IGBP_class"]
    name = meta[spot,"name"]
    variable = ["GPP"][vari]

    titlestring = "$variable in $spot $name $igbp"
    savedirstring = "$(variable)_$(spot)_$(W)d_"
    savedir = savedirname*savedirstring

    #load files
    Filename_ssa = create_file_list(outdir,"ssa",W,vari,preproc)[spot]
    Filename_nlsa = create_file_list(outdir,"diff",W,vari,preproc)[spot]
    file_ssa = load(Filename_ssa)
    rc_ssa = sum(file_ssa["RC"],dims=2)[:]
    file_nlsa = load(Filename_nlsa)
    rc_nlsa = sum(file_nlsa["RC"],dims=2)[:]

    signal = file_ssa["signal"]

    #FFT stuff
    ms = 8
    odi = 4.5 #log offset 
    normalizer(x) = x./maximum(x)
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
    freqstart = findall(x->x>=1/12,freqs)[1]
    freqend = findall(x->x>=6,freqs)[1]
    freq_domain = freqs[freqstart:freqend]


    spec_signal = (abs.(fft(signal) |> fftshift)[freqstart:freqend] |> normalizer )
    spec_ssa = (abs.(fft(rc_ssa) |> fftshift)[freqstart:freqend] |> normalizer )
    spec_nlsa = (abs.(fft(rc_nlsa) |> fftshift)[freqstart:freqend] |> normalizer )


    """
    first one signal rec + fft
    """

    #figure
    F = Figure(resolution=(800,800))


    ax_time = Axis(F[1,1],xticks=2005:3:2020,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel="time (a)",
    ylabel="GPP")

    ax_spec = Axis(F[2,1],yscale=log10,
    xlabel="frequency (/a)",
    ylabel="relative power",)

    scatter!(ax_time,years,signal,color=color_signal,markersize = ms,marker=:x)
    lines!(ax_time,years,signal,color=color_signal,markersize = 1,label = "signal")
    lines!(ax_time,years,rc_ssa,color=color_ssa,linewidth=lw,label="SSA")
    lines!(ax_time,years,rc_nlsa,color=color_nlsa,linewidth=lw,label="NLSA")

    lines!(ax_spec,freq_domain,spec_signal,color=color_signal,markersize = 1,label = "signal")
    lines!(ax_spec,freq_domain,spec_ssa,color=color_ssa,linewidth=lw,label="SSA")
    lines!(ax_spec,freq_domain,spec_nlsa,color=color_nlsa,linewidth=lw,label="NLSA")

    axislegend(ax_spec)

    #hideydecorations!(ax)
    #hidespines!(ax,:t,:r)

    #axislegend(ax)

    Label(F[0, 1], text = "reconstructions for "*titlestring,
    fontsize = 24)

    #trim!(F.layout)

    save(savedir*"reconstructions.png",F)

    """
    EOF + FFT
    """
    #new fft stuff for the W not the N !

    xtic=1:6
    years_modes = (1:W) ./ 365
    freqstart = Int(round(W/2,digits=0))+1
    tw = t0:Ts:(t0+(W-1)*Ts)
    spot=1

    #figure
    F = Figure(resolution=(800,800))


    ax_modes = Axis(F[1,1],yticksvisible = false,
    yticklabelsvisible = false,xticks=1:7,
    xlabel="time (a)",
    ylabel="individual modes")


    ax_freq = Axis(F[1,2],limits=(1/10, 7,-0.7,37),#,yscale=log10
    xticks=xtic,yticksvisible = false,
    yticklabelsvisible = false,
    xlabel="frequency (/a)",
    ylabel="relative power")

    #modes ssa

    lambda = file_ssa["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file_ssa["EOF"][:,indices]

    for k = 1:16

        mode = Eof[:,k] 
        Four = fft(mode) |> fftshift
        freqs = fftfreq(length(tw), 1.0/Ts) |> fftshift
        Four ./= maximum(abs.(Four)) 

        lines!(ax_modes,years_modes,mode .+ (k*0.07),
        color=color_ssa)
        lines!(ax_freq,freqs[freqstart:end],abs.(Four)[freqstart:end] .+ k,
        color=color_ssa)
    end

    #modes nlsa

    lambda = file_nlsa["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file_nlsa["EOF"][:,indices]

    for k = 1:16

        mode = Eof[:,k] 
        Four = fft(mode) |> fftshift
        freqs = fftfreq(length(tw), 1.0/Ts) |> fftshift
        Four ./= maximum(abs.(Four)) 

        lines!(ax_modes,years_modes,mode .+ ((18+k)*0.07),
        color=color_nlsa)
        lines!(ax_freq,freqs[freqstart:end],abs.(Four)[freqstart:end] .+ k .+ 18,
        color=color_nlsa)
    end

    #colgap!(F.layout,1,0)

    #hidespines!(ax_freq, :t, :r, :l)
    #hidespines!(ax_modes, :t, :r, :l)

    Label(F[0, 1:2], text = "EOF with FFT for "*titlestring,
    fontsize = 24)

    #trim!(F.layout)

    save(savedir*"modeshape.png",F)

    """
    seasons + FFT
    """

    #sort by threshold



    # ssa reconstruction: freq harmonics

    file = load(Filename_ssa)
    signal = file["signal"]
    lambda = file["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file["EOF"][:,indices]
    PC = file["PC"][:,indices]
    RC = file["RC"][:,indices]
    lambda = lambda[indices]

    li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)
    val = (round(sum(lambda[li_harmonics]),digits=3))
    le = length(li_harmonics)
    #val_m = (round(sum(lambda[li_mixed]),digits=3))
    #le_m = length(li_mixed)
    ssa_trend_harm = sum(RC[:,li_harmonics],dims=2)[:]

    # nlsa reconstruction: freq harmonics

    file = load(Filename)
    signal = file["signal"]
    lambda = file["lambda"]
    indices = sortperm(lambda,rev=true)
    Eof = file["EOF"][:,indices]
    PC = file["PC"][:,indices]
    RC = file["RC"][:,indices]
    lambda = lambda[indices]

    li_harmonics,li_mixed,li_h_freq,li_m_freq,li_residual = mode_harmonicity_estimation_gauss(Eof,threshold)
    val = (round(sum(lambda[li_harmonics]),digits=3))
    le = length(li_harmonics)
    #val_m = (round(sum(lambda[li_mixed]),digits=3))
    #le_m = length(li_mixed)
    nlsa_trend_harm = sum(RC[:,li_harmonics],dims=2)[:]

    #the fundamental

    li_harmonics=[1,2]
    val = (round(sum(lambda[li_harmonics]),digits=3))
    le = length(li_harmonics)
    #val_m = (round(sum(lambda[li_mixed]),digits=3))
    #le_m = length(li_mixed)
    trend_first = sum(RC[:,li_harmonics],dims=2)[:]

    #the figure

    offset = 2.5
    lw = 3
    lw_s = 1
    ms = 4

    F = Figure(resolution=(800,800))


    ax_time = Axis(F[1,1],
    xticks = 2005:3:2020,
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(3),
    xlabel = "time (a)",
    ylabel = "GPP"
    )


    #signal background

    scatter!(ax_time,years,signal .+ offset,linewidth=lw,
    color=color_signal,marker=:x,markersize=ms)
    lines!(ax_time,years,signal .+ offset,linewidth=lw_s,
    color=color_signal,marker=:x,markersize=ms)
    #scatter!(ax_time,years,signal .+ 12.0,linewidth=lw,
    #color=color_signal,marker=:x,markersize=ms,label="signal")

    #ssa
    lines!(ax_time,years,ssa_trend_harm .+ offset,linewidth=lw,
    color=color_ssa,linestyle=:solid,label="ssa")
    lines!(ax_time,years,signal .- ssa_trend_harm .+ 0,linewidth=lw_s,
    color=color_ssa,linestyle=:solid,label="ssa")

    #nlsa
    lines!(ax_time,years,nlsa_trend_harm .+ offset,linewidth=lw,
    color=color_nlsa,linestyle=:solid,label="nlsa")
    lines!(ax_time,years,signal .- nlsa_trend_harm,linewidth=lw_s,
    color=color_nlsa,linestyle=:solid,label="nlsa")

    #hideydecorations!(ax_time)
    #hidespines!(ax_time, :t, :r, :l) # only top and right
    #axislegend(ax_time)

    #save(savedir*"$(spot)_comp.png",F)

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


    ax_spec = Axis(F[2,1],
    xminorticksvisible = true,
    xminorgridvisible = true,
    xminorticks = IntervalsBetween(7),
    yscale = log10,
    xlabel = "frequency (/a)",
    ylabel = "relative power")

    #signal
    #scatter!(ax_spec,freq_domain,spec_signal,linewidth=lw,
    #color=color_signal,marker=:x,markersize=ms)
    scatter!(ax_spec,freq_domain,spec_signal .*10^odi,linewidth=lw,
    color=color_signal,marker=:x,markersize=ms,label="signal")
    #lines!(ax_spec,freq_domain,spec_signal,linewidth=lw_s,
    #color=color_signal,linestyle=:solid)
    lines!(ax_spec,freq_domain,spec_signal .* 10^odi,linewidth=lw_s,
    color=color_signal,linestyle=:solid)

    #ssa
    lines!(ax_spec,freq_domain,spec_ssa,linewidth=lw,
    color=color_ssa,linestyle=:solid)
    lines!(ax_spec,freq_domain,spec_res_ssa,linewidth=lw,
    color=color_ssa,linestyle=:solid)

    #nlsa
    lines!(ax_spec,freq_domain,spec_nlsa,linewidth=lw,
    color=color_nlsa,linestyle=:solid)
    lines!(ax_spec,freq_domain,spec_res_nlsa,linewidth=lw,
    color=color_nlsa,linestyle=:solid)

    #hideydecorations!(ax_spec)
    #hidespines!(a_spec, :t, :r, :l) # only top and right
    #axislegend(ax_spec)

    Label(F[0, 1], text = "seasons + FFT for "*titlestring,
    fontsize = 24)

    #trim!(F.layout)

    save(savedir*"seasonal.png",F)

    """
    the run on the years
    """


    year_ind =
    [
    1:365,
    366:730,
    731:1095,
    1096:1460,
    1462:1826, #89 1 offset
    1827:2191, 
    2192:2556,
    2557:2921,
    2923:3287, # 93 1 offset
    3288:3652,
    3653:4017,
    4018:4382,
    4384:4748, #97 1 offset
    4749:5113,
    5114:5478,
    ]

    years_plot = 2005 .+ (1:15)

    color_signal_grad =  cgrad([:grey20,:grey80],15,categorical=true,rev=true)
    color_ssa_grad =  cgrad([:darkgreen,:lightgreen],15,categorical=true,rev=true)
    color_nlsa_grad =  cgrad([:purple,:plum],15,categorical=true,rev=true)

    lw = 1


    fund = hcat([trend_first[i] for i in year_ind]...)

    ssa_harm = hcat([ssa_trend_harm[i] for i in year_ind]...)
    nlsa_harm = hcat([nlsa_trend_harm[i] for i in year_ind]...)


    f = Figure(resolution=(800,800),)

    ax_years = Axis(f[1, 1],
        yticks = ((1:15) ./ 2 .-1,  string.(years_plot)),
        xticks = ([365/2, 400 + 365/2, 800 + 365/2],
        ["fund","SSA harm","NLSA harm"]),
        ylabel="years"
    )


    for i in 1:15
        d = lines!(ax_years,1:365,fund[:,i].+(0.5*i),color = color_signal_grad[i],linewidth=lw)#,offset = i / 4,)
        d = lines!(ax_years,400 .+ (1:365),ssa_harm[:,i].+(0.5*i),color = color_ssa_grad[i],linewidth=lw)
        #d = lines!(800 .+ (1:365),nlsa_fund[:,i].+(0.5*i),color = colors[i])
        d = lines!(ax_years,800 .+ (1:365),nlsa_harm[:,i].+(0.5*i),color = color_nlsa_grad[i],linewidth=lw)
    end

    ax_years_2 = Axis(f[2,1],
    xticks = ([365/2, 400 + 365/2, 800 + 365/2],
        ["fund","SSA harm","NLSA harm"]),
        ylabel="GPP"
        )

    for i in 1:15
        d = lines!(ax_years_2,1:365,fund[:,i],color = color_signal_grad[i],linewidth=lw)#,offset = i / 4,)
        d = lines!(ax_years_2,400 .+ (1:365),ssa_harm[:,i],color = color_ssa_grad[i],linewidth=lw)
        #d = lines!(800 .+ (1:365),nlsa_fund[:,i].+(0.5*i),color = colors[i])
        d = lines!(ax_years_2,800 .+ (1:365),nlsa_harm[:,i],color = color_nlsa_grad[i],linewidth=lw)
    end

    Label(f[0, :], text = "changes in seasonal behavior for "*titlestring,
    fontsize = 24)

    #trim!(f.layout)

    save(savedir*"years.png",f)

end

"""
do it
"""

for spot = [1,3,9]
    savedirname = dir*"W1a/$spot/"
    #mkdir(savedirname)
    plotting_routine_1a(spot,savedirname,20)
end

"""
fluxnet reading scripts
"""

filelist = readdir("/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/")[5:end-3]

#fullset (tower)
#daily resolution
#period >= 20 (from 2020 on)

L = String[]
for file in filelist
    filenamesplit = split(file,"_")
    years = split(filenamesplit[end-1],"-")
    period = parse(Int,years[2])-parse(Int,years[1])
    met = filenamesplit[4]
    res = filenamesplit[5]
    if met == "FULLSET" && res == "DD" && period >= 15
        println(file)
        L = push!(L,String(file))
    end
end

#choose variables
"GPP_DT_VUT_USTAR50","NEE_VUT_USTAR50_DAY"
#create overlap of gapfree variables in 2000-2020
df = CSV.read(datadir*L[1],DataFrame)
n_l = names(df)

for filename in L[1:end-1]
    df = CSV.read(datadir*filename,DataFrame)
    names_i = names(df)
    #2000-2020
    starting = findall(x->x==20000101, df[!,"TIMESTAMP"])[1]
    ending = findall(x->x==20201231, df[!,"TIMESTAMP"])[1]
    # no gaps
    gapless = isempty.([findall(x->x==-9999,df[starting:ending,name]) for name in names(df)])
    names_i = names_i[gapless]
    #inside the large list
    for (i,name) in enumerate(n_l)
        if !(name in names_i)
            deleteat!(n_l,i)
        end
    end
end

# look at GPP
gpp_ind = Int64[]
for (i,n) in enumerate(n_l)

    if occursin("GPP",n)
        println(i,"\t",n)
        gpp_ind = push!(gpp_ind,i)
    end

end

"""
114     GPP_NT_VUT_USTAR50
115     GPP_NT_VUT_MEAN
125     GPP_NT_CUT_USTAR50
126     GPP_NT_CUT_MEAN
158     GPP_DT_VUT_USTAR50
159     GPP_DT_VUT_MEAN
169     GPP_DT_CUT_USTAR50
170     GPP_DT_CUT_MEAN
"""


# look at NEE
nee_ind = Int64[]
for (i,n) in enumerate(n_l)

    if occursin("NEE",n)
        println(i,"\t",n)
        nee_ind = push!(nee_ind,i)
    end

end

"""
32      NEE_VUT_USTAR50
33      NEE_VUT_USTAR50_QC
34      NEE_VUT_MEAN
35      NEE_VUT_MEAN_QC
57      NEE_VUT_USTAR50_NIGHT
58      NEE_VUT_USTAR50_NIGHT_SD
59      NEE_VUT_USTAR50_NIGHT_QC
60      NEE_VUT_USTAR50_DAY
61      NEE_VUT_USTAR50_DAY_SD
62      NEE_VUT_USTAR50_DAY_QC
"""

# look at RECO
reco_ind = Int64[]
for (i,n) in enumerate(n_l)

    if occursin("RECO",n)
        println(i,"\t",n)
        reco_ind = push!(reco_ind,i)
    end

end

"""
92      RECO_NT_VUT_USTAR50
93      RECO_NT_VUT_MEAN
103     RECO_NT_CUT_USTAR50
104     RECO_NT_CUT_MEAN
136     RECO_DT_VUT_USTAR50
137     RECO_DT_VUT_MEAN
147     RECO_DT_CUT_USTAR50
148     RECO_DT_CUT_MEAN
"""

#the big 3:
ind = sort([gpp_ind...,nee_ind...,reco_ind...])
n_l = deleteat!(n_l,ind)

for (i,n) in enumerate(n_l)
    println(i,"\t",n)
end

"""
1       TIMESTAMP
2       TA_F_MDS_QC
3       TA_F_MDS_NIGHT_SD
4       TA_F_MDS_DAY
5       TA_F_MDS_DAY_QC
6       TA_F                        air temperate
7       TA_F_QC
8       TA_F_NIGHT
9       TA_F_NIGHT_SD               standard deviation
10      TA_F_NIGHT_QC               quality flag
11      TA_F_DAY
12      TA_F_DAY_SD
13      TA_F_DAY_QC
14      SW_IN_POT                   shortwave radiation potential top of atmosphere
15      SW_IN_F_MDS_QC
16      SW_IN_F                     shortwave radiation gapfilled
17      SW_IN_F_QC
18      LW_IN_F                     longwave
19      LW_IN_F_QC
20      LW_IN_JSB_QC                longwave calculated from model
21      LW_IN_JSB_F
22      LW_IN_JSB_F_QC
23      VPD_F_MDS_QC                vapour pressure deficit
24      VPD_F
25      VPD_F_QC
26      LE_F_MDS                    latent heat flux
27      H_F_MDS                     sensible heat flux
28      DAY_D
29      DAY_RANDUNC_N
"""


#choose variables
#"NEE_VUT_USTAR50_DAY","RECO_NT_VUT_USTAR50",,"VPD_F","LE_F_MDS","H_F_MDS",]
variables = ["GPP_DT_VUT_USTAR50","TA_F_DAY","SW_IN_F","LW_IN_F","SWC_F_MDS_1","TS_F_MDS_1"]

L = String[]
for file in filelist
    filenamesplit = split(file,"_")
    years = split(filenamesplit[end-1],"-")
    period = parse(Int,years[2])-parse(Int,years[1])
    met = filenamesplit[4]
    res = filenamesplit[5]
    if met == "FULLSET" && res == "DD" && period >= 15
        println(file)
        L = push!(L,String(file))
    end
end

for filename in L
    df = CSV.read(datadir*filename,DataFrame)
    names_i = names(df)
    #check for existing variables
    if all([(vari in names_i) for vari in variables])
        #check for 2005-2020
        starting = findall(x->x==20050101, df[!,"TIMESTAMP"])[1]
        ending = findall(x->x==20201231, df[!,"TIMESTAMP"])[1]
        # no gaps
        if all(isempty.([findall(x->x==-9999,df[starting:ending,name]) for name in variables]))
            println(filename[5:10])
        end
    end

end

"""
BE-Lon
BE-Vie
CH-Lae
CH-Oe2
DE-Geb
DE-Hai
DE-Tha
DK-Sor
FI-Hyy
FR-Aur
IT-BCi
IT-Lav
"""

# these should have it!

#now lets make a spots x variables Float32 table and its good to go

spotslist = [
"BE-Lon",
"BE-Vie",
"CH-Lae",
"CH-Oe2",
"DE-Geb",
"DE-Hai",
"DE-Tha",
"DK-Sor",
"FI-Hyy",
"FR-Aur",
"IT-BCi",
"IT-Lav"
]

#create large spots x variables tensor
#iterate them, fill them

datadir = "/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/"

function choose_incomplete(spot)
    filename_incomplete = "FLX_"*spot*"_FLUXNET2015_FULLSET_DD"
    return findall(x->x[1:33]=="FLX_"*spot*"_FLUXNET2015_FULLSET_DD",filelist)
end

n_spots = length(spotslist)
n_variables = length(variables)
N = Int(floor(365.25 * 16))

large_data_matrix = zeros(Float32,N,n_spots,n_variables)


for (i,spot) in enumerate(spotslist)
    listnumber = choose_incomplete(spot)[1]
    filename = filelist[listnumber]
    df = CSV.read(datadir*filename,DataFrame)
    starting = findall(x->x==20050101, df[!,"TIMESTAMP"])[1]
    ending = findall(x->x==20201231, df[!,"TIMESTAMP"])[1]
    for (j,variable) in enumerate(variables)
        large_data_matrix[:,i,j] = Float32.(df[starting:ending,variable])
        println(spot,"\t",variable)
    end
end
    

savedir="/net/scratch/lschulz/fluxdata_march23/"
savedirname = savedir*"fluxdata_raw.jld2"
jldsave(svedirname,
data=large_data_matrix,spots=spotslist,variables=variables)
