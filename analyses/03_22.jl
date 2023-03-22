"""
put the meta in a file

took out FR-Aur because there is no meta file
jldsave("/net/scratch/lschulz/fluxfullset_march23/meta.jld2",data=meta[spots_numbers_list,:])

"""
spotslist = [
"BE-Lon","BE-Vie","CH-Lae","CH-Oe2","DE-Geb","DE-Hai","DE-Tha","DK-Sor","FI-Hyy","IT-BCi","IT-Lav"
]

spots_numbers_list = [
28,29,56,58,73,75,85,89,94,108,117,
]

variables = ["GPP_DT_VUT_USTAR50","TA_F_DAY","SW_IN_F","LW_IN_F","SWC_F_MDS_1","TS_F_MDS_1"]


"""
required stuff
"""

include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")
using CairoMakie

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

"""
plotting function
    for individual W,spot,variable
"""

#main plotting parameters
outdir="/net/scratch/lschulz/fluxfullset_march23/"
meta = load("/net/scratch/lschulz/fluxfullset_march23/meta.jld2")["data"]
variables_list = variables
L = N = 5844
startyear = 2004
preproc = "raw." 
kappa= 48
yearsamples=365.25
bins = 10 .^Array(-2:0.01:2)
Ts = 1 / 365.25
t0 = 0 
tmax = t0 + (N-1) * Ts
t = t0:Ts:tmax
tw = t0:Ts:(t0+(W-1)*Ts)

color_ssa = "darkgreen"
color_nlsa = "purple"
color_signal = "grey50"
lw = 3
ms = 5

years = ((1:L) ./ 365) .+ startyear

# first one: create reconstruction picture
function plot_reconstruction(F,time_series_arguments)
    spot,W,vari = time_series_arguments

    """
    first one signal rec + fft
    """

    #overkill?
    Filename_ssa = create_file_list(outdir,"ssa",W,vari,preproc)[spot]
    Filename_nlsa = create_file_list(outdir,"diff",W,vari,preproc)[spot]
    file_ssa = load(Filename_ssa)
    rc_ssa = sum(file_ssa["RC"],dims=2)[:]
    file_nlsa = load(Filename_nlsa)
    rc_nlsa = sum(file_nlsa["RC"],dims=2)[:]
    signal = file_ssa["signal"]
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

    #figure

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

    #Label(F[0, 1], text = "reconstructions for "*titlestring,
    #fontsize = 24)

    #trim!(F.layout)

    return F
end

#second one: create modes picture
function plot_modeshapes(F,time_series_arguments)
    spot,W,vari = time_series_arguments

    """
    EOF + FFT
    """


    #overkill?
    Filename_ssa = create_file_list(outdir,"ssa",W,vari,preproc)[spot]
    Filename_nlsa = create_file_list(outdir,"diff",W,vari,preproc)[spot]
    file_ssa = load(Filename_ssa)
    rc_ssa = sum(file_ssa["RC"],dims=2)[:]
    file_nlsa = load(Filename_nlsa)
    rc_nlsa = sum(file_nlsa["RC"],dims=2)[:]
    signal = file_ssa["signal"]
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

    #new fft stuff for the W not the N !

    xtic=1:6
    years_modes = (1:W) ./ 365
    freqstart = Int(round(W/2,digits=0))+1
    tw = t0:Ts:(t0+(W-1)*Ts)
    spot=1




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

    #Label(F[0, 1:2], text = "EOF with FFT for "*titlestring,
    #fontsize = 24)

    #trim!(F.layout)



    return F
end

#third one: create seasonality picture
function plot_seasonality(F,time_series_arguments)
    spot,W,vari = time_series_arguments

    """
    seasons + FFT
    """


    #overkill?
    Filename_ssa = create_file_list(outdir,"ssa",W,vari,preproc)[spot]
    Filename_nlsa = create_file_list(outdir,"diff",W,vari,preproc)[spot]
    file_ssa = load(Filename_ssa)
    rc_ssa = sum(file_ssa["RC"],dims=2)[:]
    file_nlsa = load(Filename_nlsa)
    rc_nlsa = sum(file_nlsa["RC"],dims=2)[:]
    signal = file_ssa["signal"]
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

    file = load(Filename_nlsa)
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

    #Label(F[0, 1], text = "seasons + FFT for "*titlestring,
    #fontsize = 24)

    #trim!(F.layout)


    return F
end

#fourth one: create annual comparison picture
function plot_annuality(F,time_series_arguments)
    spot,W,vari = time_series_arguments

    """
    the run on the years
    """


    #overkill?
    Filename_ssa = create_file_list(outdir,"ssa",W,vari,preproc)[spot]
    Filename_nlsa = create_file_list(outdir,"diff",W,vari,preproc)[spot]
    file_ssa = load(Filename_ssa)
    rc_ssa = sum(file_ssa["RC"],dims=2)[:]
    file_nlsa = load(Filename_nlsa)
    rc_nlsa = sum(file_nlsa["RC"],dims=2)[:]
    signal = file_ssa["signal"]
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

    file = load(Filename_nlsa)
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




    ax_years = Axis(F[1, 1],
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

    ax_years_2 = Axis(F[2,1],
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

    #Label(F[0, :], text = "changes in seasonal behavior for "*titlestring,
    #fontsize = 24)

    #trim!(f.layout)


    return F
end



function plotter(spot,
savedirname,
W = 2556,
vari = 1,
)

    !isdir(savedirname) ? mkdir(savedirname) : "go" 


    # plotting parameters


    igbp = meta[spot,"IGBP_class"]
    name = meta[spot,"name"]
    variable = variables_list[vari]

    titlestring = "$variable in $spot $name $igbp"
    savedirstring = "$(variable)_$(spot)_"
    savedir = savedirname*savedirstring

    time_series_arguments = [spot,W,vari]


    F = Figure(resolution=(1600,1600))

    #first one
    
    plot_reconstruction(F[1,1],time_series_arguments)

    #second one
    plot_modeshapes(F[1,2],time_series_arguments)

    #third one
    plot_seasonality(F[2,1],time_series_arguments)

    #fourth one
    plot_annuality(F[2,2],time_series_arguments)


    save(savedir*".png",F)

end

#for spot = 1:11, vari = 1:6, W = Int.(floor.(365 .* [2,3,4,5,6,7,8]))
for spot = 1:11,W = Int.(floor.(365 .* [2,3,4,5,6,7,8]))
    plotter(
    spot,
    dir*"test_$W/",
    W,
    vari,
    )
end


"""
spot
spotname
IGBPclass
W
L
variablesnumber
variablename

signal
ssa_eof
ssa_pc
ssa_rc
ssa_lambda
nlsa_eof
nlsa_pc
nlsa_rc
nlsa_lambda
"""
