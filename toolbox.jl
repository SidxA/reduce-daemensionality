#----------------------------------------------------------------------
#activate the Pkg
#----------------------------------------------------------------------
using Pkg
Pkg.activate("reduce_daemensionality_julia")
#----------------------------------------------------------------------
#time stamp for saving analysis results
#----------------------------------------------------------------------
using Dates
function init_logging()
    weekno = week(unix2datetime(time()))
    datestring = string("KW_",lpad(weekno,2,"0"),"/")
    workdir = "/net/home/lschulz/logs/"
    dir = workdir*datestring
    if isdir(dir)==false
        mkdir(dir)
    end
    return dir
end
dir = init_logging()
#----------------------------------------------------------------------
#data tools
#----------------------------------------------------------------------
using Statistics
using StatsBase
using MultivariateStats
using Random
using CSV
using DataFrames
using FFTW
using LinearAlgebra
#using Plots
using ManifoldLearning
using Peaks
using NPZ
using FourierAnalysis
#using Plots
using FFTW
using EmpiricalModeDecomposition
using ProgressMeter
using JLD2
using ProfileSVG
using LsqFit
using CairoMakie
using MathTeXEngine
using FileIO
using ImageIO
using NearestNeighbors
using SharedArrays

""" signal stuff"""
#----------------------------------------------------------------------
#centralize
#embed lag
#reconstructor

#----------------------------------------------------------------------

    #centralize
    function centralizer(data::Vector{Float32})
        m = mean(data)
        cv = std(data)
        return  (data.-m)./cv
    end


    #embed lag
    function embed_lag(data::Vector{Float32},W::Int64)
        Y = []
        for i=1:length(data)-W+1
            Y = append!(Y,data[i:i+W-1])
        end
        return reshape(float.(Y),W,length(data)-W+1)::Matrix{Float32}
    end

    #reconstruction ssa
    function reconstructor(A,rho,N,M)
        P = N-M+1
        #A is the k_th PC (projection) of length P
        #rho is the k_th eigenvector (EOF) of length M

        #M is window length W
        R(t,M_t,L_t,U_t) = M_t*sum([A[t-j+1]*rho[j] for j in L_t:U_t])
        function choicer(t)
            if 1<=t<=M-1
                return 1/t,1,t
            elseif M<=t<=P
                return 1/M,1,M
            elseif P+1<=t<=N
                return 1/(N-t+1),t-N+M,M
            end
        end
        return [R(t,choicer(t)...) for t in 1:N]
    end






""" proto stuff"""
#----------------------------------------------------------------------
# protophase        signal , num_it
# protofrequency    signal, year=samples per year
# count sign flips
# count maxima
#iterated hilbert
#----------------------------------------------------------------------


    #count sign flips
    count_sign_flip(signal) = sum([(sign.(signal[i])!=sign.(signal[i+1])) ? 1 : 0 for i in 1:length(signal)-1])
    #count maxima
    count_maxima(signal) = length(findmaxima(signal)[1])


    function iterated_hilbert(mode::Vector,l)
        ll = Int64[]
        for i = 1:l
            HilbertT = hilbert_transform(Float64.(mode))
            protophases = atan.(imag(HilbertT),real(HilbertT))
            s = count_sign_flip(protophases)
            #println(s)
            ll = append!(ll,s)
            mode = protophases
        end
        return median(ll),var(ll)/sqrt(l)
    end


    #function that spits out the indeces of the modes according to the domains in all locations as a matrix of lists [bandxloc]
    function protofrequency(signal::Vector{Float32},yearsamples,Nit)
        zeros = []
        for i in 1:Nit
            HilbertT = hilbert_transform(Float64.(signal))
            protophases = atan.(imag(HilbertT),real(HilbertT))
            zeros = append!(zeros,count_sign_flip(protophases))
            signal = protophases
        end
        T = 2 * length(signal) / median(zeros)
        return Float32(yearsamples / T )
    end

    function protophase(signal::Vector{Float32},Nit)
        protoph = atan.(imag(hilbert_transform(Float64.(signal))),real(hilbert_transform(Float64.(signal))))
        for i=1:Nit-1
            protoph = atan.(imag(hilbert_transform(Float64.(protoph))),real(hilbert_transform(Float64.(protoph))))
        end
        return protoph
    end



""" proto stuff"""
#----------------------------------------------------------------------
#entropy signal
#correlation_coefficient signal1 signal2
#rec_phase_complete RC::M, lam::v
#rec_protophase RC::M Num_it        <==
#phase_diff_hist mode1,mode2 gives bins,cauchy_weights
#----------------------------------------------------------------------


    #entropy
    function entropy(signal::Vector{Float32})
        signal = abs.(centralizer(signal))
        return -sum(signal .* log.(signal))
    end

    #crosscorrelation coefficient
    function correlation_coefficient(x::Vector{Float32},y::Vector{Float32})
        #lags = 1:Int(length(x)-1)
        #c = crosscor(x,y,lags)
        return cov(x,y) / var(x) / var(y)
    end


    function rec_phase_complete(RC::Matrix{Float32},lam::Vector{Float32})
        ph = [protophase(RC[:,i]) for i=1:size(RC)[2]]
        pha = hcat(ph...)
        phas = vec(sum(pha'.*lam,dims=1))
        return centralizer(phas)
    end

    rec_protophase(RC::Matrix{Float32},Nit) = Float32.(hcat([protophase(RC[:,i],Nit) for i=1:size(RC)[2]]...))

    function phase_diff_hist(mode1::Vector{Float32},mode2::Vector{Float32}) #takes in protophases
        bins = (-pi:0.01:pi)
        phasediff = mode1 .- mode2
        h = fit(Histogram,phasediff,bins)
        hist = normalize(h,mode=:pdf).weights
        p0 = ones(2)
        c = coef(curve_fit(cauchy,Array(bins)[1:end-1],hist,p0))
        return bins[1:end-1],hist,c[1],c[2]
    end

"""fluxnet simple measures"""
#----------------------------------------------------------------------
#spec entropy
#diff eps
#autocorr
#rec_protophase RC::M Num_it        <==
#phase_diff_hist mode1,mode2 gives bins,cauchy_weights
#----------------------------------------------------------------------

    #fluxnet entropy based on lambda
    function spectral_entropy(datadir,var_ind,loc_ind,rec_m,season_ind)
        if rec_m == "ssa"
            lambda = load(datadir*"ssa_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["lambda"]
        elseif rec_m == "diff"
            lambda = load(datadir*"diff_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["lambda"]
        end
        return entropy(lambda)
    end

    #diffusion eps
    function diff_eps(datadir,var_ind,loc_ind)
        season_ind = 1
        eps = load(datadir*"diff_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["eps"]
        return eps
    end

    #autocorrelation strength complete reconstructed signal
    function autocorr(datadir,var_ind,loc_ind,rec_m)
        season_ind = 1
        if rec_m == "ssa"
            RC = load(datadir*"ssa_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"]
        elseif rec_m == "diff"
            RC = load(datadir*"diff_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"]
        end
        return correlation_coefficient(sum(RC,dims=2)[:,1],sum(RC,dims=2)[:,1])
    end

    #autocorrelation strength raw signal
    function autocorr_raw(datadir,var_ind,loc_ind)
        season_ind = 1
        signal = load(datadir*"diff_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["signal"]
        return correlation_coefficient(signal,signal)
    end



    """
    coupling whole signal
        cauchy fits
    """
    #readout signal and perform single hilbert transform
    function protophase_by_varind_recm(datadir,loc_ind,season_ind,var_ind,rec_m)
        if rec_m == "raw"
            signal = load(datadir*"ssa_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["signal"]
        elseif rec_m == "ssa"
            signal = sum(load(datadir*"ssa_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"],dims=2)[:,1]
        elseif rec_m == "diff"
            signal = sum(load(datadir*"diff_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"],dims=2)[:,1]
        end
        return atan.(imag(hilbert_transform(Float64.(signal))),real(hilbert_transform(Float64.(signal))))
    end

    function cauchy(x,p) # fit cauchy / delta to a normalized phase difference
        x0 = p[1]
        gamma = p[2]
        @. return 1/(pi*gamma*(1+((x-x0)/gamma)^2))
    end

    function fit_cauchy(onebins) # perform the fit to the histogram list
        bins = (-pi:0.01:pi)
        bin_N = length(bins)-1
        p0 = ones(2)
        return coef(curve_fit(cauchy,Array(bins)[1:end-1],onebins,p0))
    end

    #complete signal phase difference histogram between two spots
    function varmap(loc_ind_1,loc_ind_2,datadir,season_ind,rec_m,var_ind)
        bins = (-pi:0.01:pi)
        bin_N = length(bins)-1
        phasediff = protophase_by_varind_recm(datadir,loc_ind_1,season_ind,var_ind,rec_m) .- protophase_by_varind_recm(datadir,loc_ind_2,season_ind,var_ind,rec_m)
        h = fit(Histogram,phasediff,bins)
        h = normalize(h,mode=:pdf).weights
        return h
    end

    #cauchy fit between two spots for complete reconstructior or raw signal
    function phase_diff_complete(datadir,loc_ind_1,loc_ind_2,rec_m,var_ind)
        season_ind = 1
        hist = varmap(loc_ind_1,loc_ind_2,datadir,season_ind,rec_m,var_ind)
        cauchypar = fit_cauchy(hist)
        x0 = cauchypar[1]
        gamma = cauchypar[2]
        return x0,gamma
    end

    """
    coupling index based
        cauchy fits
    """

    #index based readout signal and perform single hilbert transform
    function protophase_by_varind_recm_modeind(datadir,loc_ind,mode_ind,season_ind,var_ind,rec_m)
        if rec_m == "ssa"
            signal = sum(load(datadir*"ssa_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"][:,mode_ind],dims=2)[:,1]
        elseif rec_m == "diff"
            signal = sum(load(datadir*"diff_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"][:,mode_ind],dims=2)[:,1]
        end
        return atan.(imag(hilbert_transform(Float64.(signal))),real(hilbert_transform(Float64.(signal))))
    end

    #index based signal phase difference histogram between two spots
    function varmap_ind(loc_ind_1,loc_ind_2,mode_ind_1,mode_ind_2,datadir,season_ind,rec_m,var_ind)
        bins = (-pi:0.01:pi)
        bin_N = length(bins)-1
        phasediff = protophase_by_varind_recm_modeind(datadir,loc_ind_1,mode_ind_1,season_ind,var_ind,rec_m)
            .- protophase_by_varind_recm_modeind(datadir,loc_ind_2,mode_ind_2,season_ind,var_ind,rec_m)
        h = fit(Histogram,phasediff,bins)
        h = normalize(h,mode=:pdf).weights
        return h
    end

    #index of Tband based protophase difference cauchy fit # needs to be tested not to contain nothing
    function phase_diff_ind(datadir,loc_ind_1,loc_ind_2,mode_ind_1,mode_ind_2,rec_m,var_ind)
        season_ind = 1
        bins = (-pi:0.01:pi)
        bin_N = length(bins)-1
        hist = varmap_ind(loc_ind_1,loc_ind_2,mode_ind_1,mode_ind_2,datadir,season_ind,rec_m,var_ind)
        cauchypar = fit_cauchy(hist)
        x0 = cauchypar[1]
        gamma = cauchypar[2]
        return x0,gamma
    end

    """
    coupling
        crosscorrelation coefficient index based
    """

    #index based readout signal
    function reconstruction_by_varind_recm_modeind(datadir,loc_ind,mode_ind,season_ind,var_ind,rec_m)
        if rec_m == "ssa"
            signal = sum(load(datadir*"ssa_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"][:,mode_ind],dims=2)[:,1]
        elseif rec_m == "diff"
            signal = sum(load(datadir*"diff_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"][:,mode_ind],dims=2)[:,1]
        end
        return signal
    end

    #index based correlation crosscorrelation coefficient - can also be used for autocorrelation
    function crosscor_ind(loc_ind_1,loc_ind_2,mode_ind_1,mode_ind_2,datadir,season_ind,rec_m,var_ind)
        season_ind = 1
        signal_1 = reconstruction_by_varind_recm_modeind(datadir,loc_ind_1,mode_ind_1,season_ind,var_ind,rec_m)
        signal_2 = reconstruction_by_varind_recm_modeind(datadir,loc_ind_2,mode_ind_2,season_ind,var_ind,rec_m)
        return correlation_coefficient(signal_1,signal_2)
    end



"""the heavy weights"""

function phase_rec_fluxer(datadir,loc_ind,mode_indices,rec_m,var_ind)
    RC = load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"][:,mode_indices]
    lam = [norm(load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["PC"][:,m]) for m=mode_indices]
    return rec_phase_complete(RC,lam)
end

function mode_protofrequencies(datadir::String,loc_ind,mode_indices,rec_m,var_ind)
    RC = load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"][:,mode_indices]
    lam = [norm(load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["PC"][:,m]) for m=mode_indices]
    return [protofrequency(RC[:,i]) for i=1:size(RC)[2]],lam
end

function all_mode_protofrequencies(datadir::String,rec_m,var_ind,loc_indices,mode_indices,season_ind)
    f = Matrix{Float32}(undef,length(mode_indices),length(loc_indices))
    l = Matrix{Float32}(undef,length(mode_indices),length(loc_indices))
    for loc_ind in loc_indices
        RC = load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"][:,mode_indices]
        lam = [norm(load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["PC"][:,m]) for m=mode_indices]
        f[:,loc_ind] = [protofrequency(RC[:,i]) for i=1:size(RC)[2]]
        l[:,loc_ind] = lam
    end
    return f,l
end

function reconstruction_error(datadir,var_ind,loc_ind,rec_m,season_ind)
    original = load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_3.jld2")["signal"]
    deseason = load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["signal"]
    season = original .- deseason
    rec_deseas = vec(sum(load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"],dims=2))
    return norm(original.-season),norm(deseason .- rec_deseas)
end

function phase_diff_spots_by_modeind(datadir,loc_ind_1,loc_ind_2,mode_indices_1,mode_indices_2,rec_m,var_ind,season_ind)
    bins = (-pi:0.01:pi)
    bin_N = length(bins)-1
    phasediff = phase_rec_fluxer(datadir,loc_ind_1,mode_indices_1,rec_m,var_ind) .- phase_rec_fluxer(datadir,loc_ind_2,mode_indices_2,rec_m,var_ind)
    h = fit(Histogram,phasediff,bins)
    hist = normalize(h,mode=:pdf).weights
    cauchypar = fit_cauchy(hist)
    x0 = cauchypar[1]
    gamma = cauchypar[2]
    return x0,gamma
end

function crosscorr_spots_by_modeind(datadir,loc_ind_1,loc_ind_2,mode_indices_1,mode_indices_2,rec_m,var_ind,season_ind)
    if !isempty(mode_indices_1) && !isempty(mode_indices_2)
        RC1 = vec(sum(load(datadir*"$(rec_m)_$(loc_ind_1)_$(var_ind)_$(season_ind).jld2")["RC"][:,mode_indices_1],dims=2))
        RC2 = vec(sum(load(datadir*"$(rec_m)_$(loc_ind_2)_$(var_ind)_$(season_ind).jld2")["RC"][:,mode_indices_2],dims=2))
        return correlation_coefficient(RC1,RC2)
    else
        return 0.0
    end
end
#cauchy
function coupling_strength_cauchy(coupling_inds,mode_inds,datadir,rec_m,var_ind,season_ind)
    function strengthperind(inds)
        loc_ind_1 = inds[1]
        loc_ind_2 = inds[2]
        mode_indices_1 = mode_inds[loc_ind_1]
        mode_indices_2 = mode_inds[loc_ind_2]
        if !isempty(mode_indices_1) && !isempty(mode_indices_2)
            return phase_diff_spots_by_modeind(datadir,loc_ind_1,loc_ind_2,mode_indices_1,mode_indices_2,rec_m,var_ind,season_ind)[2]
        else
            return 0.0
        end
    end
    m = map(strengthperind,coupling_inds)
    return m
end

function coupling_phase_cauchy(coupling_inds,mode_inds,datadir,rec_m,var_ind,season_ind)
    function strengthperind(inds)
        loc_ind_1 = inds[1]
        loc_ind_2 = inds[2]
        mode_indices_1 = mode_inds[loc_ind_1]
        mode_indices_2 = mode_inds[loc_ind_2]
        if !isempty(mode_indices_1) && !isempty(mode_indices_2)
            return phase_diff_spots_by_modeind(datadir,loc_ind_1,loc_ind_2,mode_indices_1,mode_indices_2,rec_m,var_ind,season_ind)[1]
        else
            return 0.0
        end
    end
    m = map(strengthperind,coupling_inds)
    return m
end

#crosscorr coefficient for coupling strength
function coupling_strength_crosscorr(coupling_inds,mode_inds,datadir,rec_m,var_ind,season_ind)

    function strengthperind(inds)
        loc_ind_1 = inds[1]
        loc_ind_2 = inds[2]
        mode_indices_1 = mode_inds[loc_ind_1]
        mode_indices_2 = mode_inds[loc_ind_2]
        return crosscorr_spots_by_modeind(datadir,loc_ind_1,loc_ind_2,mode_indices_1,mode_indices_2,rec_m,var_ind,season_ind)
    end
    m = map(strengthperind,coupling_inds)
    return m
end


"""frequency bands"""

function list_harmonics(hmax)
    L = Float32[]
    for h=1:hmax
        for l=1:hmax
            L = append!(L,l/h)
            L = append!(L,h/l)
        end
    end
    return sort(unique(L))
end

function list_pure_harmonics(hmax)
    L = Float32[]
    for h=1:hmax
        L = append!(L,1/h)
        L = append!(L,h)
    end
    return sort(unique(L))
end

harmonic_bands(hmax,eps1) = [[i-eps1,i+eps1] for i in list_harmonics(hmax)]

anharmonic_bands(har_bands,eps2) = [[har_bands[i][2]-eps2,har_bands[i+1][1]+eps2] for i = 1:length(har_bands)-1]
simplebands(eps) = [[10^-3,0.375+eps],[0.375-eps,1.25+eps],[1.25-eps,10^3]]
simplebands() = [[10^-3,0.375],[0.375,1.25],[1.25,10^3]]

function Tbands(freq_domains,var_ind,loc_inds,datadir,rec_m,season_ind)
    ind = Matrix{Vector{Int64}}(undef,length(freq_domains),length(loc_inds))
    lambi = Matrix{Vector{Float32}}(undef,length(freq_domains),length(loc_inds))
    for (i,loc_ind) in enumerate(loc_inds), (j,domain) in enumerate(freq_domains)
        f_low = domain[1]
        f_up = domain[2]
        mode_indices=1:48
        p_frequencies,lam = mode_protofrequencies(datadir::String,loc_ind,mode_indices,rec_m,var_ind)
        #mode_ind = findall(freq->f_up>freq>f_low,p_frequencies)
        mode_ind = Array(1:48)[f_low .< p_frequencies .< f_up]
        ind[j,i] = vec(mode_ind)
        lambi[j,i] = lam[vec(mode_ind)]
    end
    return ind,lambi
end

function matrix_Tbands(freq_domains,protofrequencies)
    ind = Matrix{Vector{Int64}}(undef,length(freq_domains),size(protofrequencies)[2])
    for i in 1:size(protofrequencies)[2], (j,domain) in enumerate(freq_domains)
        f_low = domain[1]
        f_up = domain[2]
        p_frequencies = protofrequencies[:,i]
        mode_ind = Array(1:48)[f_low .< p_frequencies .< f_up]
        ind[j,i] = Int64.(vec(mode_ind))
    end
    return ind
end

function calculate_bands(datadir,rec_m,var_ind,season_ind)



    f,l = all_mode_protofrequencies(datadir,rec_m,var_ind,1:71,1:48,season_ind)

    b1 = simplebands(0.02)

    b2h = harmonic_bands(8,0.02)
    b2a = anharmonic_bands(b2h,0.01)
    b2 = vcat(b2h,b2a)

    i1 = matrix_Tbands(b1,f)
    i2 = matrix_Tbands(b2,f)

    B1 = [Int64.(vcat([vcat([Int64(j) for j in i]...) for i in unique(i2[1:43,k])]...)) for k in 1:71]
    B2 = [Int64.(vcat([vcat([Int64(j) for j in i]...) for i in unique(i2[44:end,k])]...)) for k in 1:71]

    return l,[i1[1,:],i1[2,:],i1[3,:]],[B1,B2]

end


function extract_from_files(filename_list,N=11322,k = 48) #with .jld2 both in files and savename
    year=365
    L = length(filename_list)
    
    lambda = Array{Float32}(undef,k,L)
    protophases = Array{Float32}(undef,N,k,L)
    protofreq = Array{Float32}(undef,k,L)
    RC = Array{Float32}(undef,N,k,L)
    for (i,filename) in enumerate(filename_list)
        File = load(filename)
        lambda[:,i] = File["lambda"]
        RC[:,:,i] = File["RC"]
        protophases[:,:,i] = rec_protophase(RC[:,:,i],1)
        protofreq[:,i] = hcat([protofrequency(protophases[:,kk,i],year) for kk in 1:k]...)
    end
    return lambda,protofreq,RC
    
end #return lambda(k,L),protophases(N,k,L),protofreq(k,L),RC(N,k,L)

function extract_from_single_file(filename,yearsamples,N=11322,k = 48)
    File = load(filename)
    lambda = File["lambda"]
    EOF = File["EOF"]
    RC = File["RC"]
    protophases = rec_protophase(RC,1)
    protofreq = hcat([protofrequency(EOF[:,kk],yearsamples,1) for kk in 1:k]...)
    return lambda,protofreq,RC
end


function create_file_list(outdir,method,W,vari,preproc) #at the moment this does diffusion
    allfiles = readdir(outdir)
    filenamesplit = split.(allfiles,"_")
    method_inds = findall(x->( x[1]==method && x[2] == "$W" && x[4] == "$vari" && x[end][1:4] == preproc),[filenamesplit[i] for i in 1:length(filenamesplit)])
    files = outdir .* allfiles[method_inds]
    return files
end #return files

"""
analyze the phase-coupling of two modes
"""
function phase_synchronization(mode1::Vector{Float32},mode2::Vector{Float32}) #takes in protophases
    bins = (-pi:0.01:pi)
    phasediff = mode1 .- mode2
    h = fit(Histogram,phasediff,bins)
    hist = normalize(h,mode=:pdf).weights
    p0 = ones(2)
    c = coef(curve_fit(cauchy,Array(bins)[1:end-1],hist,p0))
    return Float32.([max(min(c[1],10),-10),c[2]])
end #return x0,gamma

"""
analyze phase coupling along mode_indices inside a file where the bands are already put together
"""
function cauchy_by_modes(protobands)
    T,spots,bands = size(protobands)
    gamma = ones(Float32,spots,spots,bands,bands)
    x0 = ones(Float32,spots,spots,bands,bands)

    for band1 in 1:bands, spot1 in 1:spots, band2 in 1:bands, spot2 in 1:spots
        if spot1 == spot2
            g,x = 0,0
        else
            g,x = phase_synchronization(protobands[:,spot1,band1],protobands[:,spot2,band2])
        end
        gamma[spot1,spot2,band1,band2] = g
        x0[spot1,spot2,band1,band2] = x
    end

    #gamma[findall(x->x>10,gamma)] = 10 
    #x0[findall(x->x>10,x0)] = 10 
    return Float32.(gamma), Float32.(x0)
    
end #return gamma(L,L,bands,bands) , x0(L,L,bands,bands)

"""
put the modes together by band
"""
function band_indices(freq_domains,protofrequencies) #read in a matrix of the frequencies
    #need some feeding of locking up all the protofrequencies of all the spots we want
    #should be easily done
    
    n_bands = length(freq_domains)
    n_spots =size(protofrequencies)[2]
    ind = Matrix{Vector{Int64}}(undef,n_bands,n_spots)

    for i in 1:n_spots, (j,domain) in enumerate(freq_domains)
        f_low = domain[1]
        f_up = domain[2]
        p_frequencies = protofrequencies[:,i]
        mode_ind = Array(1:48)[f_low .< p_frequencies .< f_up]
        ind[j,i] = Int64.(vec(mode_ind))
    end
    return ind
end #return ind (bands, L) vec

function combine_by_bands(indices,RC,lambdas)
    T, n_comp, n_spots = size(RC)
    n_bands = size(indices)[1]
    #then we iterate over all of them getting to know the modes we want to have
    #then one large file of the modes! which is only N bands spots shoudl be possible for all data size, this is great!
    protobands = Array{Float32}(undef,T,n_spots,n_bands)
    lambda_bands = Array{Float32}(undef,n_spots,n_bands)

    for band in 1:n_bands, spot in 1:n_spots
        mode_indices = indices[band,spot]
        protobands[:,spot,band] = Float32.(10^-5 .+ protophase(sum(RC[:,mode_indices,spot],dims=2)[:],1))
        lambda_bands[spot,band] = Float32.(sum(10^-5 .+ lambdas[mode_indices,spot]))
        if isempty(mode_indices)
            protobands[:,spot,band] = ones(Float32,T)*10^-5
        end
    end
    return protobands,lambda_bands
end #return protobands(N,L,bands),lambda_bands(N,L,bands)


#----------------------------------------------------------------------
# modelling stuff
#----------------------------------------------------------------------

