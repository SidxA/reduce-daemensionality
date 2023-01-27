


"""
little figure for the frequency bands
"""
function mode_times()
    F = Figure(resolution=(1800,900))

    #all freq

    rec_m = "ssa"
    var_ind=1
    loc_ind=rand(vec(1:71))
    season_ind = rand(vec(1:3))

    #this does not work yet for full times
    allfreq,lam = [mode_protofrequencies(datadir,loc_ind,1:48,rec_m,var_ind) for loc_ind in 1:71]
    freq_ind = sortperm(allfreq)

    ax=Axis(F[1,1],xlabel = "log10 freq",ylabel = "variance")
    xlims!(ax,(-2,2))
    barplot!(ax,log10.(sort(allfreq)),lam[freq_ind],color=:black,strokewidth=1)


    save(dir*"test.png",F)
end

#gamma
function heatmap_coupling_by_band(name,mode_inds,mode_lam,var_ind,datadir,rec_m,season_ind)
    loc_strength = [sum(mode_lam[:,k][mode_inds[k]]) for k=1:71]
    loc_str_ind = sortperm(loc_strength)

    coupling_inds = list_all_links(loc_str_ind)
    #mode inds come by band in the form of a vector of length locations with the indices as a vector inside
    c_cauch = coupling_strength_cauchy(coupling_inds,mode_inds,datadir,rec_m,var_ind,season_ind)
    c_cauchy = reshape(c_cauch,71,71)
    c_cauchy_var = [c_cauchy[i,j] * loc_strength[loc_str_ind[i]] * loc_strength[loc_str_ind[i]] for i=1:71,j = 1:71]
    c_cross = coupling_strength_crosscorr(coupling_inds,mode_inds,datadir,rec_m,var_ind,season_ind)
    c_cros = reshape(c_cross,71,71)
    F = Figure(resolution=(1200,900))
    colorlimits=(0,2)
    ax1 = Axis(F[1,1])
    hm1 = heatmap!(ax1,c_cauchy,colormap=:Spectral)#,colorrange=colorlimits)
    #ax2 = Axis(F[2,1])
    #hm2 = heatmap!(ax2,c_cauchy_var,colormap=:Spectral)#,colorrange=colorlimits)
    #ax3 = Axis(F[3,1])
    #hm3 = heatmap!(ax3,c_cros,colormap=:Spectral)#,colorrange=colorlimits)
    Colorbar(F[1,end+1],hm1)
    #Colorbar(F[2,end+1],hm2)
    #Colorbar(F[3,end+1],hm3)

    save(dir*"$(name).png",F)
end


function perform_gammacoupling()#ssa taken out
    var_ind=1
    datadir="/net/scratch/lschulz/fluxdata/"

    for rec_m = ["diff"], season_ind = 1:4
        mode_lam, simple, harm = calculate_bands(datadir,rec_m,var_ind,season_ind)
        modilist = [simple[1],simple[2],simple[3]]
        for (i,band_name) in enumerate(["slow","annualseasonal","fast"])#,"harmonic","anharmonic"])
            mode_inds = modilist[i]
            name = "VAR_$(rec_m)_seas$(season_ind)_$(band_name)"
            coupling_inds = list_all_links(loc_str_ind)
            c_cauch = coupling_strength_cauchy(coupling_inds,mode_inds,datadir,rec_m,var_ind,season_ind)
            x_cauch = coupling_phase_cauchy(coupling_inds,mode_inds,datadir,rec_m,var_ind,season_ind)
            F = Figure(resolution=(900,1800))
            colorlimits=(0,2)
            ax1 = Axis(F[1,1])
            hm1 = heatmap!(ax1,reshape(c_cauch,71,71),colormap=:Spectral)
            Colorbar(F[1,end+1],hm1)
            ax2 = Axis(F[2,1])
            hm2 = heatmap!(ax2,reshape(x_cauch,71,71),colormap=:Spectral)
            Colorbar(F[2,end+1],hm2)

            save(dir*"$(name).png",F)
            println(name)
        end
    end
end

function perform_gammacoupling_by_specent()#ssa taken out
    var_ind=1
    datadir="/net/scratch/lschulz/fluxdata/"

    for rec_m = ["diff"], season_ind in 1:4
        ent_ind = sortperm([spectral_entropy(datadir,var_ind,k,rec_m,season_ind) for k=1:71])
        coupling_inds = list_all_links(ent_ind)
        mode_lam, simple, harm = calculate_bands(datadir,rec_m,var_ind,season_ind)
        modilist = [simple[1],simple[2],simple[3]]
        for (i,band_name) in enumerate(["slow","annualseasonal","fast"])#,"harmonic","anharmonic"])
            mode_inds = modilist[i]
            name = "ENT_$(rec_m)_seas$(season_ind)_$(band_name)"
            c_cauch = coupling_strength_cauchy(coupling_inds,mode_inds,datadir,rec_m,var_ind,season_ind)
            x_cauch = coupling_phase_cauchy(coupling_inds,mode_inds,datadir,rec_m,var_ind,season_ind)
            F = Figure(resolution=(900,1800))
            colorlimits=(0,2)
            ax1 = Axis(F[1,1])
            hm1 = heatmap!(ax1,reshape(c_cauch,71,71),colormap=:Spectral)
            Colorbar(F[1,end+1],hm1)
            ax2 = Axis(F[2,1])
            hm2 = heatmap!(ax2,reshape(x_cauch,71,71),colormap=:Spectral)
            Colorbar(F[2,end+1],hm2)

            save(dir*"$(name).png",F)
            println(name)
        end
    end
end


#the required stuff
begin
        
    """
    if we do different datasets we need to formalize this approach to one that is feasable also for large datasets!
        for this we need to be sure about the lambda!
        need to do some testing...
        maybe just build the read-in into one file for now, this might later get a split along the first axis or something
        
    """
    #for 70 spots 48 modes this is 10 MB!


    function coupling_analysis(outdir,W,vari,method)
        #files from the diff,ssa runs
        #outdir = "/net/scratch/lschulz/ta_erai_11a/"
        #savedir = dir

        #one method, one seasonality, one variable
        filename_list = create_file_list(outdir,method,W,vari,preproc)

        #extract protostuff
        #protoname = savedir*"ta_erai_11a_proto.jld2"
        lambda,protofreq,RC = extract_from_files(filename_list,5478,48)
        #jldsave(protoname,lambda = lambda,protofreq = protofreq, RC = RC)

        #simple bands ( with strict boundaries )
        freq_domains = simplebands()
        #bands_name = savedir*"ta_erai_11a_simplebands.jld2"
        indices = band_indices(freq_domains,protofreq)
        protobands,lambda_bands = combine_by_bands(indices,RC,lambda)

        
        """
        jldsave(bands_name,ind=ind, protobands=protobands,lambda_bands = lambda_bands)

        #coupling in the simple bands
        coupling_name = dir*"ta_erai_11a_simplebands_coupling.jld2"
        gamma, x0 = cauchy_by_modes(protobands)
        jldsave(coupling_name,gamma = Float32.(gamma), x0 = Float32.(x0))

        #harmonic bands
        b2h = harmonic_bands(8,0.02)
        b2a = anharmonic_bands(b2h,0.01)
        freq_domains = vcat(b2h,b2a)
        bands_name = savedir*"ta_erai_11a_harmonicbands.jld2"
        ind,protobands,lambda_bands = combined_bands(freq_domains,protofreq,protophases,lambda)
        jldsave(bands_name,ind=ind, protobands=protobands,lambda_bands = lambda_bands)

        #coupling in the harmonic bands
        coupling_name = dir*"ta_erai_11a_harmonicbands_coupling.jld2"
        gamma, x0 = cauchy_by_modes(protobands)
        jldsave(coupling_name,gamma = Float32.(gamma), x0 = Float32.(x0))
        """
        return lambda_bands
    end


    function coupling_fluxnet(outdir,W,vari,method)
        #files from the diff,ssa runs
        #outdir = "/net/scratch/lschulz/ta_erai_11a/"
        #savedir = dir
        freq_domains = simplebands()
        outdir="/net/scratch/lschulz/fluxfullset/"
        for method in ["ssa","diff"],preproc in ["gSSA","lSSA","win.","raw."]
            #GPP
            vari=1
            filename_list = create_file_list(outdir,method,W,vari,preproc)
            lambda,protofreq,RC = extract_from_files(filename_list,5478,48)
            indices = band_indices(freq_domains,protofreq)
            protobands_1,lambda_bands_1 = combine_by_bands(indices,RC,lambda)

            #SOIL TEMPERATURE
            vari=4
            filename_list = create_file_list(outdir,method,W,vari,preproc)
            lambda,protofreq,RC = extract_from_files(filename_list,5478,48)
            indices = band_indices(freq_domains,protofreq)
            protobands_4,lambda_bands_4 = combine_by_bands(indices,RC,lambda)

            for spot in [1,3,5,7,18,20]
                name = "$(method)_$(preproc)_$(spot)_$(W)"

                F = Figure(resolution=(1200,400))
                for (band,bname) in enumerate(["slow","seas/ann","fast"])
                    signal1 = protobands_1[:,spot,band]
                    signal2 = protobands_4[:,spot,band]

                    di = phase_diff_hist(signal2,signal1)
                    ax = Axis(F[1,band],title =bname,xlabel=L"\Delta \phi")
                    ylims!(-1/2,2)
                    s = series!(ax,Array(-pi:0.01:pi)[1:end-1],di[2]')
                end
                save(dir*"test/"*name*".png",F)
            end
        end
        """
        jldsave(bands_name,ind=ind, protobands=protobands,lambda_bands = lambda_bands)

        #coupling in the simple bands
        coupling_name = dir*"ta_erai_11a_simplebands_coupling.jld2"
        gamma, x0 = cauchy_by_modes(protobands)
        jldsave(coupling_name,gamma = Float32.(gamma), x0 = Float32.(x0))

        #harmonic bands
        b2h = harmonic_bands(8,0.02)
        b2a = anharmonic_bands(b2h,0.01)
        freq_domains = vcat(b2h,b2a)
        bands_name = savedir*"ta_erai_11a_harmonicbands.jld2"
        ind,protobands,lambda_bands = combined_bands(freq_domains,protofreq,protophases,lambda)
        jldsave(bands_name,ind=ind, protobands=protobands,lambda_bands = lambda_bands)

        #coupling in the harmonic bands
        coupling_name = dir*"ta_erai_11a_harmonicbands_coupling.jld2"
        gamma, x0 = cauchy_by_modes(protobands)
        jldsave(coupling_name,gamma = Float32.(gamma), x0 = Float32.(x0))
        """
        return lambda_bands
    end

    function bandstrength()
        f = Figure(resolution=(1800,3600))
        plantinds = [4,9,19,10,12,16,7,11,14,18,20,1,3,5]
        ind = ele
        Walist = [5,5.5,6,6.5,7]
        varlist = ["GPP","NEE","RECO","TS"]

        for vari in 1:4,Wa in Walist
            W =  Int.(floor.(Wa .* 365.25))
            ax1 = Axis(f[end+1,1],title="ssa $Wa years $(varlist[vari])")
            series!(ax1,hcat(coupling_analysis(outdir,W,vari,"ssa")[ind,:],sum(coupling_analysis(outdir,W,vari,"ssa"),dims=2)[ind])',labels=["slow","an/sea","fast","sum"])
            ax2 = Axis(f[end,2],title="diff $Wa years $(varlist[vari])")
            s = series!(ax2,hcat(coupling_analysis(outdir,W,vari,"diff")[ind,:],sum(coupling_analysis(outdir,W,vari,"diff"),dims=2)[ind])',labels=["slow","an/sea","fast","sum"])
            axislegend(ax1)
        end
        save(dir*"lambda_bands_elevation.png",f)
    end

    function calculate_bands(datadir,rec_m,var_ind,season_ind)



        f,l = all_mode_protofrequencies(datadir,rec_m,var_ind,1:71,1:48,season_ind)

        b1 = simplebands(0.02)



        i1 = matrix_Tbands(b1,f)
        i2 = matrix_Tbands(b2,f)

        B1 = [Int64.(vcat([vcat([Int64(j) for j in i]...) for i in unique(i2[1:43,k])]...)) for k in 1:71]
        B2 = [Int64.(vcat([vcat([Int64(j) for j in i]...) for i in unique(i2[44:end,k])]...)) for k in 1:71]

        return l,[i1[1,:],i1[2,:],i1[3,:]],[B1,B2]

    end

    """
    different f
    """

    function different_f() #approximate frequency by maximizing correlation to sinus ON ICE
        function protofrequency_by_corr(signal::Vector{Float32},year_sample)
            sinus_f_sampling = 10 .^ (-2:0.05:2)
            sinus_table = hcat([Float32.(sin.(range(0,2*pi,length(signal)) * f *year_sample)) for f in sinus_f_sampling]...)
            maxi = [cor(sinus_table[:,i],signal) for i in 1:length(sinus_f_sampling)]
            return maximum(maxi)
        end

        data = load("/net/scratch/lschulz/LAIFssa_288_3599_1_lSSA.jld2")

        year_sample=900/24
        pf_c = [protofrequency_by_corr(centralizer(data["RC"][:,k]),year_sample) for k=1:48]
        pf_s = [protofrequency(data["RC"][:,k],24) for k=1:48]

        k=2

        s_s = sin.(range(0,2*pi,900) .* pf_s[k] .* year_sample) .+3
        s_c = sin.(range(0,2*pi,900) .* pf_c[k] .* year_sample) .+6
        f,ax,s = series([centralizer(data["RC"][:,k]) s_s s_c]')

        save(dir*"test.png",f)
    end

    # build the variance plot
    #lets do 10 spots
    function varstuff()
        method = "diff"
        Walist = [5,5.5,6,6.5,7]

        vari = 2
        varlist = ["GPP","NEE","RECO","TS"]

        function variance_plot(method,Wa,vari)
            W =  Int.(floor.(Wa .* 365.25))

            bins = 10 .^Array(-2:0.01:2)
            spots = 20
            WW = ones(Float64,length(bins)-1,spots) .*10^-7

            for spot in 1:spots
                data = load("/net/scratch/lschulz/fluxfull/$(method)_$(W)_$(spot)_$(vari)_raw.jld2")
                lambda = data["lambda"]
                pp = Float32.(rec_protophase(data["RC"],1))
                pf = [protofrequency(pp[:,kk]) for kk in 1:48]
                w =fit(Histogram,pf,weights(lambda),bins).weights
                WW[:,spot] = w
            end
            ind=plantinds
            WW = WW[:,ind]
            f,ax,h = heatmap(1:spots,bins[1:end-1],log10.(abs.((WW'))),colormap=:heat,axis = (
                yscale = log10,
                title = "$(method) $Wa years $(varlist[vari])",))
            Colorbar(f[2,1],h,vertical=false,label="log10 var of frequencies")#,ticks = 0:0.1:1)
            ax2 = Axis(f[1,end+1],yscale=log10,title="total spatial variance")
            lines!(ax2,sum(abs.((WW')),dims=1)[:].+10^-7,bins[1:end-1],color="black")
            hlines!(ax2, 365.25/W)

            for band in harmonic_bands(8,0)
                hlines!(ax2,band[1],alpha=0.5,xmax=0.05,color="red")
            end
            linkyaxes!(ax, ax2)


            colsize!(f.layout, 1, Relative(2/3))
            rowsize!(f.layout, 1, Aspect(1, 1))
            save(dir*"var/igbp_var_$(method)_$(W)_$(vari).png",f)
        end

        for  Wa in Walist, method in ["diff","ssa"],vari in 1:4
            variance_plot(method,Wa,vari)
        end
    end

    #build coupling for 2 locations: signal and complete reconstruction


    function individual_coupling_stuff()
        function coupling_plot(Wa,spot1,spot2,vari1,vari2)
            W = Int.(floor.(Wa .* 365.25))

            f=Figure(resolution=(1800,900))

            for (i,method) in enumerate(methods)
                data1 = load("/net/scratch/lschulz/fluxfull2/$(method)_$(W)_$(spot1)_$(vari1)_raw.jld2")
                lambda1 = data1["lambda"]
                signal1 = data1["signal"]
                rc1 = sum(data1["RC"],dims=2)[:]
                pp1 = rec_protophase(sum(data1["RC"],dims=2),1)[:]

                data2 = load("/net/scratch/lschulz/fluxfull2/$(method)_$(W)_$(spot2)_$(vari2)_raw.jld2")
                lambda2 = data2["lambda"]
                signal2 = data2["signal"]
                rc2 = sum(data2["RC"],dims=2)[:]
                pp2 = rec_protophase(sum(data2["RC"],dims=2),1)[:]

                #coupling signal
                coupling_signal = phase_diff_hist(signal1,signal2)

                #coupling reconstruction
                coupling_rc = phase_diff_hist(rc1,rc2)

                #coupling protophases
                coupling_pp = phase_diff_hist(pp1,pp2)

                #plotting
                
                ax1 = Axis(f[i,1],title = method*" decomposition $(sum(lambda1)) $(sum(lambda2))",xlabel = L"t")
                s = series!(ax1,hcat(signal1,signal2.+5,rc1.+10,rc2.+15,pp1.+20,pp2.+25)',labels=["s1","s2","rc1","rc2","pp1","pp2"])
                axislegend(ax1)

                ax2 = Axis(f[i,2],title = method*" coupling between signals",xlabel=L"\Delta \phi")
                ylims!(-1/2,2)
                s = series!(ax2,Array(-pi:0.01:pi)[1:end-1],hcat(coupling_signal[2],coupling_rc[2],coupling_pp[2])',labels=["signal","rc","pp"])
                axislegend(ax2)

            end
            Label(f[0, :], text = "$Wa years $spot1 $(varlist[vari1]) $spot2 $(varlist[vari2])", textsize = 20)
            name = "$(Wa)_$(spot1)_$(vari1)_$(spot2)_$(vari2)"
            save(dir*"cop/"*name*".png",f)
        end

        for Wa in [5,5.5,6,6.5,7], vari1 in 1:4, vari2 in 1:4
            spot1 = 4
            spot2 = 19
            coupling_plot(Wa,spot1,spot2,vari1,vari2)
            spot1 = 4
            spot2 = 4
            coupling_plot(Wa,spot1,spot2,vari1,vari2)
            spot1 = 19
            spot2 = 19
            coupling_plot(Wa,spot1,spot2,vari1,vari2)
        end
    end

    """
    hierarchy based robustness
    """

    #this just sorts the RC by the variance and gives the first kappa #MOVE TO SINGLE SPOT SINGLE TIMESERIES ANALYSIS
    function robustness(method,W,vari,preproc,kappa,spots)
        N = 5478
        k = 48

        outdir="/net/scratch/lschulz/fluxfullset/"
        filename_list = create_file_list(outdir,method,W,vari,preproc)
        lambda,protofreq,RC = extract_from_files(filename_list,N,k)
        lambda_l = hcat([(sort(lambda[:,spot],rev=true)[1:kappa]) for spot in 1:spots]...) #kappa x spots
        RC_l = reshape(hcat([RC[:,sortperm(lambda[:,spot],rev=true)[1:kappa],spot] for spot in 1:spots]...),N,kappa,spots)
        return RC_l,lambda_l
    end

        #this just sorts the RC by the variance and gives the first kappa
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


    function strongmodes_by_W(outdir,method,Ws,vari,preproc,spot,kappa,bins,year)
        WW = ones(Float32,length(bins)-1,length(Ws)) .*10^-7

        for (i,W) in enumerate(Ws)
            RC_l, lambda_l,pf_l = single_robustness(outdir,method,W,vari,preproc,kappa,spot,year)
            w =fit(Histogram,pf_l,weights(lambda_l),bins).weights
            WW[:,i] = Float32.(w)
        end
        return WW
    end

end


function readout_and_combine()
    spot = 2
    vari = 1
    kappa = 48
    year_samples=24
    bins = 10 .^Array(-2:0.05:2)
    Walist = [5,5+1/4,5+1/3,5+1/2,5+2/3,5+3/4,6,6+1/4,6+1/3,6+1/2,6+2/3,6+3/4,7]
    #Ws = Int.(floor.(Walist .* 365.25))
    Ws = 24 .*Array(11:17)
    outdir = "/net/scratch/lschulz/"
    f = Figure(resolution=(1800,600))
    for (i,preproc) in enumerate(["raw.","gSSA","lSSA"])#,"win."])
        for (j,method) in enumerate(["LAIFssa","LAIFdiff"])
            modi_kappa = strongmodes_by_W(outdir,method,Ws,vari,preproc,spot,kappa,bins,year_samples)
            ax = Axis(f[j,i],title = method*" "*preproc*" ",xlabel="W",yscale = log10,ylabel="f")
            h = heatmap!(ax,Ws,bins[1:end-1],log10.(abs.((modi_kappa'))),colormap=:heat,colorlimits=(-4,1))
        #Colorbar(f[:,end+1],h,vertical=false,label="log10 var of frequencies")#,ticks = 0:0.1:1)
        end
    end

    save(dir*"test.png",f)
end

LL = Matrix{Float64}(undef,length(bins)-1,700)
for spot = 1:700
    LL[:,spot] = sum(strongmodes_by_W(outdir,method,Ws,vari,preproc,spot,kappa,bins,year_samples),dims=2)[:]
end

function put_to_map()
    #img = FileIO.load("logs/KW_29/"*"worldmap.png")
    F = Figure()
    a = Axis(F[1,1],limits=(-5, 20, 40, 60))
    #image!(a,[-180,180],[-90,90],rotr90(img),yflip=true)

    """
    bbox-east-long 	180.0
    bbox-north-lat 	90.0
    bbox-south-lat 	-90.0
    bbox-west-long 	-180.0
    """

    #4320 longitude (degrees east) x 2160 latitude degrees north
    lon = Array(range(-180,180,4320))
    lat = Array(range(90,-90,2160))
    XX = Matrix{Union{Float32,Missing}}(undef,4320,2160)
    


    using Zarr

    zz  = zopen("/net/data/LAI/LAI_AVHRR.zarr",fill_as_missing=true)
    x = zz.arrays["layer"]

    longitude_w = 2101:2290
    latidute_w = 401:590

    data = @view x[longitude_w,latidute_w,:]
    datacols = reshape(data[:],size(data,1)*size(data,2),912)

    function sectionlength(coli::SubArray{Union{Missing,Float32}})
        L = 0
        l = 0
        for i in coli
            l +=1
            L = (l>L) ? l : L
            if ismissing(i)
                l=0
            end
        end
        return L
    end

    l = [sectionlength(@view datacols[i,13:end]) for i=1:size(datacols,1)]

    filler = Matrix{Union{Missing,Float32}}(undef,length(longitude_w),length(latidute_w))
    filler[findall(x-> x== 900,reshape(l,length(longitude_w),length(latidute_w)))[1:5:end]] .= 1.0
    XX[longitude_w,latidute_w] = filler
    heatmap!(lon,lat,XX)
    save(dir*"map.png",F)

    sery = Float32.(datacols[findall(x-> x== 900,l)[1:5:end],13:end])

    jldsave("/net/scratch/lschulz/LAIF/data.jld2",data=sery,filler=filler)

end

#make some large distribution plot for each preproc/method
#both in linear and in log10 scale
#show the harmonics
#showcase some interesting frequency bands
bins = 10 .^Array(-2:0.01:2)
outdir="/net/scratch/lschulz/LAIF/"
Ws = [2556]
preproc = ["raw.","lSSA"][1]
method = ["diff","ssa"][2]
vari = 4
kappa = 48
year_samples= 365.25

"""
LAIF
"""
Ws = Int.(floor.([15,16,17]*24))
LL = Array{Union{Missing,Float32}}(undef,length(bins)-1,3,4080)
for spot = 1:4080
    LL[:,:,spot] = strongmodes_by_W(outdir,method,Ws,vari,preproc,spot,kappa,bins,year_samples)
end

for spot in 1:4080, W in Ws, preproc = ["raw","lSSA"],method = ["diff","ssa"]
    if !isfile(outdir*method*"_$(W)_$(spot)_1_$(preproc).jld2")
        println(outdir*method*"_$(W)_$(spot)_1_$(preproc).jld2")
    end
end

"""
FLUX
"""
LL = Array{Union{Missing,Float32}}(undef,2,2,length(bins)-1,20)
for spot = 1:20, (i,preproc) = enumerate(["raw.","lSSA"]), (j,method) = enumerate(["diff","ssa"])
    LL[i,j,:,spot] = strongmodes_by_W(outdir,method,Ws,vari,preproc,spot,kappa,bins,year_samples)
end
"""
FLUX only raw BUT different Ws
"""
Ws = Int.(floor.([5,5+1/4,5+1/3,5+1/2,5+2/3,5+3/4,6,6+1/4,6+1/3,6+1/2,6+2/3,6+3/4,7].*365.25))
LL = Array{Union{Missing,Float32}}(undef,2,length(Ws),length(bins)-1,20)
for (i,spot) = enumerate(1:20), (m,method) = enumerate(["ssa","diff"]), (w,W) in enumerate(Ws)
    LL[m,w,:,i] = strongmodes_by_W(outdir,method,W,vari,preproc,spot,kappa,bins,yearsamples)
end
"""
LAI only raw BUT different Ws
"""
Ws = Int.(floor.([5,5+1/4,5+1/3,5+1/2,5+2/3,5+3/4,6,6+1/4,6+1/3,6+1/2,6+2/3,6+3/4,7].*365.25))
LL = Array{Union{Missing,Float32}}(undef,2,3,length(bins)-1,20)
spots = rand(1:4000,20)
for (i,spot) = enumerate(spots), (m,method) = enumerate(["ssa","diff"]), (w,W) in enumerate(Ws)
    LL[m,w,:,i] = strongmodes_by_W(outdir,method,W,vari,preproc,spot,kappa,bins,year_samples)
end
#========================================================================================#

"""
 complete variance LAIF
"""
begin
    lon = Array(range(-180,180,4320))
    lat = Array(range(90,-90,2160))
    longitude_w = 2101:2290
    latidute_w = 401:590
    f = load("/net/scratch/lschulz/LAIF/data.jld2")["filler"]
    indices = findall(!ismissing,f)



    variance_full = sum(LL,dims=1)[:]



    XX = Matrix{Union{Float32,Missing}}(undef,4320,2160)
    filler = Matrix{Union{Missing,Float32}}(undef,length(longitude_w),length(latidute_w))
    #maybe even add another heatmap? ON ICE

    filler[indices] .= variance_full ./3#rand(length(indices))
    XX[longitude_w,latidute_w] = filler

    F = Figure(resolution=(900,1200))
    a = Axis(F[1,1],limits=(-7, 12, 40, 58))
    heatmap!(lon,lat,ones(length(lon),length(lat)),color="black")
    h = heatmap!(lon,lat,XX,colormap = cgrad(:lajolla, 10, categorical = true),colorrange=(0,1))

    Colorbar(F[end+1,:],h,vertical = false,label="test")
    save(dir*"map.png",F)

    #histogram
    F = Figure(resolution=(1800,900))
    for (j,method) = enumerate(["diff","ssa"]),(k,Ws) =enumerate([15,16,17]), (i,preproc) = enumerate(["raw","lSSA"])
        name = "laif $(method) $(Ws) $(preproc)"
        LL = load(dir*"$(i)_$(j).jld2")["data"]
        histogram_full = sum(LL,dims=3)[:,k,1]
        ax = Axis(F[k,(i!=j) ? j+1 : i*j],title=name,xscale=log10,yscale=log10,xlabel=L"f",ylabel=L"\sum \lambda")
        series!(ax,bins[1:end-1],(histogram_full.+10^-5)')
        vlines!(ax,hcat(harmonic_bands(7,0.0)...)[1,:],ymin=0.0,ymax=0.05)
        vlines!(ax, 1/Ws)
    end
    save(dir*"hist.png",F)
end
#=============================================================================#


#band plot
b1_ind = 21:24
ranges=[(0,0.0000001),(0.00001)]

F = Figure(resolution=(1800,1200))
hs = []
for (j,method) = enumerate(["diff","ssa"]),(k,Ws) =enumerate([15,16,17]), (i,preproc) = enumerate(["raw","lSSA"])
    name = "laif $(method) $(Ws) $(preproc)"
    LL = load(dir*"$(i)_$(j).jld2")["data"]
    ax = Axis(F[k,((i!=j) ? j+1 : i*j)],title=name,limits=(-7, 12, 40, 58))

    XX = Matrix{Union{Float32,Missing}}(undef,4320,2160)
    filler = Matrix{Union{Missing,Float32}}(undef,length(longitude_w),length(latidute_w))
    
    
    filler[indices] .= sum(LL[b1_ind,k,:],dims=1)[1,:,][:]
    XX[longitude_w,latidute_w] = filler

    heatmap!(ax,lon,lat,ones(length(lon),length(lat)),color="black")
    h = heatmap!(ax,lon,lat,XX,colormap = cgrad(:lajolla, 20, scale=log10,categorical = true),colorrange=ranges[j])
    hs=push!(hs,h)
end

save(dir*"map_seven.png",F)

#=============================================================================#


"""
BOTH
variance stuff: plot both lambda over f
"""

outdir="/net/scratch/lschulz/fluxfullset/"
meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
W = 2556
preproc = "raw."
kappa= 48
yearsamples=365.25
vari = 1
currentclass = "test"


for spot in sortperm([meta[spot,"IGBP_class"] for spot in 1:20])

    F = Figure(resolution=(1200,500))
    class = meta[spot,"IGBP_class"]
    tit = "$(spot)_$(class)_raw_GPP"

    ax = Axis(F[1,1],yscale=log10,xscale=log10,ylabel=L"\lambda",xlabel=L"f",limits=(10^-1.5,10^1.5,10^-5,10^0.5),title=tit)

    for method in ["ssa","diff"]
        r_ssa = single_robustness(outdir,method,W,vari,preproc,kappa,spot,yearsamples)
        scatter!(ax,r_ssa[3],r_ssa[2],marker=:+,markersize=30,markerstrokewidth=10,label=method)
    end
    axislegend(ax)
    vlines!(ax,[harmonic_bands(4,0)[i][1] for i in 1:11],color="grey")
    text!(ax,string.(round.([harmonic_bands(4,0)[i][1] for i in 1:11],digits=2)),
    position=[(harmonic_bands(4,0)[i][1],1.0) for i in 1:11],textsize=10)
    save(dir*"both_GPP/$tit.png",F)
end


#=============================================================================#

"""
individual ecosystem combined pictures
"""

F = Figure(resolution=(4*1200,1*700))
for (i,spot) in enumerate(sortperm([meta[spot,"IGBP_class"] for spot in 1:20])[20])
    method = "ssa"
    r_ssa = single_robustness(outdir,method,W,vari,preproc,kappa,spot,yearsamples)
    class = meta[spot,"IGBP_class"]
    tit = "$(spot)_$(class)_$(method)_raw_GPP"


    ax1 = Axis(F[i,1],yscale=log10,title = tit,limits=(-1,50,10^-4,10^0.3),xlabel=L"k",ylabel=L"\lambda")
    ax2 = Axis(F[i,2],title=tit,xlabel="T",ylabel=L"GPP")

    scatter!(ax1,r_ssa[2],marker=:x,markersize=20,)
    series!(ax2,hcat([r_ssa[1][:,1:10][:,k] .+ 2*k for k=1:10]...)',color=:viridis,)
    circlepoints=Array(2:2:20)
    text!(ax2,
    circlepoints,
    text = string.(round.(r_ssa[3][1:10],digits=2)),
    rotation = 0,#LinRange(0, 2pi, 16)[1:end-1],
    #align = (:right, :baseline),
    #color = cgrad(:Spectral)[LinRange(0, 1, 15)]
    )

    ax3 = Axis(F[i,3],xscale=log10,yscale=log10,xlabel=L"f",ylabel=L"\lambda",limits=(10^-1,10^1,10^-5,10^0.5),)
    scatter!(ax3,r_ssa[3],r_ssa[2],marker=:+,markersize=20,color=:black)
    vlines!(ax3,[harmonic_bands(8,0)[i][1] for i in 1:43],color="red")


    method = "diff"
    r_ssa = single_robustness(outdir,method,W,vari,preproc,kappa,spot,yearsamples)
    class = meta[spot,"IGBP_class"]
    tit = "$(spot)_$(class)_$(method)_raw_GPP"
    

    ax1 = Axis(F[i,4],yscale=log10,title = tit,limits=(-1,50,10^-4,10^0.3),xlabel=L"k",ylabel=L"\lambda")
    ax2 = Axis(F[i,5],title=tit,xlabel="T",ylabel=L"GPP")
    
    scatter!(ax1,r_ssa[2],marker=:x,markersize=20,)
    series!(ax2,hcat([r_ssa[1][:,1:10][:,k] .+ 2*k for k=1:10]...)',color=:viridis,)
    circlepoints=Array(2:2:20)
    text!(ax2,
    circlepoints,
    text = string.(round.(r_ssa[3][1:10],digits=2)),
    rotation = 0,#LinRange(0, 2pi, 16)[1:end-1],
    #align = (:right, :baseline),
    #color = cgrad(:Spectral)[LinRange(0, 1, 15)]
    )

    ax3 = Axis(F[i,6],xscale=log10,yscale=log10,xlabel=L"f",ylabel=L"\lambda",limits=(10^-1,10^1,10^-5,10^0.5),)
    scatter!(ax3,r_ssa[3],r_ssa[2],marker=:+,markersize=28,color=:black)
    vlines!(ax3,[harmonic_bands(8,0)[i][1] for i in 1:43],color="red")
    
end

save(dir*"opensavannah.png",F)

#=============================================================================#

"""
bins (saved in a file as above
"""

#size LL 2x2x80x20
#preproc method bins spots


for spot = 1:20
    F = Figure(resolution=(1800,1200))
    for (i,preproc) = enumerate(["raw.","lSSA"]), (j,method) = enumerate(["diff","ssa"])
        class = meta[spot,"IGBP_class"]
        tit = "$(spot)_$(class)_$(method)_$(preproc)_GPP"
        hist = LL[i,j,:,spot]
        ax = Axis(F[i,j],yscale=log10,xscale=log10,ylabel=L"\lambda",xlabel=L"f",limits=(10^-2,10^2,10^-5.1,10^0.5),title=tit)
        series!(ax,bins[1:end-1],hist' .+ 10^-5)
        vlines!(ax,[harmonic_bands(4,0)[i][1] for i in 1:11],color="red")
    end
    save(dir*"bands/$spot.png",F)
end



"""
does not seem to amount to much
"""

#=============================================================================#

"""
band-coupling GPP TS
"""

for spot = 1:20
    F = Figure(resolution=(1800,1200))
    number=1
    for (i,preproc) = enumerate(["raw.","lSSA"]), (j,method) = enumerate(["diff","ssa"])
        class = meta[spot,"IGBP_class"]
        tit = "$(spot)_$(class)_$(method)_$(preproc)_GPP"

        #simple bands
        freq_domains = simplebands(0.0)


        #first read in GPP
        r_GPP = single_robustness(outdir,method,W,1,preproc,kappa,spot,yearsamples) #RC,lambda,freq
        protofrequencies = reshape(r_GPP[3],48,1)
        RC = reshape(r_GPP[1],5478,48,1)
        lambdas = reshape(r_GPP[2],48,1)
        indices = band_indices(freq_domains,protofrequencies) #return ind (bands, L) vec
        protobands_GPP, lambda_bands_GPP = combine_by_bands(indices,RC,lambdas) #return protobands(N,L,bands),lambda_bands(N,L,bands)

        #then read in TS
        r_TS = single_robustness(outdir,method,W,4,preproc,kappa,spot,yearsamples) #RC,lambda,freq
        protofrequencies = reshape(r_GPP[3],48,1)
        RC = reshape(r_TS[1],5478,48,1)
        lambdas = reshape(r_TS[2],48,1)
        indices = band_indices(freq_domains,protofrequencies) #return ind (bands, L) vec
        protobands_TS, lambda_bands_TS = combine_by_bands(indices,RC,lambdas) #return protobands(N,L,bands),lambda_bands(N,L,bands)

        #phase difference standard ssa diff
        #phase difference bands ssa diff
        for (k,band) in enumerate(["slow","annual seasonal","fast"])
            ax = Axis(F[number,k],limits=(-pi,pi,-1,2),xlabel=L"\delta \phi",ylabel=L"p",title=tit*" $band")
            bins,hist = phase_diff_hist(protobands_GPP[:,1,k],protobands_TS[:,1,k]) #gives bins,cauchy_weights
            series!(ax,bins[1:end],hist')
        end
        number+=1

    end
    save(dir*"bands_coupling/$spot.png",F)
end

#=============================================================================#

"""
the reconstruction vs signal quality plots
"""

outdir="/net/scratch/lschulz/fluxfullset/"
meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
W = 2556
preproc = "raw."
kappa= 48
yearsamples=365.25
vari = 1
currentclass = "test"
L = 5478
startyear = 2005
years = startyear .+ ((1:L) ./ 365.25)
limits = (startyear,startyear+15,-1.5,4.0)
function plot_rec(spot,savedirname)
    F = Figure(resolution=(600,400))
    axis_list = [Axis(F[1,1],ylabel = "GPP",limits=limits,xticks = LinearTicks(15),yticks=[0.0],xticklabelrotation=45),
                Axis(F[2,1],ylabel="GPP by SSA",limits=limits,xticks = LinearTicks(15),yticks=[0.0],xticklabelrotation=45),
                Axis(F[3,1],ylabel="GPP by NLSA",xlabel="t [a]",limits=limits,xticks = LinearTicks(15),yticks=[0.0],xticklabelrotation=45)]

    for (i,method) = enumerate(["ssa","diff"])

        Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
        Savefile = load(Filename)
        signal = Savefile["signal"]
        r_GPP = single_robustness(outdir,method,W,1,preproc,kappa,spot,yearsamples) #RC,lambda,freq
        RC = r_GPP[1]
        lambda = r_GPP[2]
        protofreq = r_GPP[3]


        lines!(axis_list[i+1],years,sum(RC,dims=2)[:],label=method,color="black")
        if method=="ssa" lines!(axis_list[1],years,signal,label="signal",color="black") end

        rc = sum(RC,dims=2)[:]
        #abs_norm = string.(round.(norm(rc-signal),digits=2))
        #cos_an = string.(round.(cos((dot(rc-signal,rc-signal)^2 / norm(rc) / norm(signal))),digits=2))
        #allvar = string.(round.(1-var(rc),digits=2))
        #mse = string.(round.(dot(rc-signal,rc-signal),digits=2) / L )
        #txt = "rec norm " * abs_norm * " rec angle " * cos_an * " rec var " * allvar * " mse " * mse
        #text!(axis_list[i+1],txt,position=(years[1],3),textsize=16)
        
    end

    linkxaxes!(axis_list[1],axis_list[2])
    linkxaxes!(axis_list[2],axis_list[3])
    hidespines!(axis_list[1], :t, :r, :b, :l)
    hidespines!(axis_list[2], :t, :r, :b, :l)
    hidespines!(axis_list[3], :t, :r, :b, :l)

    hidexdecorations!(axis_list[1],grid=false)
    hidexdecorations!(axis_list[2],grid=false)
    hidexdecorations!(axis_list[3],ticks=false,label=false,ticklabels=false,grid=false)


    save(savedirname*".png",F)
end

for spot in [3,6,19]
    plot_rec(spot,dir*"pic/rec_$spot")
end
#=============================================================================#

"""
the big component investigation plot
"""

spot = 19
outdir="/net/scratch/lschulz/fluxfullset/"
meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
W = 2556
preproc = "raw."
kappa= 48
yearsamples=365.25
vari = 4
L = 5478
bins = 10 .^Array(-2:0.01:2)


function big_figure(spot,savedirname)
    F = Figure(resolution= (800,700))
    ax_hist = Axis(F[1,1],xscale=log10,ylabel=L"\sum \lambda",xlabel="f [a⁻¹]",limits=(10^-1.5,10^1.5,-.1,0.8))
    ax_text = Axis(F[1,2])
    ax_both = Axis(F[2,1],yscale=log10,xscale=log10,ylabel=L"\lambda",xlabel="f [a⁻¹]",limits=(10^-1.5,10^1.5,10^-5,10^0.5))
    ax_spec = Axis(F[2,2],xlabel="k",ylabel=L"\lambda",yscale=log10,limits=(0,49,10^-4,10^0.5))


    for (i,method) = enumerate(["ssa","diff"])


        r_GPP = single_robustness(outdir,method,W,vari,preproc,kappa,spot,yearsamples) #RC,lambda,freq
        RC = r_GPP[1]
        lambda = r_GPP[2]
        protofreq = r_GPP[3]
        scatter!(ax_both,protofreq,lambda,marker=[:cross,:xcross][i],markersize=30,markerstrokewidth=10,label=method)

        h_weights =fit(Histogram,protofreq,weights(lambda),bins).weights
        lines!(ax_hist,bins[1:end-1],10^-5 .+ h_weights)

        lines!(ax_spec,1:48,lambda)
    end

    hidespines!(ax_both, :t, :r)
    axislegend(ax_both)

    harm = list_pure_harmonics(7)

    vlines!(ax_both,harm,color="grey")
    vlines!(ax_hist,harm,color="grey")
    text!(ax_both,string.(round.(harm,digits=2)),
    position=[(i,10^-4) for i in harm],textsize=10)

    hidespines!(ax_hist, :t, :r, :b)
    hidexdecorations!(ax_hist,ticks=false,grid=false)


    #linkxaxes!(ax_hist,ax_both)
    #linkyaxes!(ax_both,ax_spec)

    colsize!(F.layout, 1, Fixed(400))
    colsize!(F.layout, 2, Fixed(200))

    rowsize!(F.layout, 1, Fixed(200))
    rowsize!(F.layout, 2, Fixed(300))
    save(savedirname*".png",F)
end

for spot in 1:20
    big_figure(spot,dir*"bigf/$spot")
end
#=============================================================================#

"""
individual component annual investigation plot
first find out the components in question  || manual replaced by automatic
then plot them to see whats going on
"""


spot = 14
outdir="/net/scratch/lschulz/fluxfullset/"
meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
W = 2556
preproc = "gSSA"
kappa= 48
yearsamples=365.25
vari = 1
N = 5478
k=48
bins = 10 .^Array(-2:0.02:2)

function component_annual_investigation(spot,savedirname)
    F = Figure()
    axi1 = [Axis(F[1,1],title = "SSA RC"),Axis(F[1,2],title = "NLSA RC")]
    axi2 = [Axis(F[2,1],title = "SSA EOF"),Axis(F[2,2],title = "NLSA EOF")]
    axi3 = [Axis(F[3,1],title = "SSA PC"),Axis(F[3,2],title = "NLSA PC")]



    for (m,method) = enumerate(["ssa","diff"])

        N = 5478
        k = 48
        Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
        lambda,protofreq,RC = extract_from_single_file(Filename,yearsamples,N,k)
        protofreq = protofreq[:]
        
        indices = findall(x->abs(x-1.0)<= 0.05, protofreq)


        for (i,index) in enumerate(indices)
            lines!(axi1[m],RC[:,index] .+ i*2,color="black")
            text!(axi1[m],string.(round.(protofreq[index],digits=2)),position=(0,i*2))
        end

        lines!(axi1[m],sum(RC[:,indices],dims=2)[:] .- 5)


        Savefile = load(Filename)
        EOF = Savefile["EOF"]


        for (i,index) in enumerate(indices)
            lines!(axi2[m],EOF[:,index].+ i* 0.1,color="black")
        end

        lines!(axi2[m],sum(EOF[:,indices] .- 0.2,dims=2)[:])

        PC = Savefile["PC"]


        for (i,index) in enumerate(indices)
            lines!(axi3[m],PC[:,index].+ i* 60,color="black")
        end

        lines!(axi3[m],sum(PC[:,indices] .- 40,dims=2)[:])

    end

    for a in [axi1...,axi2...,axi3...]
        hidexdecorations!(a)
        hideydecorations!(a)
        hidespines!(a)
    end

    save(savedirname*".png",F)
end

for spot in 1:20
    component_annual_investigation(spot,dir*"annual_pp/$spot")
end

#=============================================================================#

"""
TABLEMAKING

IGBP
Height
locNr
annual strength SSA Diff
harmonic power SSA Diff
Annual Seasonal Power SSA Diff
Diff eps

this absolute error makes it look like everything is in the harmonics
"""

meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]

spot = 14

function bands_investigation(spot)

    outdir="/net/scratch/lschulz/fluxfullset/"
    W = 2556
    preproc = "raw."
    kappa= 48
    yearsamples=365.25
    vari = 1
    L = 5478
    bins = 10 .^Array(-2:0.02:2)
    N = 5478
    k = 48
    eps = 0.001


    harmonics = list_pure_harmonics(7)
    function harm_indi(protofreq,eps) L = []; [append!(L,findall(x->abs(x-i)<= eps, protofreq)) for i in harmonics][:]; return unique(L) end
    band_indi(protofreq,eps,low,up) = findall(x->low<=x<=up, protofreq)


    slow_bounds = [10^-2,0.2]
    seasonal_bounds = [0.2,1.2]
    fast_bounds = [1.2,10^2]


    function bands_invest(method)


        Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
        lambda,protofreq,RC = extract_from_single_file(Filename,yearsamples,N,k)
        protofreq = protofreq[:]
        #total coverage
        total_var = sum(lambda)
        #annual power
        indices_a = findall(x->abs(x-1.0)<= 0.02, protofreq)
        annual = sum(lambda[indices_a])
        #harmonic power
        indices_h = harm_indi(protofreq,eps)
        harmonics = sum(lambda[indices_h])
        #band power
        slow = sum(lambda[band_indi(protofreq,eps,slow_bounds...)])
        seasonal = sum(lambda[band_indi(protofreq,eps,seasonal_bounds...)])
        fast = sum(lambda[band_indi(protofreq,eps,fast_bounds...)])
        return total_var,annual,harmonics,slow,seasonal,fast
    end

    #IGBP
    class = meta[spot,"IGBP_class"]
    #Height
    Height = meta[spot,"elevation"]

    eps = load(create_file_list(outdir,"ssa",W,vari,preproc)[spot],"eps")

    return class,Height,eps,round.(bands_invest("ssa"),digits=4)...,round.(bands_invest("diff"),digits=4)...
end

function bands_print()
    attr = ["class","height","kernel","var s", "an s","harmo s","slow s","seas s","fast s",
                "var d", "an d","harmo d","slow d","seas d","fast d"]
    for i in attr
        print(i*"\t")
    end
    print("\n")
    for spot = sortperm([meta[spot,"IGBP_class"] for spot in 1:20])
        for i in bands_investigation(spot)
            print("$i\t")
        end
        print("\n")
    end
end

#=============================================================================#

"""
relaitve harmonics : not anymore
"""
eps_rel = 0.01

function table_content(method,spot,W)

    function freq_power(protofreq,harmonic,eps)

        indices_a = findall(x->harmonic*(1-eps)<=x<=harmonic*(1+eps), protofreq)
        return sum(lambda[indices_a])
    end

    function freq_power_uplow(protofreq,up,low)
        indices_a = findall(x->low<x<=up, protofreq)
        return sum(lambda[indices_a])
    end

    Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
    lambda,protofreq,RC = extract_from_single_file(Filename,yearsamples,N,k)
    signal = load(Filename)["signal"]
    protofreq = protofreq[:]

    #band strength harmonics

    harmonics = unique(round.(vcat([i for i=1:7],[1/i for i=1:7]),digits=6))
    overtones = unique(round.([i!=j ? (k*i)/j : 0 for i=1:8 for j=1:8 for k=1:3],digits=6))
    harmonic_power = sum([freq_power(protofreq,harmonic,eps_rel) for harmonic = harmonics])
    overtone_power = sum([freq_power(protofreq,harmonic,eps_rel) for harmonic = overtones]) - harmonic_power

    #band strength fast seas slow

    fast = freq_power_uplow(protofreq,10^2,2)
    seas = freq_power_uplow(protofreq,2,1/2)
    slow = freq_power_uplow(protofreq,1/2,10^-2)

    #complete lambda

    signal = load(Filename)["signal"]
    rc = sum(RC,dims=2)[:]
    e = rc .- signal

    rc_l = var(rc)
    e_l = var(e)
    mse = dot(e,e) / N
    bias = mean(e).^2

    return [rc_l,e_l,mse,bias,fast,seas,slow,harmonic_power,overtone_power]
end # [rc_l,e_l,mse,e_bias,fast,seas,slow,harmonic_power,overtone_power]

table_results = Array{Float32}(undef,20,2,9)

for spot=1:20,(i,method) = enumerate(["ssa","diff"])
    table_results[spot,i,:] = table_content(method,spot)
end

"""
table plot
"""
begin

    for spot in sortperm([meta[spot,"IGBP_class"] for spot in 1:20]), (j,method) in enumerate(["ssa","diff"])
        class = meta[spot,"IGBP_class"]
        print("$class\t$spot\t$method\t")
        [print("$(round.(i,digits=3))\t") for i in table_results[spot,j,:]]
        println()
    end

    f,ax,s1 = scatter(table_results[:,1,8],table_results[:,1,1],label="ssa rc lambda",axis=(xlabel="harmonic power",))
    s2 = scatter!(ax,table_results[:,2,8],table_results[:,2,1],label="diff rc lambda")

    s3 = scatter!(ax,table_results[:,1,8],table_results[:,2,3],label="ssa mse")
    s4 = scatter!(ax,table_results[:,2,8],table_results[:,2,3],label="diff mse")

    s5 = scatter!(ax,table_results[:,1,8],table_results[:,2,4],label="ssa bias")
    s6 = scatter!(ax,table_results[:,2,8],table_results[:,2,4],label="diff bias")

    #s7 = scatter!(ax,table_results[:,1,8],table_results[:,2,5],label="ssa slow")
    #s8 = scatter!(ax,table_results[:,2,8],table_results[:,2,5],label="diff slow")

    #s9 = scatter!(ax,table_results[:,1,8],table_results[:,2,9],label="ssa overtones")
    #s10 = scatter!(ax,table_results[:,2,8],table_results[:,2,9],label="diff overtones")

    Legend(f[1,2],ax)

    indi=[8,15,2,6,17,4,9,19,13,10,12,20,18,14,11,16,5,3,1]
    F = Figure(resolution=(1200,900))
    a = Axis(F[1,1],title="ssa",xticks=(1:19,meta[:,"IGBP_class"][indi]),xticklabelrotation=45.0)
    series!(a,table_results[indi,1,:]',color=:Set1,labels=["lam rc" "lam e" "mse" "bias" "fast" "seas" "slow" "harm" "overt"][:])
    b = Axis(F[2,1],title="diff",xticks=(1:19,meta[:,"IGBP_class"][indi]),xticklabelrotation=45.0)
    series!(b,table_results[indi,2,:]',color=:Set1,labels=["lam rc","lam e","mse","bias","fast","seas","slow","harm","overt"])
    Legend(F[1,2],a)
    Legend(F[2,2],b)
    save(dir*"test.png",F)
end
#=============================================================================#

"""
components with the strongest variance
choose a number kappa
"""
method="ssa"
spot = 19
outdir="/net/scratch/lschulz/fluxfullset/"
meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
W = 2556
preproc = "raw."
kappa= 16
yearsamples=365.25
vari = 1
L = 5478
bins = 10 .^Array(-2:0.02:2)

e= Array(1:W) ./ yearsamples
p = Array(1:(L-W+1)) ./ yearsamples
xpos = -1.2

function component_variance_investigation(spot,savedirname)
    F = Figure(resolution=(840,700))

    axi1_EOF = [Axis(F[1,1],title = "EOF",titlealign = :center,xlabel="t [a]",xticks=Array(1:7)),Axis(F[2,1], title=" ",titlealign = :center,xlabel="t [a]",xticks=Array(1:7))]
    axi2_PC = [Axis(F[1,2],title="PC",xlabel="t [a]",xticks=Array(1:7)),Axis(F[2,2],xlabel="t [a]",title=" ",xticks=Array(1:7))]
    #axi3_RC = [Axis(F[3,1],title = "rc ssa", limits=(0,5478,-2,2),titlealign = :left,subtitlefont = "TeX Gyre Heros Italic Makie",),Axis(F[3,2], limits=(0,5478,-2,2),title = "rc nlsa", titlealign = :left,subtitlefont = "TeX Gyre Heros Italic Makie",)]
    for (m,method) = enumerate(["ssa","diff"])

        N = 5478
        k = 48
        Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
        lambda,protofreq,RC = extract_from_single_file(Filename,yearsamples,N,k)
        protofreq = protofreq[:]
        
        indices = sortperm(lambda,rev=true)[1:kappa]

        Savefile = load(Filename)
        EOF = Savefile["EOF"]
        PC = Savefile["PC"]

        for (i,index) in enumerate(indices)
            lines!(axi1_EOF[m],e,EOF[:,index] .- i*0.2,color="black")
            lines!(axi2_PC[m],p,PC[:,index].- i* 100,color="black")
            #lines!(axi3_RC[m],sum(RC[:,index],dims=2)[:],color="blue")

            p_fr = string.(round.(protofreq[index],digits=2))
            l = string.(round.(lambda[index],digits=3))

            text!(axi1_EOF[m],p_fr,position=(xpos,-i*0.2 -0.1),textsize=12)
            text!(axi2_PC[m],l,position=(xpos,-i*100-20),textsize=12)

        end

        lines!(axi1_EOF[m],e,sum(EOF[:,indices],dims=2)[:])
        lines!(axi2_PC[m],p,sum(PC[:,indices],dims=2)[:])
        #lines!(axi3_RC[m],sum(RC[:,indices],dims=2)[:])

        #text!(axi1_EOF[m],"EOF",position=(xpos,0.2),textsize=16)
        #text!(axi2_PC[m],"PC",position=(xpos,140),textsize=16)

        text!(axi1_EOF[m],L"f",position=(xpos,0),textsize=20)
        text!(axi2_PC[m],L"\lambda",position=(xpos,40),textsize=20)


        #text!(axi1_EOF[m],string.(round.(sum(protofreq[index]),digits=2)),position=(-200,0),textsize=16)
        text!(axi2_PC[m],string.(round.(sum(lambda[indices]),digits=3)),position=(xpos,-30),textsize=12,color=:blue)

    end


    #colsize!(F.layout, 1, Fixed(450))
    #colsize!(F.layout, 2, Fixed(450))

    #rowsize!(F.layout, 1, Fixed(300))
    #rowsize!(F.layout, 2, Fixed(300))
    #rowsize!(F.layout, 3, Fixed(100))



    for a in [axi1_EOF...,axi2_PC...]#axi3_RC...
        #hidexdecorations!(a)
        hideydecorations!(a)
        hidespines!(a)
    end

    hidexdecorations!(axi1_EOF[1],grid=false)
    hidexdecorations!(axi2_PC[1],grid=false)
    text!(axi1_EOF[1],"SSA",position=(-1.4,-2),textsize=16,rotation=pi/2)
    text!(axi1_EOF[2],"NLSA",position=(-1.4,-2),textsize=16,rotation=pi/2)

    text!(axi1_EOF[1],"(a)",position=(xpos+0.8,0))
    text!(axi1_EOF[2],"(b)",position=(xpos+0.8,0))
    text!(axi2_PC[1],"(c)",position=(xpos+0.8,0))
    text!(axi2_PC[2],"(d)",position=(xpos+0.8,0))



    save(savedirname*".png",F)
end


for spot in [3,6,19]
    component_variance_investigation(spot,dir*"pic/mostvar_$spot")
end


#=============================================================================#

"""
quick FFT comparison
"""
begin
    # Number of points
    N= 5478 
    # Sample period
    Ts = 1 / 365.25
    # Start time 
    t0 = 0 
    tmax = t0 + (N-1) * Ts
    # time coordinate
    t = t0:Ts:tmax

    # signal 
    signal = load(create_file_list(outdir,method,W,vari,preproc)[spot])["signal"]


    # Fourier Transform of it 
    F = fft(signal) |> fftshift
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

    # plots 
    #f,ax,s = plot(t, signal, title = "Signal")
    f,ax,s = lines(freqs[2740:end], abs.(F)[2740:end], color=:black,axis = (limits=(10^-1.5, 10^2.5,10^-0.5,10^4),xscale=log10,yscale=log10,))

    harm = list_pure_harmonics(7)

    vlines!(ax,harm,color="grey")


end

#=============================================================================#

"""
the big component investigation plot WITH FFT !
"""

spot = 19
outdir="/net/scratch/lschulz/fluxfullset/"
meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
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

xtic = round.(Array([1/7:1/7:5/7...,[1,2,3,4,7] ...]),digits=2)

function big_figure_fourier(spot,savedirname,name)
    F = Figure(resolution= (640,440))
    ax_hist = Axis(F[1,1],xscale=log10,ylabel=L"\sum \lambda",xlabel="f [a⁻¹]",yscale=log10, limits=(1/9, 10^1.5,10^-4,0.9),xticks=xtic)
    #,subtitle = name,titlealign = :left,
    #ax_text = Axis(F[1,2])
    ax_both = Axis(F[2,1],yscale=log10,xscale=log10,ylabel=L"\lambda",xlabel="f [a⁻¹]",limits=(1/9,10^1.5,10^-4,0.9),xticks=xtic,xticklabelrotation=pi/2)
    ax_spec = Axis(F[2,2],xlabel="k",ylabel=L"\lambda",yscale=log10,limits=(0,49,10^-4,0.9),xticklabelrotation=pi/2)

    colsize!(F.layout, 1, Fixed(320))
    colsize!(F.layout, 2, Fixed(80))

    rowsize!(F.layout, 1, Fixed(100))
    rowsize!(F.layout, 2, Fixed(200))

    harm = list_pure_harmonics(7)

    vlines!(ax_both,harm[end-6:end],color="grey")
    #vlines!(ax_hist,harm,color="grey")


    signal = load(create_file_list(outdir,"ssa",W,vari,preproc)[spot])["signal"]
    Four = fft(signal) |> fftshift
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
    Four ./= maximum(abs.(Four))
    lines!(ax_hist,freqs[2740:end], abs.(Four)[2740:end], color=:grey,label = "normalized FFT")


    for (i,method) = enumerate(["ssa","diff"])


        r_GPP = single_robustness(outdir,method,W,vari,preproc,kappa,spot,yearsamples) #RC,lambda,freq
        RC = r_GPP[1]
        lambda = r_GPP[2]
        protofreq = r_GPP[3]
        scatter!(ax_both,protofreq,lambda,marker=[:cross,:xcross][i],markersize=16,markerstrokewidth=4)

        h_weights =fit(Histogram,protofreq,weights(lambda),bins).weights
        lines!(ax_hist,bins[1:end-1],10^-5 .+ h_weights,label=["SSA","NLSA"][i])

        lines!(ax_spec,1:48,lambda)
    end


    hidespines!(ax_both, :t, :r)

    Legend(F[1,2],ax_hist,"Methods")



    hidespines!(ax_hist, :t, :r, :b)
    hidespines!(ax_spec, :t, :r)

    hidexdecorations!(ax_hist,ticks=false,grid=false)


    linkxaxes!(ax_hist,ax_both)
    linkyaxes!(ax_both,ax_spec)

    txt1 = "(a)"
    txt2 = "(b)"
    txt3 = "(c)"

    text!(ax_both,txt2,position=(1/7,0.15))
    text!(ax_spec,txt3,position=(5,0.15))
    text!(ax_hist,txt1,position=(1/7,0.15))


    save(savedirname*".png",F)
end

for (i,spot) in enumerate([3,6,19])
    name = ["Vielsalm MF","Oensingen CRO","Monte Bondone GRA"][i]
    big_figure_fourier(spot,dir*"pic/bigf_$spot",name)
end

#=============================================================================#

"""
the W stability plot
"""

function W_stab()
    
    #need LL 
    #need table_content function
    #needs proper kappa
    Ws = [5,5+1/4,5+1/3,5+1/2,5+2/3,5+3/4,6,6+1/4,6+1/3,6+1/2,6+2/3,6+3/4,7]
    f = Array(1:7)
    tabulus = [table_content(method,spot,W) for method=["ssa","diff"], spot=[3,19,6],W in Int.(floor.(Ws .*365.25))]
    #[rc_l,e_l,mse,e_bias,fast,seas,slow,harmonic_power,overtone_power]

    ytick_values = [Array(1:4)./ 5 ...,[1,2,3,4,7]...]
    ytick_labels = ["1/W", "2/W", "3/W", "4/W", "1","2","3","4","7"]
    cmap = cgrad(:summer,rev=true)

    for (s,spot) in enumerate([3,19,6])
        name = ["Vielsalm \n MF","Oensingen \n CRO","Monte \n Bondone \n GRA"][s]
        F = Figure(resolution=(850,400))

        for (i,method) in enumerate(["ssa","diff"])
            GL = F[1,1+i] = GridLayout()


            data = LL[i,:,:,spot]
            tit = ["SSA","NLSA"][i]
            ax = Axis(GL[1,1],limits=(Ws[1],Ws[end],1/7,7),yscale=log10,yticks=(ytick_values,ytick_labels),xlabel="W [a]",ylabel="f [a⁻¹]",xticks=round.(Ws,digits=2),title=tit,titlealign=:left)  
            #ax2 = Axis(F[2,i],xlabel=L"W [a]",limits=(Ws[1],Ws[end],0.4,0.8),ylabel=L" \lambda")
            ax3 = Axis(GL[2,1],xlabel="W [a]",limits=(Ws[1],Ws[end],0,1),xticks=round.(Ws,digits=2),xticklabelrotation=45,ylabel=L" \lambda") 


            hidespines!(ax, :t, :r, :b)
            #hidespines!(ax2, :t, :r,:b)
            hidespines!(ax3, :t, :r)


            #hidexdecorations!(ax2)
            hidexdecorations!(ax,ticks=false,grid=false)


            [lines!(ax,Ws,i ./Ws,color=:grey) for i=1:6]

            hlines!(ax,1,color=:grey)
            hlines!(ax,2,color=:grey)
            hlines!(ax,3,color=:grey)
            hlines!(ax,4,color=:grey)

            hm = contourf!(ax,Ws,bins[1:end-1],data .+ 0.0,colormap=cmap,lowclip=:white,highclip=:black,levels = [5*10^-4,10^-3,5*10^-3,10^-2,5*10^-2,10^-1,0.2,0.3,0.4,0.5])

            
            rc_l    =[i[1] .+ 10^-7 for i in tabulus[i,s,:]]
            #rc_l = sum(data,dims=2)[:]
            fast    =[i[5] .+ 10^-7 for i in tabulus[i,s,:]]
            seas    =[i[6] .+ 10^-7 for i in tabulus[i,s,:]]
            slow    =[i[7] .+ 10^-7 for i in tabulus[i,s,:]]
            harm    =[i[8] .+ 10^-7 for i in tabulus[i,s,:]]
            overt   =[i[9] .+ 10^-7 for i in tabulus[i,s,:]]

            #lines!(ax2,Ws,rc_l,label=L"\sum \lambda",color=:red)
            #lines!(ax2,Ws,fast,label="fast")
            #lines!(ax2,Ws,seas,label="seas")
            #lines!(ax2,Ws,slow,label="slow")
            lines!(ax3,Ws,harm + overt,label=L"\sum \lambda",color=:red)
            lines!(ax3,Ws,harm,label="harmonic")
            lines!(ax3,Ws,overt,label="residual")

            #Legend(F[2,3],ax2,"total")
            if i == 2
                Legend(GL[2,2],ax3),#"Variance")
                Colorbar(GL[1,2],hm,label = L"\sum \lambda")
                #text!(a,name,position = (7,0.1))
            end

            rowsize!(GL, 1,200)

            rowsize!(GL, 2, 50)

            colsize!(GL, 1, [270,280][i])


            colgap!(GL, 10)
            rowgap!(GL, 20)

            linkxaxes!(ax,ax3)
        end

        #a = Axis(F[1,1],title=name,titlealign=:right)
        #hidedecorations!(a)
        #hidespines!(a)

    save(dir*"/pic/Wstab_$spot.png",F)
    end

end


#=============================================================================#

"""
the W stability plot for LAIF
"""

function W_stab_LAI()
    
    #need LL 
    #need table_content function
    Ws = [15,16,17]
    f = Array(1:8)
    tabulus = [table_content(method,spot,W) for method=["ssa","diff"], spot=spots,W in Int.(floor.(Ws .*24))]
    #[rc_l,e_l,mse,e_bias,fast,seas,slow,harmonic_power,overtone_power]

    for (s,spot) in enumerate(spots)
        F = Figure(resolution=(1400,600))

        for (i,method) in enumerate(["ssa","diff"])
            GL = F[1,i] = GridLayout()


            data = LL[i,:,:,s]
            
            ax = Axis(GL[1,1],limits=(Ws[1],Ws[end],10^-1.5,10^1.5),yscale=log10,yticks=f,xlabel=L"W [a]",ylabel=L"f [a^{-1}]",xticks=round.(Ws,digits=2))  
            #ax2 = Axis(F[2,i],xlabel=L"W [a]",limits=(Ws[1],Ws[end],0.4,0.8),ylabel=L" \lambda")
            ax3 = Axis(GL[2,1],xlabel=L"W [a]",limits=(Ws[1],Ws[end],0,1),xticks=round.(Ws,digits=2),xticklabelrotation=45,ylabel=L" \lambda") 


            hidespines!(ax, :t, :r, :b)
            #hidespines!(ax2, :t, :r,:b)
            hidespines!(ax3, :t, :r)


            #hidexdecorations!(ax2)
            hidexdecorations!(ax,ticks=false,grid=false)


            [lines!(ax,Ws,i ./Ws,color=:grey) for i=1:6]

            hm = contourf!(ax,Ws,bins[1:end-1],data,colormap=:berlin,lowclip=:white,levels = [10^-6,10^-5,10^-4,10^-3,0.01,0.025,0.05, 0.1,0.25,0.5,0.75])

            
            rc_l    =[i[1] .+ 10^-7 for i in tabulus[i,s,:]]
            #rc_l = sum(data,dims=2)[:]
            fast    =[i[5] .+ 10^-7 for i in tabulus[i,s,:]]
            seas    =[i[6] .+ 10^-7 for i in tabulus[i,s,:]]
            slow    =[i[7] .+ 10^-7 for i in tabulus[i,s,:]]
            harm    =[i[8] .+ 10^-7 for i in tabulus[i,s,:]]
            overt   =[i[9] .+ 10^-7 for i in tabulus[i,s,:]]

            #lines!(ax2,Ws,rc_l,label=L"\sum \lambda",color=:red)
            #lines!(ax2,Ws,fast,label="fast")
            #lines!(ax2,Ws,seas,label="seas")
            #lines!(ax2,Ws,slow,label="slow")
            lines!(ax3,Ws,harm + overt,label=L"\sum \lambda",color=:red)
            lines!(ax3,Ws,harm,label="harm")
            lines!(ax3,Ws,overt,label="overt")

            #Legend(F[2,3],ax2,"total")
            if i == 2
                Legend(GL[2,2],ax3,"distribution")
                Colorbar(GL[1,2],hm,label = L"\sum \lambda")
            end
            rowsize!(GL, 1,400)
            #rowsize!(F.layout, 2, Fixed(100))
            rowsize!(GL, 2, 100)

            colsize!(GL, 1, [500,600][i])


            colgap!(GL, 10)
            rowgap!(GL, 10)
            linkxaxes!(ax,ax3)
        end



    save(dir*"/W_stab_LAI/$spot.png",F)
    end

end


#=============================================================================#

"""
the big component investigation plot WITH FFT ! LAI !
"""

spot = 19
outdir="/net/scratch/lschulz/fluxfullset/"
meta = load("/net/scratch/lschulz/fluxnetfullset/fullset_15a_gpp_nee_reco_ts.jld2")["meta"]
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

function big_figure_fourier_LAI(spot,savedirname)
    F = Figure(resolution= (1000,640))
    ax_hist = Axis(F[1,1],xscale=log10,ylabel=L"\sum \lambda",xlabel=L"f",yscale=log10, limits=(10^-1.5, 10^1.5,10^-7,10^0.5),xticks=Array(1:7))
    #ax_text = Axis(F[1,2])
    ax_both = Axis(F[2,1],yscale=log10,xscale=log10,ylabel=L"\lambda",xlabel=L"f",limits=(10^-1.5,10^1.5,10^-7,0.9),xticks=Array(1:7))
    ax_spec = Axis(F[2,2],xlabel="k",ylabel=L"\lambda",yscale=log10,limits=(0,49,10^-7,0.9))


    colsize!(F.layout, 1, Fixed(600))
    colsize!(F.layout, 2, Fixed(150))

    rowsize!(F.layout, 1, Fixed(200))
    rowsize!(F.layout, 2, Fixed(300))

    harm = list_pure_harmonics(7)

    vlines!(ax_both,harm,color="grey")
    #vlines!(ax_hist,harm,color="grey")


    signal = load(create_file_list(outdir,"ssa",W,vari,preproc)[spot])["signal"]
    Four = fft(signal) |> fftshift
    freqs = fftfreq(length(t), 1.0/Ts) |> fftshift
    Four ./= maximum(abs.(Four))
    lines!(ax_hist,freqs[2740:end], abs.(Four)[2740:end], color=:grey,label = "normalized FFT")


    for (i,method) = enumerate(["ssa","diff"])


        r_GPP = single_robustness(outdir,method,W,vari,preproc,kappa,spot,yearsamples) #RC,lambda,freq
        RC = r_GPP[1]
        lambda = r_GPP[2]
        protofreq = r_GPP[3]
        scatter!(ax_both,protofreq,lambda,marker=[:cross,:xcross][i],markersize=16,markerstrokewidth=4)

        h_weights =fit(Histogram,protofreq,weights(lambda),bins).weights
        lines!(ax_hist,bins[1:end-1],10^-7 .+ h_weights,label=["SSA","NLSA"][i])

        lines!(ax_spec,1:48,lambda)
    end


    hidespines!(ax_both, :t, :r)

    Legend(F[1,2],ax_hist,"Methods")



    hidespines!(ax_hist, :t, :r, :b)
    hidespines!(ax_spec, :t, :r)

    hidexdecorations!(ax_hist,ticks=false,grid=false)


    linkxaxes!(ax_hist,ax_both)
    linkyaxes!(ax_both,ax_spec)

    save(savedirname*".png",F)
end

for (s,spot) in enumerate(spots)
    big_figure_fourier_LAI(spot,dir*"bigf/$spot")
end

#=============================================================================#

"""
standard coupling GPP TS
"""
preproc = "raw."

for spot = 1:20
    F = Figure(resolution=(1200,400))

    for (i,method) = enumerate(["ssa","diff"])
        class = meta[spot,"IGBP_class"]
        tit = "$(spot)_$(class)_$(method)_$(preproc)_GPP"

        #first read in GPP
        r_GPP = single_robustness(outdir,method,W,1,preproc,kappa,spot,yearsamples) #RC,lambda,freq
        protofrequencies = reshape(r_GPP[3],48,1)
        RC_GPP = sum(reshape(r_GPP[1],5478,48,1),dims=2)[:]
        lambdas_GPP = reshape(r_GPP[2],48,1)
        signal_GPP = load(create_file_list(outdir,"ssa",W,1,preproc)[spot])["signal"]

        #then read in TS
        r_TS = single_robustness(outdir,method,W,4,preproc,kappa,spot,yearsamples) #RC,lambda,freq
        protofrequencies = reshape(r_GPP[3],48,1)
        RC_TS = sum(reshape(r_TS[1],5478,48,1),dims=2)[:]
        lambdas_TS = reshape(r_TS[2],48,1)
        signal_TS = load(create_file_list(outdir,"ssa",W,4,preproc)[spot])["signal"]

        ax = Axis(F[1,i],limits=(-pi,pi,-1/2,2),xlabel=L"\Delta",ylabel=L"p",title=tit*" $band")

        bins,hist = phase_diff_hist(signal_GPP,signal_TS) #gives bins,cauchy_weights
        lines!(ax,bins[1:end],hist,label="Observation",color = :red)

        bins,hist = phase_diff_hist(RC_GPP,RC_TS) #gives bins,cauchy_weights
        lines!(ax,bins[1:end],hist,label="RC",color = :blue)

        bins,hist = phase_diff_hist(Float32.(protophase(RC_GPP,1)),Float32.(protophase(RC_TS,1))) #gives bins,cauchy_weights
        lines!(ax,bins[1:end],hist,label=L"\phi",color = :black)


    end

    Legend(F[1,3],ax)

    save(dir*"cop/$spot.png",F)
end

#=============================================================================#
"""
the slow modes plot 
"""

method="diff"
preproc="raw."
spot=19
N = 5478
k = 48
tima = 2005 .+ (1:N) ./365.25
vari = 1
pos= (2005,2.5)
pos2 = (2005+7,3)
#tic = LinearTicks(7.5)
tic = 2005:2:2019
function harmonic_indices(protofreq)
    harmonics = unique(round.(vcat([i for i=1:7],[1/i for i=1:7]),digits=6))
    L = []
    epsi = 0.01
    for harmonic in harmonics
        indices = findall(x->harmonic*(1-epsi)<=x<=harmonic*(1+epsi), protofreq)
        L = append!(L,indices)
    end
    return L 
end



limits=(tima[1],tima[end],-1.5,4)


for (s,spot) in enumerate([3,19,6])

    f = Figure(resolution=(800,450))
    
    for (m,method) in enumerate(["ssa","diff"])


        Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
        lambda,protofreq,RC = extract_from_single_file(Filename,yearsamples,N,k)
        protofreq = protofreq[:]
        GL = f[1,m] = GridLayout()

        #x<0.85 this is slow!
        ax = Axis(GL[1,1],xticks=tic,limits=limits)
        hidedecorations!(ax)
        hidespines!(ax)

        lines!(ax,tima,load(Filename)["signal"] .+ 0.0,color=:grey)
        indices = findall(x->x<0.85,protofreq)
        #indices = findall(x->x<1.0,protofreq)

        lines!(ax,tima,sum(RC[:,indices],dims=2)[:] .+ 0.0 ,color=:blue)
        text!(ax,"slow",position=pos)
        text!(ax,"$(["SSA","NLSA"][m])",position=pos2)

        #0.85 < x 1.15 this is annual
        ax = Axis(GL[2,1],xticks=tic,limits=limits)
        hidedecorations!(ax)
        hidespines!(ax)

        lines!(ax,tima,load(Filename)["signal"] .- 0.0,color=:grey)
        indices = findall(x->0.85<x<1.15,protofreq)
        #indices = findall(x->x==1.0,protofreq)

        lines!(ax,tima,sum(RC[:,indices],dims=2)[:] .- 0.0 ,color=:orange)
        text!(ax,"annual",position=pos)

        #1.15 < x this is fast
        ax = Axis(GL[3,1],xticks=tic,limits=limits,xlabel="t [a]")
        hideydecorations!(ax)
        hidespines!(ax)

        lines!(ax,tima,load(Filename)["signal"] .- 0.0,color=:grey)
        indices = findall(x->x>1.15,protofreq)
        #indices = findall(x->x>1.0,protofreq)

        lines!(ax,tima,sum(RC[:,indices],dims=2)[:] .- 0.0 ,color=:red)
        text!(ax,"fast",position=pos)

        rowgap!(GL,Fixed(0))
        colgap!(GL,Fixed(0))
    end

    save(dir*"pic/fa_$spot.png",f)


#=
    #FFT
    begin
            
        N= 5478 
        # Sample period
        Ts = 1 / 365.25
        # Start time 
        t0 = 0 
        tmax = t0 + (N-1) * Ts
        # time coordinate
        t = t0:Ts:tmax

        GL = f[1,3] = GridLayout()

        Filename = create_file_list(outdir,"ssa",W,vari,preproc)[spot]
        signal = load(Filename)["signal"]
        Four = fft(signal) |> fftshift
        freqs = fftfreq(length(t), 1.0/Ts) |> fftshift

        protofreq = freqs
        stren = abs.(Four)

        #x<0.85 this is slow!
        ax = Axis(GL[1,1],xticks=1:15,limits=limits)
        hidedecorations!(ax)
        hidespines!(ax)
        lines!(ax,tima,load(Filename)["signal"] .+ 0.0,color=:grey)

        indices = findall(x->0.0<=abs(x)<0.85,protofreq[2439:end])
        brrt = copy(Four)
        brrt[Not(indices)] .= 0.0

        lines!(ax,tima, abs.(ifft(brrt)) ,color=:blue)
        text!(ax,"slow",position1=[1,1])

        #0.85 < x 1.15 this is annual
        ax = Axis(GL[2,1],xticks=1:15,limits=limits)
        hidedecorations!(ax)
        hidespines!(ax)
        lines!(ax,tima,load(Filename)["signal"] .- 0.0,color=:grey)
        indices = findall(x->0.85<abs(x)<1.15,protofreq)
        brrt = copy(Four)
        brrt[Not(indices)] .= 0.0
        lines!(ax,tima,abs.(ifft(brrt)) ,color=:orange)
        text!(ax,"season",position1=[1,1])

        #1.15 < x this is fast
        ax = Axis(GL[3,1],xticks=1:15,limits=limits)
        hidedecorations!(ax)
        hidespines!(ax)
        lines!(ax,tima,load(Filename)["signal"] .- 0.0,color=:grey)
        indices = findall(x->abs.(x)>1.15,protofreq)
        brrt = copy(Four)
        brrt[Not(indices)] .= 0.0
        lines!(ax,tima,abs.(ifft(brrt)) ,color=:red)
        text!(ax,"fast",position1=[1,1])

        # harmonic indices
        ax = Axis(GL[4,1],xticks=1:15,limits=limits)
        hidedecorations!(ax)
        hidespines!(ax)
        lines!(ax,tima,load(Filename)["signal"] .- 0.0,color=:grey)
        indices = harmonic_indices(abs.(protofreq))
        brrt = copy(Four)
        brrt[Not(indices)] .= 0.0
        lines!(ax,tima,abs.(ifft(brrt)) ,color=:green)
        text!(ax,"overT",position1=[1,1])

        # unharmonic indices
        ax = Axis(GL[5,1],xticks=1:15,limits=limits)
        hidedecorations!(ax)
        hidespines!(ax)
        lines!(ax,tima,load(Filename)["signal"] .- 0.0,color=:grey)
        indices = Array(1:5478)[Not(indices)]
        brrt = copy(Four)
        brrt[Not(indices)] .= 0.0
        lines!(ax,tima,abs.(ifft(brrt)) ,color=:purple)
        text!(ax,"wo overT",position1=[1,1])

        rowgap!(GL,Fixed(0))
        colgap!(GL,Fixed(0))
    end
=#
        """
        partials
        """
    f = Figure(resolution=(800,300))

    for (m,method) in enumerate(["ssa","diff"])


        Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
        lambda,protofreq,RC = extract_from_single_file(Filename,yearsamples,N,k)
        protofreq = protofreq[:]
        GL = f[1,m] = GridLayout()

        # harmonic indices
        ax = Axis(GL[1,1],xticks=tic,limits=limits)
        hidedecorations!(ax)
        hidespines!(ax)

        lines!(ax,tima,load(Filename)["signal"] .- 0.0,color=:grey)
        indices = harmonic_indices(protofreq)
        lines!(ax,tima,sum(RC[:,indices],dims=2)[:] .- 0.0 ,color=:green)
        text!(ax,"harmonic",position=pos)
        text!(ax,"$(["SSA","NLSA"][m])",position=pos2)
        # unharmonic indices
        ax = Axis(GL[2,1],xticks=tic,limits=limits,xlabel="t [a]")
        hideydecorations!(ax)
        hidespines!(ax)

        lines!(ax,tima,load(Filename)["signal"] .- 0.0,color=:grey)
        indices = Array(1:48)[Not(indices)]
        lines!(ax,tima,sum(RC[:,indices],dims=2)[:] .- 0.0 ,color=:purple)
        text!(ax,"residual",position=pos)

        rowgap!(GL,Fixed(0))
        colgap!(GL,Fixed(0))
    end

    save(dir*"pic/pa_$spot.png",f)


end

