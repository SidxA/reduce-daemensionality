"""
from toolbox
"""

    function pairing_full(EOFs,epsi_T,epsi_A,print)
        pairs = Tuple[]
        periods = Float64[]
        k = size(EOFs)[2]
        indices = Array(1:k)
        lags = 1:Int(size(EOFs)[1]-1)
        for i=1:k, j= 1:k
            x=EOFs[:,i]
            y=EOFs[:,j]
            c = crosscor(x,y,lags)
        if print  println("$i \t $j \t",length(argmaxima(c))) end
            if length(argmaxima(c)) >= 3
                cc = diff(argmaxima(c))
                if abs(cc[1]-cc[2])<=epsi_T && all(c[argmaxima(c)][1:2] .> epsi_A)
                    if print println("$i \t $j found T $(cc[1])") end
                    pairs = push!(pairs,(i,j))
                    periods = push!(periods,mean(cc[1:2]))
                end
            end
        end
        return pairs,periods
    end


    function interval_velocity(signal::Vector,L)
        ind = L:L:length(signal)
        n_intervals = length(ind)
        flips_zeros = Array{Float64}(undef,0)
        flips_maxima = Array{Float64}(undef,0)
        for interval in [signal[i-L+1:i] for i in ind]
            flips_zeros = push!(flips_zeros,count_sign_flip(interval)/L*2*pi)
            flips_maxima = push!(flips_maxima,2*count_maxima(interval)/L*2*pi)
        end
        return hcat(flips_zeros,flips_maxima)
    end


    #pricipal components
    function pc!(X,rho,P,M)
        return unfold([Array(sum([X[t+j-1]*rho[j,:] for j in 1:M],dims = 1)...) for t in 1:P])
    end

    #make 2d from nested 2d array
    unfold(multiArr) = reshape(hcat(multiArr...),size(multiArr)[1],size(multiArr[1])[1])

    #load
    function load_data()
        name = "testdata/FLX_US-PFa_FLUXNET2015_FULLSET_DD_1995-2014_1-4.csv"
        df = DataFrame(CSV.File(name,delim=",",header=2))
        #remove timestamp
        df = df[:,2:end]
        #remove constant rows
        il=[]
        for i in 1:ncol(df)
            if df[2,i]==df[3,i]
                il = push!(il,false)
            elseif any(i -> abs(i) > 500,df[:,i])
                il = push!(il,false)
            else
                il = push!(il,true)
            end
        end

        df = Array(df)
        il = Array(Bool.(il))
        df = df[:,il]

        return df[:,20]
    end


    # windows at W/2 along the data cuts off the rest
    function embed_lag_centered(x,win_l)
        T = length(x)
        P = win_l
        N = floor(Int,(2(T+1)/P))-1
        data = zeros(N,P)
        ind = 1
        for t in 1:N
            ind_l = ind + P - 1
            if ind_l <= T
                data[t,:] = x[ind:ind_l]
                ind += Int(P/2)
            else
                break
            end
        end
        return data
    end

    #simple enough boundary hack for local ssa
    boundary_replicate(data) = append!([vcat(data) for i=1:5]...)

    #----------------------------------------------------------------------
    # global SSA
    #----------------------------------------------------------------------
    #returns ONLY RC NOW 
    #rc, possiblysignal,eof,pc,rc,lambdas
    # this does not centralize itself now neither does it do the boundaries
    function SSA_g(data::Vector,W,k,key=:auto)

        N=length(data)
        P = N - W +1

        emb_data = embed_lag(data,W)
        pca = fit(PCA,emb_data,method=key,pratio=1.0,maxoutdim=k)

        EOFs = projection(pca)
        PC = hcat([pc!(data,EOFs[:,i],P,W) for i in 1:k]...)'
        RC = hcat([reconstructor(PC[i,:],EOFs[:,i],N,W) for i in 1:k]...)

        #return cutoff(data),EOFs,PC,RC,principalvars(pca)
        return RC
    end

    #----------------------------------------------------------------------
    # global Diffusion map
    #----------------------------------------------------------------------
    #returns rc
    # this does not centralize itself now neither does it do the boundaries

    function Diff_g(signal::Vector,W,k)
        t=1
        α=0.5
        ɛ=512

        N=length(signal)
        P = N - W +1

        diff = fit(DiffMap,embed_lag(signal,W)',maxoutdim=k,t=t, α=α, ɛ=ɛ)

        EOFs = ManifoldLearning.transform(diff)'
        PC = hcat([pc!(signal,EOFs[:,i],P,W) for i in 1:k]...)'
        RC = hcat([reconstructor(PC[i,:],EOFs[:,i],N,W) for i in 1:k]...)

        return RC
    end


    #----------------------------------------------------------------------
    # module linking the global methods and centralize and boundary
    #----------------------------------------------------------------------

    #simple enough cutoff ~ remove additional N points
    function cutoff(signal::Vector,N_bound)
        return signal[N_bound+1:end-N_bound]
    end

    #simple enough boundary hack ~ put additional N points
    function boundary(data::Vector,N_bound)
        first = data[1:N_bound]
        last = data[end-N_bound+1:end]
        return append!(first,data,last)
    end

    function components_g(signal::Vector,W,k,N_bound,Bool_centre)
        signal = boundary(signal,N_bound)
        if Bool_centre
            signal = centralizer(signal)
        end
        rc_gSSA = SSA_g(signal,W,k)
        rc_gDiff = Diff_g(signal,W,k)
        return hcat([cutoff(i,N_bound) for i in
                [signal,
                [rc_gSSA[:,c_ind] for c_ind in 1:k]...,
                [rc_gDiff[:,c_ind] for c_ind in 1:k]...]
                ]...)
    end

"""
from somehere else
"""


"""
from somehere else
"""


"""
plotting stuff
"""

begin

    """
    mapping to a picture of europe with slightly tweaked locations
    """
    function lat_long()#return x,y
        #first use original locations
        table = CSV.read("logs/KW_29/"*"fluxsites.csv",DataFrame,header=0)
        filelist = readdir("/net/data/Fluxnet/FLUXNET2020-ICOS-WarmWinter/")[5:end-3]
        nameslist = String[]
        for file in filelist
            filenamesplit = split(file,"_")
            years = split(filenamesplit[end-1],"-")
            period = parse(Int,years[2])-parse(Int,years[1])
            if filenamesplit[4] == "ERAI" && filenamesplit[5] == "DD" && filenamesplit[end] =="beta-3.csv" && period > 30
            nameslist = push!(nameslist,filenamesplit[2])
            end
        end

        loc_ind=1:71
        i = [findfirst(x->x[1:2]==name[1:2],table[:,1]) for name in nameslist]
        i = [isnothing(ix) ? 0 : ix for ix in i]
        #israel would be around 15,40
        #island would be around 5, 50
        i[55] = 121
        i[56] = 27
        #loc_ind = loc_ind[findall(!iszero,i)]
        i = i[findall(!iszero,i)]

        #the stupid manual list
        tsml = 
        reshape(
            vcat(
        [
            #belgium
            [0,0+2],
            [2,0+2],
            [0,2+2],
            [-2,0+2],
            [0,-2+2],
            [-2,2+2],
            #switz?
            [0+2,0-2],
            [2+2,0-2],
            [0+2,2-2],
            [-2+2,0-2],
            [0+2,-2-2],
            [2+2,2-2],
            #czesk
            [0+2,0-2],
            [2+2,0-2],
            [0+2,2-2],
            [-2+2,0-2],
            [0+2,-2-2],
            [-2+2,2-2],
            #ger
            [0-2,0-1],
            [2-2,0-1],
            [0-2,2-1],
            [-2-2,0-1],
            [0-2,-2-1],
            [-2-2,2-1],
            [4-2,2-1],
            [4-2,0-1],
            [0-2,-4-1],
            [-4-2,0-1],
            [-4-2,-4-1],
            [-2-2,-4-1],
            #denmark
            [-2,-2+7],
            [-3,-4+7],
            #spain
            [-1,3],
            [-5,3],
            [-3,5],
            [-5,3],
            [-3,2],
            [-5,5],
            #finland
            [0,0],
            [2,0],
            [0,2],
            [-2,0],
            [0,-2],
            [-2,-2],
            #france
            [0,0-3],
            [2,0-3],
            [0,2-3],
            [-2,0-3],
            [0,-2-3],
            [-2,2-3],
            [2,2-3],
            [3,0-3],
            #gf
            [+45,+30],
            #gl
            [+45,-6],
            #ie
            [+10,0],
            #il
            [-5,+15],
            #italy
            [0,0],
            [2,0],
            [0,2],
            [-2,0],
            [0,-2],
            [-2,2],
            [4,2],
            [4,0],
            #russland
            [-2,-2],
            [-3,-4],
            #SE
            [0-4,0-2],
            [2-4,0-2],
            [0-4,2-2],
            [-2-4,0-2],
            [0-4,-2-2],
        ]...),
        2,71)'

        x = table[i,6]
        y = table[i,5]
        x = x .+ tsml[:,1]
        y = y .+ tsml[:,2]
        x .+= rand(71).-0.5
        y .+= rand(71).-0.5
        return x,y
    end
    """
    indices pairs of links that are inside k and epsilon
    """
    function find_links(x,y,k,r)
        #neigbors
        data = hcat(x,y)'
        tree = KDTree(data)
        idxs_k,dists = knn(tree,data,k+1,true)
        #OR for balls of radius r
        idxs_r = inrange(tree, data, r)
        L = []
        for loc_ind in 1:71
            arr = intersect(idxs_k[loc_ind],idxs_r[loc_ind])
            self = findall(x->x==loc_ind,arr)
            arr = deleteat!(arr,self)
            for x in arr
                if isempty(findall(y->y==[x,loc_ind],L)) 
                    L = push!(L,[loc_ind,x])
                end
            end
        end
        return Matrix(hcat(L...)')
    end

    function find_links2(x,y,k,r)
        #neigbors
        data = hcat(x,y)'
        tree = KDTree(data)
        idxs_k,dists = knn(tree,data,k+1,true)
        #OR for balls of radius r
        idxs_r = inrange(tree, data, r)
        L = []
        for loc_ind in 1:71
            arr = intersect(idxs_k[loc_ind],idxs_r[loc_ind])
            self = findall(x->x==loc_ind,arr)
            arr = deleteat!(arr,self)
            for x in arr
                if isempty(findall(y->y==[x,loc_ind],L)) 
                    L = push!(L,[loc_ind,x])
                end
            end
        end
        return L
    end

    """
    complete individual spot plotting procedure
    for the color we need a list of the sorted indeces
    """
    img = FileIO.load("logs/KW_29/"*"worldmap.png")

    x,y =lat_long()

    function fluxscatter(x,y,loc_inds,inds,label,savedir)

        x = x[loc_inds]
        y = y[loc_inds]
        c = inds

        r = 15
        k = 7
        F = Figure()
        a = Axis(F[1,1],limits=(-10, 30, 30, 70))
        image!(a,[-180,180],[-90,90],rotr90(img),yflip=true)


        links = find_links(x,y,k,r)
        start = links[:,1]
        stop = links[:,2]
        for l in 1:size(links)[1]
            lines!(a,[x[start][l],x[stop][l]],[y[start][l],y[stop][l]],color=:black)
        end

        s = scatter!(a,vec(x),vec(y),markersize=30,color=c,colormap=:Spectral)
        Colorbar(F[1,2], s,label=label)

        save(dir*savedir*".png",F)
    end

    """
    plot coupling between spots for a single variable
    with complete phase diff we have rec_m = raw

    sort the spots
    """


    function list_all_links(loc_inds::Vector{Int64})
        L = []
        for i in loc_inds, j in loc_inds
            L = push!(L,[i,j])
        end
        return L
    end

    function later()
        H_diff = [spectral_entropy(datadir,1,loc_ind,"diff") for loc_ind in 1:71]

        x = sortperm(H_diff)

        raw = reshape(coupling_strength(list_all_links(Array(x)),datadir,"raw",1),71,71)
        diff = reshape(coupling_strength(list_all_links(Array(x)),datadir,"diff",1),71,71)
        ssa = reshape(coupling_strength(list_all_links(Array(x)),datadir,"ssa",1),71,71)

        joint_limits = (0, 2)  # here we pick the limits manually for simplicity instead of computing them

        fig, ax1, hm1 = heatmap(raw,  colorrange = joint_limits)
        ax2, hm2 = heatmap(fig[1, 2], ssa, colorrange = joint_limits)
        ax3, hm3 = heatmap(fig[2, 1], diff, colorrange = joint_limits)

        Colorbar(fig[:, end+1], hm1)

        save(dir*name*".png",fig)

    end

    """
    lets do the bands
    directly jump to diffusion
    """


    """
    plotting procedure to show the phase coupling in the whole fluxnet dataset
    by rec_m and season_ind, var_ind
    resulting plot into one directory for the overlook for me choose 1 big plot
    with 2 columns of spec_en and rec_err

    spectral entropy plot                   spec_en
    quadratic reconstruction error plot     rec_err

    each with the other sorted by its indices to show wether there is some corr

    the phasecoupling must be investigated with protophase per mode added together...
    can not put protophase on the sum of rc!
    for this we propably need some proper cutoff function, maybe even with normalization?

    ideally we would have some
    full reconstruction coupling between spots sorted by spec_en and rec_err
    Frequency histogram overview and bandstrengths plot for the 2 sortings of rec_err and spec_en
    bands coupling between spots sorted by spec_en and rec_err

    """

    function doit!()
        var_ind = 1
        season_ind = 1
        rec_m = "ssa"

        spec_en = [spectral_entropy(datadir,var_ind,loc_ind,rec_m,season_ind) for loc_ind in 1:71]
        spec_en_ind = sortperm(spec_en)
        spec_en_links = list_all_links(spec_en_ind)

        rec_err_season = [reconstruction_error(datadir,var_ind,loc_ind,rec_m,season_ind)[1] for loc_ind in 1:71]
        rec_err_deseason = [reconstruction_error(datadir,var_ind,loc_ind,rec_m,season_ind)[2] for loc_ind in 1:71]
        rec_err = rec_err_season .+ rec_err_deseason
        rec_err_ind = sortperm(rec_err)
        rec_err_links = list_all_links(rec_err_ind)

        Fig = Figure(resolution=(1800,900))

        ax1 = Axis(Fig[1,1])
        lines!(ax1,centralizer(sort(spec_en)),label="spec en")
        lines!(ax1,centralizer(rec_err_season[spec_en_ind]),label="season err")
        lines!(ax1,centralizer(rec_err_deseason[spec_en_ind]),label="deseason err")
        lines!(ax1,centralizer(rec_err[spec_en_ind]),label="rec err")
        axislegend(ax1, merge =true, unique = true,position =:rb)

        ax2 = Axis(Fig[1,2])
        lines!(ax2,centralizer(spec_en[rec_err_ind]),label="spec en")
        lines!(ax2,centralizer(rec_err_season[rec_err_ind]),label="season err")
        lines!(ax2,centralizer(rec_err_deseason[rec_err_ind]),label="deseason err")
        lines!(ax2,centralizer(sort(rec_err)),label="rec err")
        axislegend(ax2, merge =true, unique = true,position =:rb)

        """
        the heatmaps for the gamma
        """
        crange = (0,2.5)

        
        # the complete reconstructed coupling strength sorted by spec_en and rec_err
        all_modes = repeat([Array(1:48)],71)
        comp_cpl_by_spec_en = reshape(coupling_strength(spec_en_links,all_modes,datadir,rec_m,var_ind,season_ind),71,71)
        comp_cpl_by_rec_err = reshape(coupling_strength(rec_err_links,all_modes,datadir,rec_m,var_ind,season_ind),71,71)

        save(dir*"test.png",Fig)

        raw = reshape(coupling_strength(list_all_links(Array(x)),datadir,"raw",1),71,71)
        diff = reshape(coupling_strength(list_all_links(Array(x)),datadir,"diff",1),71,71)
        ssa = reshape(coupling_strength(list_all_links(Array(x)),datadir,"ssa",1),71,71)

        joint_limits = (0, 2)  # here we pick the limits manually for simplicity instead of computing them

        fig, ax1, hm1 = heatmap(raw,  colorrange = joint_limits)
        ax2, hm2 = heatmap(fig[1, 2], ssa, colorrange = joint_limits)
        ax3, hm3 = heatmap(fig[2, 1], diff, colorrange = joint_limits)

        Colorbar(fig[:, end+1], hm1)

        save(dir*name*".png",fig)

        signal_raw_1 = centralizer(load(datadir*"diff_$(loc_ind_1)_$(var_ind)_3.jld2")["signal"])
        ssa_1 = load(datadir*"ssa_$(loc_ind_1)_$(var_ind)_$(season_ind).jld2")["RC"]
        diff_1  = load(datadir*"diff_$(loc_ind_1)_$(var_ind)_$(season_ind).jld2")["RC"]

        signal_raw_2 = centralizer(load(datadir*"diff_$(loc_ind_2)_$(var_ind)_3.jld2")["signal"])
        ssa_2 = load(datadir*"ssa_$(loc_ind_2)_$(var_ind)_$(season_ind).jld2")["RC"]
        diff_2  = load(datadir*"diff_$(loc_ind_2)_$(var_ind)_$(season_ind).jld2")["RC"]

        data = hcat(cutoff(signal_raw_1,365),cutoff(sum(ssa_1,dims=2)[:,1],365),cutoff(sum(diff_1,dims=2)[:,1],365),cutoff(signal_raw_2,365),cutoff(sum(ssa_2,dims=2)[:,1],365),cutoff(sum(diff_2,dims=2)[:,1],365))
        protophase = hcat([protophase_by_varind_recm(datadir, loc_ind, season_ind, var_ind, rec_m) for loc_ind in [loc_ind_1,loc_ind_2], rec_m in ["raw","ssa","diff"]]...)

        bins = Array(-pi:0.01:pi)

        #phases

        h = [fit(Histogram,data[:,k],bins) for k in 1:6]
        [normalize(h[k],mode=:pdf) for k in 1:6]
        hist_d = Float32.(hcat([h[k].weights for k in 1:6]...))


        h = [fit(Histogram,protophase[:,k],bins) for k in 1:6]
        [normalize(h[k],mode=:pdf) for k in 1:6]
        hist_p = Float32.(hcat([h[k].weights for k in 1:6]...))

        # differences

        h = [fit(Histogram,data[:,k] .- data[:,3+k],bins) for k in 1:3]
        [normalize(h[k],mode=:pdf) for k in 1:3]
        hist_d_D = Float32.(hcat([h[k].weights for k in 1:3]...))


        h = [fit(Histogram,protophase[:,k] .- protophase[:,3+k],bins) for k in 1:3]
        [normalize(h[k],mode=:pdf) for k in 1:3]
        hist_p_D = Float32.(hcat([h[k].weights for k in 1:3]...))

        fig = Figure(resolution = (3600,1800))

        ax1 = Axis(fig[1,1],title="2 grey signals with red ssa reconstruction, blue diffmap")
        data[:,4:6] .+= 5
        sp1 = series!(ax1,data',color=[:grey,:red,:blue,:grey,:red,:blue])

        ax2 = Axis(fig[1,2],limits=(1, 1000, -10, 20),title="extracted protophases")
        protophase[:,4:6] .+= 10
        sp2 = series!(ax2,protophase',color=[:grey,:red,:blue,:grey,:red,:blue])

        ax3 = Axis(fig[2,1],title="raw signal phase histogram")
        hist_d[:,4:6] .+= 1.2*  maximum(hist_d[:,1:3])
        sp3 = series!(ax3,bins[1:end-1],hist_d',color=[:grey,:red,:blue,:grey,:red,:blue])

        ax4 = Axis(fig[2,2],title="protophase histogram")
        hist_p[:,4:6] .+= 1.2*  maximum(hist_p[:,1:3])
        sp4 = series!(ax4,bins[1:end-1],hist_p',color=[:grey,:red,:blue,:grey,:red,:blue])

        ax5 = Axis(fig[3,1],title = "raw signal phase difference histogram")
        hist_d_D[:,2] .+= 3
        hist_d_D[:,3] .+= 6
        sp5 = series!(ax5,bins[1:end-1],hist_d_D',color=[:grey,:red,:blue,:grey,:red,:blue])

        ax6 = Axis(fig[3,2],title = "protophase difference histogram")
        hist_p_D[:,2] .+= 3
        hist_p_D[:,3] .+= 6
        sp6 = series!(ax6,bins[1:end-1],hist_p_D',color=[:grey,:red,:blue,:grey,:red,:blue])

        save(dir*"test.png",fig)
    end

    """
    plot a small little thing to explain the different things that can happen in the reconstruction
    """
    function protophase()
        mode_indices=1:48
        var_ind = 1
        season_ind = 1
        rec_m = "ssa"
        loc_ind=rand(vec(1:71))

        RC = load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["RC"][:,mode_indices]
        signal = load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["signal"]
        lam = load(datadir*"$(rec_m)_$(loc_ind)_$(var_ind)_$(season_ind).jld2")["lambda"]
        ph = [protophase(RC[:,i]) for i=1:size(RC)[2]]
        pha = hcat(ph...)
        phas = centralizer(vec(sum(pha'.*lam,dims=1)))
        Sin = centralizer(sum([lam[i]*sin.(ph[i]) for i=1:48]))
        Cos = centralizer(sum([lam[i]*sin.(ph[i]) for i=1:48]))

        fig,ax,l = lines(signal,label="signal")
        lines!(ax,vec(sum(RC,dims=2)).+ 3,color=:red,label="reconstruction")
        lines!(ax,centralizer(Sin).+6,color=:green,label="sum proto-sin")
        #lines!(ax,centralizer(Cos).+9,color=:yellow,label="sum proto-cos")
        lines!(ax,centralizer(phas).+9,color=:black,label="sum proto-phase")
        axislegend(ax,position=:lb)
        save(dir*"proto.png",fig)
    end

    """
    plot 3 little figures
    to show the epsilon in the diff map
    """
    function diff_eps()
        #this needs to take out the loop!
        include("/net/home/lschulz/scripts/struct_big_flux.jl")

        N = 11895
        W = 5844
        k = 48

        d = matrixholder(N,W,k)

        put_data_by_season(d,wholedata[:,loc_ind,var_ind],season)
        season = 3

        put_epsilon(d)
        put_EOF_diff(d)
        calculate(d)

        stepnumber      = 16
        sampling_size   = 128
        it_number       = 8

        data = d.emb
        P = size(data)[2]
        object = epsilon(W,stepnumber,sampling_size)
        fit_eps_L = Vector{Float32}(undef,it_number)
        sample = rand(1:P,sampling_size,it_number)
        t = 3
        object.data_samples = data[:,sample[:,t]]
        for (i,eps) in enumerate(object.eps)
            for i in 1:sampling_size, j in 1:sampling_size
                object.Weight[i,j] = exp(- norm(object.data_samples[:,i] - object.data_samples[:,j])^2 / eps)
            end
            object.L[i] = sum(object.Weight)
        end

        p0 = ones(3)
        model(eps,p) = p[1] .* atan.(eps .- p[2]) .+ p[3]
        fit_eps_L[t] = 10^(coef(curve_fit(model, log10.(object.eps), log10.(object.L), p0))[2])

        F = Figure(resolution = (1800,900))

        #data manifold
        ax1 = Axis3(F[1,1],azimuth=0.1*pi)
        scatter!(ax1,RC[:,1],RC[:,2],RC[:,3],markersize=5)

        epsi = 10 .^ [0,1,2,3,4,5,6]
        elevate = Array(vec(1:7) .* -3)

        #epsilon saturation
        ax2 = Axis(F[1,2],xlabel="log10 eps",ylabel="log10 sum Weights")
        lines!(ax2,log10.(object.eps), log10.(object.L))
        vlines!(ax2,log10.(epsi),color=:grey)

        #different reconstructions
        ax3 = Axis(F[1,3])
        lines!(ax3,d.signal,label="signal")


        for i in 1:7
            d.eps = epsi[i]
            put_EOF_diff(d)
            calculate(d)
            lines!(ax3,vec(sum(d.RC,dims=2)).+elevate[i],label="eps $(d.eps)")
        end

        axislegend(ax3)

        save(dir*"test.png",F)
    end
end

