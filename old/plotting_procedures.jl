"""
THE RED THREAD
results of the form
one file per spot, per variable, per W

results that we want to plot
sigma / T scatterplot for a different color per spot/variable

for this we need
reading of a file with EOF PC RC gives its T and sigma per k, so this returns a matrix with EOF,PC,RC T for given iteration number
magic axis
scattering function (name:spot/variable,X,Y,axes) that has good coloring scheme
"""

#this function does not cut off the boundary but could easily
#this function uses a hardcoded variable number of hilbert transform iterations
#   name,W,variable, hcat(sigma,EOF_T,PC_T,RC_T)
function fullspectrum(filename,dir)
    filenamesplit = split(filename,"_")

    method = filenamesplit[1]
    i = parse(Int64,filenamesplit[2])
    name = filenamesplit[3]
    N = parse(Int64,filenamesplit[4])
    W = parse(Int64,filenamesplit[5])
    k = parse(Int64,filenamesplit[6])

    a = parse(Float64,filenamesplit[7])
    e = parse(Float64,filenamesplit[8])
    t = parse(Float64,filenamesplit[9])

    variable = filenamesplit[10]

    P = N - W +1
    
    dictfile = load(datadir*filename)

    EOF = dictfile["EOF"]
    PC = dictfile["PC"]
    RC = dictfile["RC"]#hcat([cutoff(dictfile["RC"][:,i],365) for i=1:k]...)

    it_num = 4

    sigma = Float32.([norm(PC[i,:]) for i=1:k])
    EOF_T = Float32.(2*W ./ [iterated_hilbert(EOF[:,i],it_num)[1] for i=1:k])
    PC_T = Float32.(2*P ./ [iterated_hilbert(PC[i,:],it_num)[1] for i=1:k])
    RC_T = Float32.(2*N ./ [iterated_hilbert(RC[:,i],it_num)[1] for i=1:k])


    return name,W,variable, hcat(sigma,EOF_T,PC_T,RC_T)
end

#axis: quick definition of xlabel, xlims, xticks, xscale, xflip, and xtickfont
#ticks: position of the ticks and the labels # with ticks=native you can let calculate
#marker:markershape,markersize,markeralpha,markercolor,markerstrokewidth,markerstrokealpha,markerstrokecolor,markerstrokestyle

# to make the days into years we need to scale the lists prior!
yearxaxis = ("T[years]", :none, Array(0:182.5:365*7),1,false,4)

function scattering(title,list_labels,list_X,list_Y,xaxis,yaxis)

    #this is to plut very many different spots on one plot  without labels
    cc = cgrad(:hawaii, length(list_labels), categorical = true)

    p = plot(xaxis=xaxis,yaxis=yaxis,dpi=800,title=title)

    marker = (:+,5)

    for (i,lab) = enumerate(list_labels)
        p = scatter!(list_X[i],list_Y[i],label=lab,marker=marker,c = cc[i])
    end
    return p
end

datadir = "/net/scratch/lschulz/ssa_sweep_results/"
#filter the names of the ssa sweep
fullnames = readdir(datadir)

Wl = 10
Wlist = Int.(floor.(Array(365.25:365.25:365.25*10)))

all_spectra= [fullspectrum(filename,datadir) for filename in fullnames] #   name,W,variable, hcat(sigma,EOF_T,PC_T,RC_T)
Wsort = [findall(i->i[2] == W,all_spectra) for W in Wlist] # this gives the list of indices for each W

for (i,W) in enumerate(Wlist)
    indices = Wsort[i]
    Dirname = dir*"spectrum_$W"
    title = "SSA for W = $W"
    list_labels = [all_spectra[indices][k][1] for k=1:length(indices)]
    list_X = [all_spectra[indices][k][4][:,4] ./365 for k=1:length(indices)]
    list_Y = [all_spectra[indices][k][4][:,1] for k=1:length(indices)]
    xaxis = ("T[years]", Array(0:1:15))
    yaxis=("sigma")
    p = scattering(title,list_labels,list_X,list_Y,xaxis,yaxis)
    p = plot!(legend=:none)
    savefig(p,Dirname)
end
"""
Wl = 10
Wlist = Int.(floor.(Array(365.25:365.25:365.25*10)))
k=32
spots = 70

M_sigma= Array{Float32}(undef,spots,Wl,k)
M_EOF_T1= Array{Float32}(undef,spots,Wl,k)
M_EOF_T4= Array{Float32}(undef,spots,Wl,k)
M_EOF_T8= Array{Float32}(undef,spots,Wl,k)
M_PC_T1= Array{Float32}(undef,spots,Wl,k)
M_PC_T4= Array{Float32}(undef,spots,Wl,k)
M_PC_T8= Array{Float32}(undef,spots,Wl,k)
M_RC_T1= Array{Float32}(undef,spots,Wl,k)
M_RC_T4= Array{Float32}(undef,spots,Wl,k)
M_RC_T8= Array{Float32}(undef,spots,Wl,k)


#loop all the names find out W and ind and put it in the matrix
for i in fullnames
    filename = i
    filenamesplit = split(filename,"_")
    method = filenamesplit[1]
    ind = parse(Int64,filenamesplit[2])
    name = filenamesplit[3]
    N = parse(Int64,filenamesplit[4])
    W = parse(Int64,filenamesplit[5])
    k = parse(Int64,filenamesplit[6])

    dictfile = load(datadir*filename)

    sigma = dictfile["sigma"]
    EOF_T1 = dictfile["EOF_T1"]
    EOF_T4 = dictfile["EOF_T4"]
    EOF_T8 = dictfile["EOF_T8"]
    PC_T1 = dictfile["PC_T1"]
    PC_T4 = dictfile["PC_T4"]
    PC_T8 = dictfile["PC_T8"]
    RC_T1 = dictfile["RC_T1"]
    RC_T4 = dictfile["RC_T4"]
    RC_T8 = dictfile["RC_T8"]

    Wi = findfirst(x->x==W,Wlist)
    for ki in 1:k
        M_sigma[ind,Wi,ki] = sigma[ki]
        M_EOF_T1[ind,Wi,ki] = EOF_T1[ki]
        M_EOF_T4[ind,Wi,ki] = EOF_T4[ki]
        M_EOF_T8[ind,Wi,ki] = EOF_T8[ki]
        M_PC_T1[ind,Wi,ki] = PC_T1[ki]
        M_PC_T4[ind,Wi,ki] = PC_T4[ki]
        M_PC_T8[ind,Wi,ki] = PC_T8[ki]
        M_RC_T1[ind,Wi,ki] = RC_T1[ki]
        M_RC_T4[ind,Wi,ki] = RC_T4[ki]
        M_RC_T8[ind,Wi,ki] = RC_T8[ki]
    end
end

#plotting best utilizes one list for x and one for y
#first plot is the list is W

function plotting_results_lists(matrix,Wl,spots,k,Wlist_short)
    Wlist = []
    var_list = []
    for ind_s in 1:spots
        for ind_W in 1:Wl
            for ind_k in 1:k
                Wlist =push!(Wlist,Wlist_short[ind_W])
                var_list = push!(var_list,matrix[ind_s,ind_W,ind_k])
            end
        end
    end
    return Wlist,var_list
end
   
macro Name(arg)
    string(arg)
end

function scatter_W_matrix(matrix,name)
    W_list,var_list = plotting_results_lists(matrix,Wl,spots,k,Wlist)
    p = scatter(W_list,var_list,legend=:none,xlabel="W[days]",ylabel=name,
    markersize=2,c=:black,marker=:x)

    savefig(p,dir*"sweep_$name")

end
Mnames = ["sigma","EOF_T1","EOF_T4","EOF_T8","PC_T1","PC_T4","PC_T8","RC_T1","RC_T4","RC_T8"]

for (i,M) in enumerate([M_sigma,M_EOF_T1,M_EOF_T4,M_EOF_T8,M_PC_T1,M_PC_T4,M_PC_T8,M_RC_T1,M_RC_T4,M_RC_T8])
    scatter_W_matrix(M,Mnames[i])
end

# this brute coloring scheme can be improved
# with findall and then put the integer numbers inside
# for (index,W) in enumerate(Wl) x[findall(w->w=W,x)] = index end
x,RCT = plotting_results_lists(M_RC_T8,Wl,spots,k,Wlist)
x,sigmal = plotting_results_lists(M_sigma,Wl,spots,k,Wlist)
cc = palette([:purple, :green], 4000)
savefig(
           scatter(RCT,sigmal,legend=:none,markersize=2,marker=:X,c=cc[x],
       ylabel="sigma",dpi=1000,alpha=0.3,xaxis=("T",Wlist)),
       dir*"test")

#what about the entropy for T and sigma

function entropy(signal::Vector{Float32})
    signal = signal ./ sum(signal)
    return - sum(signal .* log.(abs.(signal)))
end

function entropy_list(matrix,Wl,spots,k,Wlist_short)
    Wlist = []
    var_list = []
    for ind_s in 1:spots
        for ind_W in 1:Wl
                Wlist =push!(Wlist,Wlist_short[ind_W])
                var_list = push!(var_list,entropy(matrix[ind_s,ind_W,:]))

        end
    end
    return Wlist,var_list
end

function scatter_W_entropy(matrix,name)
    W_list,var_list = entropy_list(matrix,Wl,spots,k,Wlist)
    p = scatter(W_list,var_list,legend=:none,xlabel="W[days]",ylabel="$(name)_entropy",
    markersize=2,c=:black,marker=:x)

    savefig(p,dir*"sweep_$(name)_entropy")

end

for (i,M) in enumerate([M_sigma,M_EOF_T1,M_EOF_T4,M_EOF_T8,M_PC_T1,M_PC_T4,M_PC_T8,M_RC_T1,M_RC_T4,M_RC_T8])
    scatter_W_entropy(M,Mnames[i])
end

"""