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

    BW1 = parse(Bool,filenamesplit[7])
    BW2 = parse(Bool,filenamesplit[8])
    method_sea = filenamesplit[9]
    Rep_W = parse(Int64,filenamesplit[10][1:end-5])

    B = (method_sea == "lSSA") ?  BW1   :  BW2
    Wf = (method_sea == "lSSA") ?  BW2   :  BW1
    P = N - W +1
    
    dictfile = load(datadir*filename)

    EOF = dictfile["EOF"]
    PC = dictfile["PC"]
    RC = dictfile["RC"]
    
    #RC = (B == true) ? hcat([cutoff(dictfile["RC"][:,i],365) for i=1:k]...) : RC

    it_num = 16

    sigma = Float32.([norm(PC[i,:]) for i=1:k])
    EOF_T = Float32.(2*W ./ [iterated_hilbert(EOF[:,i],it_num)[1] for i=1:k])
    PC_T = Float32.(2*P ./ [iterated_hilbert(PC[i,:],it_num)[1] for i=1:k])
    #RC_T = Float32.(2*N ./ [iterated_hilbert(RC[:,i],it_num)[1] for i=1:k])


    return name,W,i,method_sea,Rep_W,B,Wf,hcat(sigma,EOF_T,PC_T)
end

#axis: quick definition of xlabel, xlims, xticks, xscale, xflip, and xtickfont
#ticks: position of the ticks and the labels # with ticks=native you can let calculate
#marker:markershape,markersize,markeralpha,markercolor,markerstrokewidth,markerstrokealpha,markerstrokecolor,markerstrokestyle


datadir = "/net/scratch/lschulz/preprocessing/"
#filter the names of the ssa sweep
fullnames = readdir(datadir)

all_spectra= [fullspectrum(filename,datadir) for filename in fullnames] #   name,W,variable, hcat(sigma,EOF_T,PC_T,RC_T)

nameslist = readdir("/net/scratch/lschulz/cropdata_deseason/")

#find out which spots there were
names_ind       = ([findall(i->i[1] == W,all_spectra) for W in nameslist],nameslist)
deleteat!(names_ind[2],findall(isempty,names_ind[1]))
deleteat!(names_ind[1],findall(isempty,names_ind[1]))

#one plotting per spot!

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

function plotting(names_ind)

    names = names_ind[2]
    names_ind = names_ind[1]

    for (i,name) in enumerate(names)
        Dirname = dir*"spectrum_"*name
        ind = names_ind[i]

        #the indices for the calculations at this spot
        B_ind           = ([findall(i->i[6] == W,all_spectra[ind]) for W in [true,false]],[true,false])
        Wf_ind          = ([findall(i->i[7] == W,all_spectra[ind]) for W in [true,false]],[true,false])
        method_sea_ind  = ([findall(i->i[4] == W,all_spectra[ind]) for W in ["gSSA","lSSA"]],["gSSA","lSSA"])
        Rep_W_ind       = ([findall(i->i[5] == W,all_spectra[ind]) for W in [0,366,416,438,456]],[0,366,416,438,456])
        W_ind           = ([findall(i->i[2] == W,all_spectra[ind]) for W in [5478,5844,5946,]],[5478,5844,5946,])

        #finding the two different W full and the largest full year
        a=ind[1]
        b = ind[end]
        fullW(list) = deleteat!(Array(a:b),findall(isempty,[in(k,Wf_ind[1][1]) && in(k,list) ? k : [] for k = a:b]))
        yearW(list) = deleteat!(Array(a:b),findall(isempty,[in(k,Wf_ind[1][2]) && in(k,list) ? k : [] for k = a:b]))

        #5 different boundary techniques: none, gssa, lssa1:4
        #each with Wf harmonic or not!
        none = ("none",fullW(B_ind[1][2] ),yearW(B_ind[1][2] ))
        gSSA = ("gSSA",fullW(Rep_W_ind[1][1]),yearW(Rep_W_ind[1][1]))
        lSSA1= ("l 1a",fullW(Rep_W_ind[1][2]),yearW(Rep_W_ind[1][2]))
        lSSA2= ("l 7a/6",fullW(Rep_W_ind[1][3]),yearW(Rep_W_ind[1][3]))
        lSSA3= ("l 6a/5",fullW(Rep_W_ind[1][4]),yearW(Rep_W_ind[1][4]))
        lSSA4= ("l 5a/4",fullW(Rep_W_ind[1][5]),yearW(Rep_W_ind[1][5]))

        

        #one plot for the T/sigma EOF
        p2 = plot(xaxis = ("T[years]", Array(0:1:15)),yaxis=("sigma"),title="spectrum EOF")
        #one plot for the T/sigma PC
        p3 = plot(xaxis = ("T[years]", Array(0:1:15)),yaxis=("sigma"))

        #create different colors for the different methods, symbols for the Wfull,Wyear
        colors = [:black,:red,:blue,:yellow,:green,:orange]
        symbols = [:x,:o]
        labeladdendum = [" F"," year"]

        #here we create the lists for the individual plots
        labels = []
        p2X = []
        p2Y = []
        p3X = []
        p3Y = []

        for (j,B_method_ind) in enumerate([none,gSSA,lSSA1,lSSA2,lSSA3,lSSA4])
            color = colors[j]
            label = B_method_ind[1]
            fullW_list = B_method_ind[2]
            yearW_list = B_method_ind[3]

            for (l,W_method_ind) in enumerate([fullW_list,yearW_list])
                symbol = symbols[l]
                label = label * labeladdendum[l]

                list_sigma = [all_spectra[ind][W_method_ind][k][8][:,1] for k=1:length(W_method_ind)]
                list_EOF_T = [all_spectra[ind][W_method_ind][k][8][:,2] ./365 for k=1:length(W_method_ind)]
                list_PC_T = [all_spectra[ind][W_method_ind][k][8][:,3] ./365 for k=1:length(W_method_ind)]

                labels =    push!(labels,label)
                p2X =       push!(p2X,list_EOF_T)
                p2Y =       push!(p2Y,list_sigma)
                p3X =       push!(p3X,list_PC_T)
                p3Y =       push!(p3Y,list_sigma)

            end
        end

        savefig(scattering("eof",labels,p2X,p2Y,"EOF T[year]","sigma"),Dirname*"eof")
        savefig(scattering("pc",labels,p3X,p3Y,"PC T[year]","sigma"),Dirname*"pc")

    end
end
plotting(names_ind)
"""
for the spots
    we have different plots

        for the boundary_methods
            we have different labels and colors
            put a label in the list
            for the W_choice
                we have different symbols
                put the symbol in the list
                put the X, Y in the list

    return list_labels, listX,listY

plotting works by list_labels list_x list_y



k = 32

p1a = scatter(1:32,list_EOF_T,title="EOF T ",label=label,xlabel="mode",ylabel="T[years]",markershape=symbol,c=color)
p1b = scatter(1:32,list_PC_T,title="PC  T ",label=label,xlabel="mode",ylabel="T[years]",markershape=symbol,c=color)
p1c = scatter(1:32,list_sigma,title="sigma ",label=label,xlabel="mode",ylabel="sigma",markershape=symbol,c=color)

p1 = plot!(p1a,p1b,p1c)

p2 = scatter!(list_EOF_T,list_sigma,label=label,markershape=symbol,c=color)
p3 = scatter!(list_PC_T,list_sigma,label=label,markershape=symbol,c=color)

savefig(p1,Dirname*"_1")
savefig(p2,Dirname*"_2")
savefig(p3,Dirname*"_3")

"""