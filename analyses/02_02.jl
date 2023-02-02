#plot all signals to find out where the amplitude are least modulated

include("/net/home/lschulz/reduce-daemensionality/toolbox.jl")

savedirname=dir*"signals/"
outdir="/net/scratch/lschulz/fluxfullset_1a/"
yearsamples=365
preproc="raw."
W = 365

# aparently there is only variable number one computed which is GPP

for method = ["ssa"], vari = 1, spot = 1:20

    filename = create_file_list(outdir,method,W,vari,preproc)[spot]

    file = load(filename)
    eof = file["EOF"]
    lambda = file["lambda"]
    pc = file["PC"]
    rc = file["RC"]
    signal = file["signal"]
    N = length(signal)
    F,ax,s = lines(1:N,signal,solid_color="black")
    save(savedirname*"$(method)_$(spot)_$(vari)_$W.png",F)
end

# 3,9,14,15,16 seem to be suitable

#now, lets do HT iterations of signcount on the first 10 modes which need to be sorted
#do print to find out how many modes are somehow somewhat stable with embedding of 1a (with 7a none were)


function protophase(signal::Vector{Float32})
    ht = hilbert(Float64.(signal))
    return Float32.(atan.(imag(ht), real(ht)))
end

count_sign_flip(signal) = sum([(sign.(signal[i])!=sign.(signal[i+1])) ? 1 : 0 for i in 1:length(signal)-1])

vari = 1
n_modes = 6
Nit = 4

for spot = [3,9,14,15,16]
    for method = ["ssa","diff"]

        filename = create_file_list(outdir,method,W,vari,preproc)[spot]

        file = load(filename)
        eof = file["EOF"]
        lambda = file["lambda"]

        indices = sortperm(lambda,rev=true)[1:n_modes]

        print(spot,"\t",method,"\n")

        for kap in 1:n_modes
            print(indices[kap],"\t")
            signal = eof[:,indices[kap]]
            for i=1:Nit
                phi = protophase(signal)
                T = count_sign_flip(phi)
                print("$T \t")
                signal = phi
            end
          
vari = 1
n_modes = 6
Nit = 4

for spot = [3,9,14,15,16]
    for method = ["ssa","diff"]

        filename = create_file_list(outdir,method,W,vari,preproc)[spot]

        file = load(filename)
        eof = file["EOF"]
        lambda = file["lambda"]

        indices = sortperm(lambda,rev=true)[1:n_modes]

        print("#\t",spot,"\t",method,"-----------------\n")

        for kap in 1:n_modes
            print("#\t",indices[kap],"|\t")
            signal = eof[:,indices[kap]]
            for i=1:Nit
                phi = protophase(signal)
                T = count_sign_flip(phi)
                print("$T \t")
                signal = phi
            end
            print("\n")
        end
    end
end

"""
results
"""

#       3       ssa-----------------
#       1|      2       2       2       4 
#       2|      2       2       2       4 
#       3|      4       4       4       4 
#       4|      4       4       4       4 
#       5|      6       8       6       6 
#       6|      6       12      10      12 
#       3       diff-----------------
#       3|      2       2       2       4 
#       2|      2       2       4       6 
#       4|      4       4       4       4 
#       5|      4       4       4       4 
#       6|      6       8       12      14 
#       7|      6       8       8       10 
#       9       ssa-----------------
#       1|      2       2       2       4 
#       2|      2       2       2       4 
#       3|      4       4       4       6 
#       4|      4       6       4       4 
#       5|      6       8       10      9 
#       6|      6       10      8       8 
#       9       diff-----------------
#       12|     12      23      36      38 
#       22|     18      22      26      32 
#       17|     16      28      22      34 
#       5|      8       10      12      10 
#       1|      7       13      13      19 
#       18|     20      34      28      36 
#       14      ssa-----------------
#       1|      2       2       2       4 
#       2|      2       2       2       4 
#       3|      4       4       4       6 
#       4|      4       4       4       4 
#       5|      8       8       8       10 
#       6|      8       10      12      18 
#       14      diff-----------------
#       3|      2       2       2       2 
#       2|      2       2       2       4 
#       4|      4       4       4       4 
#       1|      2       2       2       2 
#       5|      4       6       4       6 
#       6|      6       6       6       6 
#       15      ssa-----------------
#       1|      2       2       2       4 
#       2|      2       2       2       4 
#       3|      6       6       6       8 
#       4|      6       8       8       10 
#       5|      6       8       6       8 
#       6|      4       6       4       4 
#       15      diff-----------------
#       5|      4       2       6       12 
#       3|      12      28      32      41 
#       2|      2       4       4       4 
#       4|      4       2       2       2 
#       10|     10      10      10      12 
#       1|      6       12      8       12 
#       16      ssa-----------------
#       1|      2       2       2       4 
#       2|      2       2       2       4 
#       3|      4       4       4       4 
#       4|      4       4       4       4 
#       5|      6       6       6       6 
#       6|      6       6       6       6 
#       16      diff-----------------
#       2|      2       2       2       4 
#       1|      2       4       2       2 
#       3|      2       2       2       2 
#       4|      4       4       4       6 
#       6|      6       10      6       7 
#       5|      4       4       4       4 


#filter that finds a modes that is only 2,4,6,8 among the first 3 iterations, counting only the first 8 modes

function hard_mode_filter(one,two,three)

    if one==two==three==2
        return true,2
    elseif one==two==three==3
        return true,3
    elseif one==two==three==4
        return true,4
    elseif one==two==three==6
        return true,6
    elseif one==two==three==8
        return true,8
    else
        return false,0
    end
end


vari = 1
n_modes = 6
Nit = 3

for spot = [3,9,14,15,16]
    for method = ["ssa","diff"]

        filename = create_file_list(outdir,method,W,vari,preproc)[spot]

        file = load(filename)
        eof = file["EOF"]
        lambda = file["lambda"]

        indices = sortperm(lambda,rev=true)[1:n_modes]

        print("#\t",spot,"\t",method,"-----------------\n")

        for kap in 1:n_modes
            print("#\t",indices[kap],"|\t")
            signal = eof[:,indices[kap]]
            L = Int64[]
            for i=1:Nit
                phi = protophase(signal)
                T = count_sign_flip(phi)
                L = append!(L,Int64(T))
                print("$T \t")
                signal = phi
            end
            fil = hard_mode_filter(L...)
            if fil[1]
                print("HERE")
            end
            print("\n")
        end
    end
end



"""
results look promising
"""
#       3       ssa-----------------
#       1|      2       2       2       HERE
#       2|      2       2       2       HERE
#       3|      4       4       4       HERE
#       4|      4       4       4       HERE
#       5|      6       8       6 
#       6|      6       12      10 
#       3       diff-----------------
#       3|      2       2       2       HERE
#       2|      2       2       4 
#       4|      4       4       4       HERE
#       5|      4       4       4       HERE
#       6|      6       8       12 
#       7|      6       8       8 
#       9       ssa-----------------
#       1|      2       2       2       HERE
#       2|      2       2       2       HERE
#       3|      4       4       4       HERE
#       4|      4       6       4 
#       5|      6       8       10 
#       6|      6       10      8 
#       9       diff-----------------
#       12|     12      23      36 
#       22|     18      22      26 
#       17|     16      28      22 
#       5|      8       10      12 
#       1|      7       13      13 
#       18|     20      34      28 
#       14      ssa-----------------
#       1|      2       2       2       HERE
#       2|      2       2       2       HERE
#       3|      4       4       4       HERE
#       4|      4       4       4       HERE
#       5|      8       8       8       HERE
#       6|      8       10      12 
#       14      diff-----------------
#       3|      2       2       2       HERE
#       2|      2       2       2       HERE
#       4|      4       4       4       HERE
#       1|      2       2       2       HERE
#       5|      4       6       4 
#       6|      6       6       6       HERE
#       15      ssa-----------------
#       1|      2       2       2       HERE
#       2|      2       2       2       HERE
#       3|      6       6       6       HERE
#       4|      6       8       8 
#       5|      6       8       6 
#       6|      4       6       4 
#       15      diff-----------------
#       5|      4       2       6 
#       3|      12      28      32 
#       2|      2       4       4 
#       4|      4       2       2 
#       10|     10      10      10 
#       1|      6       12      8 
#       16      ssa-----------------
#       1|      2       2       2       HERE
#       2|      2       2       2       HERE
#       3|      4       4       4       HERE
#       4|      4       4       4       HERE
#       5|      6       6       6       HERE
#       6|      6       6       6       HERE
#       16      diff-----------------
#       2|      2       2       2       HERE
#       1|      2       4       2 
#       3|      2       2       2       HERE
#       4|      4       4       4       HERE
#       6|      6       10      6 
#       5|      4       4       4       HERE

"""
lets do the reconstruction using them against rest
"""



savedirname=dir*"hard_filter/"
vari = 1
n_modes = 8
Nit = 3

for spot = [3,9,14,15,16]
    for method = ["ssa","diff"]

        filename = create_file_list(outdir,method,W,vari,preproc)[spot]

        file = load(filename)
        eof = file["EOF"]
        lambda = file["lambda"]

        indices = sortperm(lambda,rev=true)[1:n_modes]
        L_modes = Int64[]

        for kap in 1:n_modes

            signal = eof[:,indices[kap]]
            L = Int64[]
            for i=1:Nit
                phi = protophase(signal)
                T = count_sign_flip(phi)
                L = append!(L,Int64(T))

                signal = phi
            end
            fil = hard_mode_filter(L...)
            if fil[1]
                L_modes = append!(L_modes,indices[kap])
            end
        end

        #plot the results
        signal = file["signal"]
        rc = file["RC"]
        trend = 0 .+ sum(rc[:,L_modes],dims=2)[:]
        residual = signal .- trend
        combined = hcat(signal,trend .+3 ,residual .+ 6)
        F,ax,s = series(1:length(signal),combined')
        save(savedirname*"_$(method)_$(spot)_trend.png",F)
    end
end
