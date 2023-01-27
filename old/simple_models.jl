include("/net/home/lschulz/scripts/toolbox.jl")
# we have three different models
simple = npzread(dir*"simple")


#with each 4 different scenarios of w inside them simple[i,:] of 8000 length
# w1 = 10, w2 = 60, w3 = 300, w4 = 800

#global diffusion returns 8000x20
#parameter variation: W : 8 x s.t. 8x8000x20
glo_diff=zeros(4,8,8000,20)
for index = 1:4
    for (i,W) in enumerate(Int.(100:100:800))
        glo_diff[index,i,:,:] = global_diffmap(centralizer(simple[index,:]),W)
    end
end
npzwrite(dir*"glo_diff_simple",glo_diff)



overlay = npzread(dir*"overlay")

glo_diff=zeros(4,8,8000,20)
for index = 1:4
    for (i,W) in enumerate(Int.(100:100:800))
        glo_diff[index,i,:,:] = global_diffmap(centralizer(overlay[index,:]),W)
    end
end
npzwrite(dir*"glo_diff_overlay",glo_diff)




noise = npzread(dir*"noise")

glo_diff=zeros(4,8,8000,20)
for index = 1:4
    for (i,W) in enumerate(Int.(100:100:800))
        glo_diff[index,i,:,:] = global_diffmap(centralizer(noise[index,:]),W)
    end
end
npzwrite(dir*"glo_diff_noise",glo_diff)

