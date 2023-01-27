"""
first step : make a 3d blob
"""

gaussian(x,sigma) = exp(- x^2 / sigma)

p = 1000
#xs = cos.(1:0.5:20)
#ys = sin.(1:0.5:20)
xs = gaussian.(rand(p),1)
ys = gaussian.(rand(p),1)
zs = gaussian.(rand(p),1/5)


f,ax,p = meshscatter(xs, ys, zs, markersize = 0.01, color = zs,limits=(-1,1))
save(dir*"test.png",f)
