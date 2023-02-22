"""
single signals for the presentation
    """

F = Figure(resolution=(800,800))
ax = Axis(F[1,1],
yticksvisible = false,yticklabelsvisible = false,xlabel="t (d)",ylabel="GPP")

for (i,spot) = enumerate(1:20)
#for (i,spot) in enumerate([2,4,8,10,13,14,15,17])
Filename = create_file_list(outdir,method,W,vari,preproc)[spot]
file = load(Filename)
signal = file["signal"] .+ i*6

series!(ax,1:N,signal',solid_color="black",linewidth=1)

end

save(savedir*"spots.png",F)