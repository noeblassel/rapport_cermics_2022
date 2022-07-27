using Plots, LinearAlgebra

norton_dat=readlines("norton_data.txt")
nemd_dat=readlines("nemd_data.txt")

norton_dat=map(l->split(l," "), norton_dat[2:end])
nemd_dat=map(l->split(l, " "),nemd_dat[2:end])

norton_color=reduce(hcat,[parse.(Float64,l[2:end]) for l in filter(l->l[1]=="COLOR",norton_dat)])
norton_single=reduce(hcat,[parse.(Float64,l[2:end]) for l in filter(l->l[1]=="SINGLE",norton_dat)])

nemd_color=reduce(hcat,[parse.(Float64,l[2:end]) for l in filter(l->l[1]=="COLOR",nemd_dat)])
nemd_single=reduce(hcat,[parse.(Float64,l[2:end]) for l in filter(l->l[1]=="SINGLE",nemd_dat)])

λ_color=norton_color[1,:]
λ_single=norton_single[1,:]

η_color=nemd_color[1,:]
η_single=nemd_single[1,:]

R_norton_color=norton_color[2,:]
R_norton_single=norton_single[2,:]

R_nemd_color=nemd_color[2,:]
R_nemd_single=nemd_single[2,:]

std_norton_color= sqrt.( norton_color[3,:] ./ norton_color[5,:])
std_norton_single= sqrt.( norton_single[3,:] ./ norton_single[5,:])

std_nemd_color= sqrt.( nemd_color[3,:] ./ nemd_color[5,:])
std_nemd_single= sqrt.( nemd_single[3,:] ./ nemd_single[5,:])

norton_color_perm=sortperm(λ_color)
λ_color=λ_color[norton_color_perm]
R_norton_color=R_norton_color[norton_color_perm]
std_norton_color=std_norton_color[norton_color_perm]

norton_single_perm=sortperm(λ_single)
λ_single=λ_single[norton_single_perm]
R_norton_single=R_norton_single[norton_single_perm]
std_norton_single=std_norton_single[norton_single_perm]

nemd_color_perm=sortperm(η_color)
η_color=η_color[nemd_color_perm]
R_nemd_color=R_nemd_color[nemd_color_perm]
std_nemd_color=std_nemd_color[nemd_color_perm]

nemd_single_perm=sortperm(η_single)
η_single=η_single[nemd_single_perm]
R_nemd_single=R_nemd_single[nemd_single_perm]
std_nemd_single=std_nemd_single[nemd_single_perm]

n_regr=10
n_plot=12
n_nl_plot=30    

ρ_norton_color=inv(dot(λ_color[1:n_regr],λ_color[1:n_regr]))*dot(λ_color[1:n_regr],R_norton_color[1:n_regr])
ρ_norton_single=inv(dot(λ_single[1:n_regr],λ_single[1:n_regr]))*dot(λ_single[1:n_regr],R_norton_single[1:n_regr])


x_lims=(0,max(λ_color[n_plot],λ_single[n_plot]))

norton_plot=plot(legend=:topleft,xlabel="⟨λ⟩",ylabel="Mobility response",xlims=xlims)
full_plot=plot(legend=:topleft,xlabel="Forcing",ylabel="Mobility response")

scatter!(norton_plot,λ_color[1:n_plot],R_norton_color[1:n_plot],markershape=:xcross,color=:blue,label="color drift (ρ≈$(round(ρ_norton_color,digits=3)))",xerr=std_norton_color)
scatter!(norton_plot,λ_single[1:n_plot],R_norton_single[1:n_plot],markershape=:xcross,color=:red,label="single drift (ρ≈$(round(ρ_norton_single,digits=3)))",xerr=std_norton_single)
plot!(norton_plot,t->ρ_norton_color*t,label="",linestyle=:dot,color=:blue)
plot!(norton_plot,t->ρ_norton_single*t,label="",linestyle=:dot,color=:red)

scatter!(full_plot,λ_color,R_norton_color,markershape=:xcross,label="norton color")
scatter!(full_plot,λ_single,R_norton_single,markershape=:xcross,label="norton single")
scatter!(full_plot,η_color,R_nemd_color,markershape=:xcross,label="thevenin color")
scatter!(full_plot,η_single,R_nemd_single,markershape=:xcross,label="thevenin single")

savefig(norton_plot,"norton_mobility_plot.pdf")
savefig(full_plot,"norton_mobility_full.pdf")


vars_color=norton_color[3,1:n_regr] ./ norton_color[5,1:n_regr]
vars_single=norton_single[3,1:n_regr] ./ norton_single[5,1:n_regr]

Rλ_color=dot(λ_color[1:n_regr],R_norton_color[1:n_regr])
λ2_color=dot(λ_color[1:n_regr],λ_color[1:n_regr])

Rλ_single=dot(λ_single[1:n_regr],R_norton_single[1:n_regr])
λ2_single=dot(λ_single[1:n_regr],λ_single[1:n_regr])

est_var_color=dot(vars_color, ((R_norton_color[1:n_regr]*λ2_color-2Rλ_color*λ_color[1:n_regr]) / λ2_color^2) .^ 2)
est_var_single=dot(vars_single, ((R_norton_single[1:n_regr]*λ2_single-2Rλ_single*λ_single[1:n_regr]) / λ2_single^2) .^ 2)


α=1/999
println("ρ_norton_color: $ρ_norton_color , std $(sqrt(est_var_color)), mobility: $(inv(1+α)*(ρ_norton_color+α))")
println("ρ_norton_single : $ρ_norton_single , std $(sqrt(est_var_single)), mobility: ")