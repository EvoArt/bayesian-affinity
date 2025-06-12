using Pkg
Pkg.activate(".")
using Phylo,CSV,DataFrames,DelimitedFiles

PA_df = CSV.read("global_ants.csv",DataFrame)
taxa = PA_df.Column1
PA= transpose(Array(PA_df[:,2:end]))
tree = open(parsenewick, "15k_NCuniform_crown_mcc.tre")
in_tree = hasnode.(tree,taxa)
taxa = taxa[in_tree]
PA = PA[:,in_tree]
n_sample, n_obs = size(PA)
phydist = zeros(n_obs,n_obs)
for j in 1:n_obs-1
    for i in j+1:n_obs
        d =distance(tree,taxa[i],taxa[j])
        phydist[i,j] = d
    end
end
writedlm("ant_dist.csv",phydist)