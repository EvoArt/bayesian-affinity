using Pkg
Pkg.activate(".")
using Turing, Random, DelimitedFiles, CSV, DataFrames,Phylo
import ADTypes, Mooncake
include("utils.jl")
PA_df = CSV.read("global_ants.csv",DataFrame) 
taxa = PA_df.Column1
gens = [split(taxon,".")[1] for taxon in taxa]
PA= transpose(Array(PA_df[:,2:end]))

tree = open(parsenewick, "15K_NCuniform_crown_mcc.tre")
in_tree = hasnode.(tree,taxa)
taxa = taxa[in_tree]
PA = PA[:,in_tree]
gens = gens[in_tree]
ugens = unique(gens)

nodes = collect(nodenameiter(tree))
gen_nodes = [mrca(tree,filter(contains(gen),nodes)) for gen in ugens]
n_gen = length(gen_nodes)
n_samp =size(PA,1)

phydist = zeros(n_gen,n_gen)
for j in 2:n_gen
    for i in 1:j-1
        d =distance(tree,gen_nodes[i],gen_nodes[j])
        phydist[i,j] = d
        phydist[j,i] = d
    end
end

GPA = zeros(Int,n_samp,n_gen)
for (i,gen) in enumerate(ugens)
    GPA[:,i] = Bool.(maximum.(eachrow(PA[:,gens .== gen])))
end
n_sample, n_obs = size(GPA)

mask = (size(GPA,1)-0) .> sum.(eachcol(GPA)) .> 0
GPA = GPA[:,mask]
phydist = scale_center(phydist[mask,mask])
n_pars = sum([sum(gens .== gen) for gen in ugens])+1
df_name = "ants_genus_chn_$(ARGS[1]).csv" 

k = zeros(Int64, n_obs,n_obs)
for i in 2:n_obs
    for j in 1:i-1
        k[i,j] = sum(GPA[:,i] .* GPA[:,j])
    end
end

pars = phydist,k, n_obs,n_sample, vec(sum(GPA, dims = 1))
@model function genus_model(X,k, n_obs,N,m)
    a ~ Normal(0.0,4.0)
    λ ~ filldist(truncated(Normal(0.0,1.0)),n_obs)
    β ~ truncated(Normal(0.0,1.0))
    lp = 0.0
    for j in 1:n_obs-1
        mB = m[j]
        for i in j+1:n_obs
            mA = m[i]
            μ = a + β * X[i,j] + λ[i] + λ[j]
            lp +=  FNCH_logpdf(mA,N-mA,mB,exp(μ),k[i,j])
        end
    end
    Turing.@addlogprob! lp
end


m = genus_model(pars...)
n_genera = size(phydist,1)
chn = sample(m,NUTS(200,0.85; adtype=ADTypes.AutoMooncake(config=nothing)),1000,initial_params = randn(n_genera+2));
CSV.write(df_name,DataFrame(chn))
