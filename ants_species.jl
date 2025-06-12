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
phydist = CSV.File("ant_dist.csv",header = false) |> Tables.matrix
n_sample, n_obs = size(PA)

mask = (size(PA,1)-0) .> sum.(eachcol(PA)) .> 0
PA = PA[:,mask]
phydist = phydist[mask,mask]
gens = gens[mask]
ugens = unique(gens)

c = sampler_ind = 3
c_name = c+1
counts = [sum(gens .== gen) for gen in ugens]
mask = Inf .> counts .>c
ugens = ugens[mask]
n_pars = sum([sum(gens .== gen) for gen in ugens])+1

P= []
for gen in ugens
    inds = gens .== gen
    pa = PA[:,inds]
    dist = scale_center(phydist[inds,inds])
    d = dist[dist .!=0]
    if length(unique(d)) > 0.25* length(d) # filter on number of unique distances
        n_sample, n_obs = size(pa)
        k = zeros(Int64, n_obs,n_obs)
        for i in 2:n_obs
            for j in 1:i-1
                k[i,j] = sum(pa[:,i] .* pa[:,j])
            end
        end
        params = dist,k, n_obs,n_sample, vec(sum(pa, dims = 1))
        push!(P,params)
    end
end
P = Tuple.(P)

function sub_model(X,k, n_obs,N,m,β,λ,a)
    lp = 0.0
    for j in 1:n_obs-1
        mB = m[j]
        for i in j+1:n_obs
            mA = m[i]
            μ = a + β*X[i,j] + λ[i] +λ[j]
            lp +=  FNCH_logpdf(mA,N-mA,mB,exp(μ),k[i,j])
        end
    end
    lp
end

@model function hier(P,n_tot,n = length(P))
    a ~ Normal(0.0,4.0)
    β ~ filldist(Normal(0.0,1.0),n)
    mu ~ Normal(0.0,1.0)
    sig ~ Gamma(2.0,0.1)
    λ ~ filldist(Normal(0,1.0),n_tot)

    λ_offset =0
    for group in 1:n    
        Turing.@addlogprob!  sub_model(P[group]...,(β[group]*sig)+mu,λ[1+λ_offset:P[group][3]+λ_offset],a)
        λ_offset += P[group][3]
    end
end

@model function pool(P,n_tot,n = length(P))
    a ~ Normal(0.0,4.0)
    β ~ Normal(0.0,1.0)
    λ ~ filldist(Normal(0.0,1.0),n_tot)

    λ_offset =0
    for group in 1:n    
        Turing.@addlogprob!  sub_model(P[group]...,β,λ[1+λ_offset:P[group][3]+λ_offset],a)
        λ_offset += P[group][3]
    end
end


n_species = sum([p[3] for p in P])
β = randn(length(P))
a = 4 .*randn(1)
mu = randn()
sig = rand(truncated(Gamma(2.0,0.1),lower = 0, upper = 3))
λ = randn(n_species)
hier_ip = vcat(a...,β...,mu,sig,λ...)
pool_ip = vcat(a...,mu,λ...)

pool_mod = pool(P,n_species)
pool_chn = sample(pool_mod,NUTS(200,0.85; adtype=ADTypes.AutoMooncake(config=nothing)),1000,initial_params=pool_ip);
CSV.write("results/ants_species_pool_chn_$(ARGS[1]).csv",DataFrame(pool_chn))

hier_mod = hier(P,n_species)
hier_chn = sample(hier_mod,NUTS(200,0.85; adtype=ADTypes.AutoMooncake(config=nothing)),1000,initial_params=hier_ip);
CSV.write("results/ants_species_hier_chn_$(ARGS[1]).csv",DataFrame(hier_chn))

