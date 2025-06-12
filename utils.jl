using Statistics, Memoization, Turing

# Data Transformations and "getters"

function scale_center(x)
    y = x/std(x[x .!=0])
    y[y .!=0] .-= mean(y[y .!=0])
    y
end

function condense(X,cols,func = any)
    n,m = size(X)
    Y = reduce(hcat,[func(X[:,i:i+cols-1],dims=2) for i in 1:cols:m-cols+1])
    Y
end

s(x,y,R) = x == y ? 0.0 : 1.0 
s(x::Real,y::Real,R) = abs(x-y)/R

function gower(x...; weights = ones(length(x[1]),length(x[1]),length(x)))
    n = length(x[1])
    m = length(x)
    R = ones(length(x))
    
    for i in 1:m
        if all(isa.(x[i],Number))
            R[i] = maximum(x[i]) - minimum(x[i])
        end
    end

    S = zeros(n,n)
    for j in 1:n-1
        for i in j+1:n
            dist = sum(weights[i,j,:] .* [s(x[k][i],x[k][j],R[k]) for k in 1:m] )/sum(weights[i,j,:])
            S[i,j] = dist
            S[j,i] = dist
        end
    end
    S
end

gower(X::Array; weights = ones(size(X,1),size(X,1),size(X,2))) = gower(eachcol(X)...,weights = weights)

check_seq(x) = (typeof(x) == BitVector) || (typeof(x) == Vector{Bool}) || (typeof(x) == Vector{Union{Missing, Bool}})
check_sum(x) = !check_seq(x) || (0 < sum(x) < length(x))
remove_empties(df) = df[:,check_sum.(eachcol(df))]

check_sum2(x) = !check_seq(x) || (1 < sum(x) < (length(x))-1)
remove_nearly_empties(df) = df[:,check_sum2.(eachcol(df))]

function remove_empties!(df)
    df = remove_empties(df)
end

get_part(x,i) = [split(name,"_genus_")[i] for name in x]
get_gens(df) = get_part(DataFrames.names(df)[check_seq.(eachcol(df))],2)
get_seqs(df) = get_part(DataFrames.names(df)[check_seq.(eachcol(df))],1)

# Working with Biom files

function obs_mat(biom)
    obs_ids = biom["observation/ids"]
    samp_ids = biom["sample/ids"]
    obs_data = biom["observation/matrix/data"]
    obs_inds = 1 .+ biom["observation/matrix/indices"]
    samp_inds = 1 .+biom["sample/matrix/indices"]
    n_samp,n_obs = length(samp_ids),length(obs_ids)
    mat = zeros(n_samp,n_obs)
    for i in 1:length(obs_data)
        mat[obs_inds[i],samp_inds[i]] = obs_data[i]
    end
    sum.(eachrow(mat)),mat .>0
end

function obs_df(biom,gens)
    has_gen = (gens .!="") .& (.!ismissing.(gens))
    obs_ids = biom["observation/ids"]
    samp_ids = biom["sample/ids"]
    depth,pa = obs_mat(biom)
    obs_names = obs_ids[has_gen] .* "_genus_" .* gens[has_gen]
    df = DataFrame(pa[:,has_gen],obs_names)
    df.sample_name = samp_ids
    df,depth 
end

function genus_df(df,gens)
    has_gen = (gens .!="") .& (.!ismissing.(gens))
    ugens = unique(gens[has_gen])
    n_samp = size(df,1)
    umat = zeros(Bool,n_samp,length(ugens))
    pa = Array(df[:,1:end-1])
    for (i,gen) in enumerate(ugens)
        umat[:,i] = Bool.(any.(eachrow(pa[:,gens[has_gen] .== gen])))
    end
    umat = Bool.(umat)
    seq_inds = [findfirst(x -> x == gen,gens[has_gen]) for gen in ugens]
    obs_names = DataFrames.names(df)[1:end-1][seq_inds] 
    udf =DataFrame(umat,obs_names)
    udf.sample_name = df.sample_name
    udf
end

# Fishers noncentral hypergeometric Likelihood function

function _FNCH_mode(ns,nf,n,ω)
    A = ω - 1
    B = n - nf - (ns + n + 2)*ω
    C = (ns + 1)*(n + 1)*ω
    return -2C / (B - sqrt(B^2 - 4A*C))
end

FNCH_mode(ns,nf,n,ω) = floor(Int, _FNCH_mode(ns,nf,n,ω))

function FNCH_pdf(ns,nf,n,ω, k)

    ω, _ = promote(ω, float(k))
    -100000000000<ω<100000000000 ? nothing : return zero(ω)
    l = max(0, n - nf)
    u = min(ns, n)
    η = _FNCH_mode(ns,nf,n,ω)
    isfinite(η) ? η=floor(Int, η) : return zero(ω)
    s = one(ω)
    fᵢ = one(ω)
    fₖ = one(ω)
    for i in (η + 1):u
        rᵢ = (ns - i + 1)*ω/(i*(nf - n  + i))*(n - i + 1)
        fᵢ *= rᵢ
        # break if terms no longer contribute to s
        sfᵢ = s + fᵢ
        if sfᵢ == s && i > k
            break
        end
        s = sfᵢ

        if i == k
            fₖ = fᵢ
        end
    end
    fᵢ = one(ω)
    for i in (η - 1):-1:l
        rᵢ₊ = (ns - i)*ω/((i + 1)*(nf - n  + i + 1))*(n - i)
        fᵢ /= rᵢ₊
        # break if terms no longer contribute to s
        sfᵢ = s + fᵢ
        if sfᵢ == s && i < k
            break
        end
        s = sfᵢ

        if i == k
            fₖ = fᵢ
        end
    end

    return fₖ/s
end

FNCH_logpdf(ns,nf,n,ω, k) = log(FNCH_pdf(ns,nf,n,ω, k))

# Network functions

@model function affinity_model(N,mA,mB,k,sig)
    mu ~ Normal(0.0,sig)
    try
    Turing.@addlogprob! FNCH_logpdf(mA,N-mA,mB,exp(mu),k)
    catch
     Turing.@addlogprob! -Inf
     end
end

infer(mod,alg,iter) = sample(mod,alg,iter)
@memoize function conf_and_est(N,mA,mB,k,sig)
    chn = infer(affinity_model(N,mA,mB,k,sig),NUTS(100,0.65),1000)
        x=vec(chn[:mu])
        sum(x .>0.0)/length(x),median(x)
        end

function individual_alpha(PA::Array;sig = 4.0)
    n_sample, n_obs = size(PA)
    res = Array{Float64}(undef,n_obs,n_obs)
    est = Array{Float64}(undef,n_obs,n_obs)
    N = n_sample
    sig = 4.0
    for j in 2:n_obs
        mB = sum(PA[:,j])
        for i in 1:j-1
            mA = sum(PA[:,i])
            k = sum(PA[:,i] .* PA[:,j])
            ma,mb = min(mA,mB), max(mA,mB)
            alpha,val = conf_and_est(N,ma,mb,k,sig)
            res[i,j] = alpha
            res[j,i] = alpha
            est[i,j] = val
            est[j,i] = val
            end
    end
    return res,est
end

function make_graphs(res,est,cutoff = 0.05)
    n_obs = size(res,1)
    neg = falses(n_obs,n_obs)
    pos = falses(n_obs,n_obs)
    res_both = zeros(n_obs,n_obs)
    for j in 2:n_obs
        for i in 1:j-1
            p = res[i,j]
            if p > 1.0-cutoff/2
                pos[i,j] = true
                pos[j,i] = true
                res_both[j,i] = res[i,j]
                res_both[i,j] = res[i,j]
            elseif p<cutoff/2
                neg[i,j] = true
                neg[j,i] = true
                res_both[j,i] = 1-res[i,j]
                res_both[i,j] = 1-res[i,j]
        end
        end
    end
    both = pos .| neg
    return SimpleGraph(neg),SimpleGraph(pos),SimpleGraph(both),1 .-res[LowerTriangular(neg)],res[LowerTriangular(pos)],res_both[LowerTriangular(both)],
    est[LowerTriangular(neg)],est[LowerTriangular(pos)],est[LowerTriangular(both)]
end


function plot_prep(g,classes;alpha =0.4)
    nodecolors = [(clr_dict[class],alpha) for class in classes]
    empties = degree(g) .==0
    g_copy = copy(g)
    rem_vertices!(g_copy,vertices(g)[empties],keep_order=true)
    nodecolors=nodecolors[.!empties]
    return g_copy,nodecolors
end
