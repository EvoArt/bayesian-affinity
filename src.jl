using Turing, Optim, RCall, Random, CairoMakie, GLMakie
Random.seed!(123)
R"library(CooccurrenceAffinity)"

# Helper functions

# Calculate errors
rmse(x,x̂) = √mean((x .- x̂) .^2)
# Alpha_mle from integers
aints(N,mA,mB,k) = rcopy(R"ML.Alpha($k,c($mA,$mB,$N), lev=0.9)")[:est]
# Turing map estimates
tmap(N,mA,mB,k,prior) = optimize(alpha_bayes(N,mA,mB,k,prior),MAP()).values[:alpha]
@model function alpha_bayes(N,mA,mB,k,prior)
alpha ~ Normal(0.0,prior)
k ~ FisherNoncentralHypergeometric(mA,N-mA,mB,exp(alpha))
end

# Generate data for Figure 1
ry = Float64[]
ty = Float64[]
x = Float64[]
mb = Float64[]
ma = Float64[]
N = 30
mA = 15
for mA in [5,15,25]
    for mB in [1,5,10,15,25,29]
        for a in -2:0.1:2
            for rep in 1:10
                try
                    k = rand(FisherNoncentralHypergeometric(mA,N-mA,mB,exp(a)))
                    push!(ty,tmap(N,mA,mB,k,3.0))
                    push!(ry,aints(N,mA,mB,k))
                    push!(x,a)
                    push!(mb,mB)
                    push!(ma,mA)
                    catch e
                end
            end
        end
    end
end
c3 ="#b2d5eb"
c3 ="#41b199"
c1 ="#f5c12f"
c2 ="#e13228"
# Plot Figure 1
fig = Figure()
for (j,mB) in enumerate([1,5,10,15,25,29])
    for (i,mA) in enumerate([5,15,25])
        mask = (ma .== mA) .& (mb .== mB)
        rmser = round(rmse(x[mask],ry[mask]), sigdigits=3)
        rmset = round(rmse(x[mask],ty[mask]), sigdigits=3) 
        rp = randperm(2sum(mask))   

    ax = Axis(fig[i,j], title = "mA = $mA, mB = $mB\nRMSEmle = $rmser\nRMSEmap = $rmset",
                xlabel = "Affinity", ylabel = "Estimated affinity")
                clr = vcat(fill(c1,sum(mask))...,fill(c2,sum(mask))...)
                scatter!(ax,vcat(x[mask]...,x[mask]...)[rp] ,vcat(ry[mask]...,ty[mask]...)[rp], alpha = 0.5, color = clr[rp],
        strokecolor = clr[rp],strokewidth = 0.5, markersize = 7)
    lines!(-2:2,-2:2,color = :black,linewidth = 4)
    if i < 3
        hidexdecorations!(ax)
    else
        hidexdecorations!(ax, ticks = false, label = false,ticklabels = false)
    end
    if j >1
        hideydecorations!(ax)
    else
        hideydecorations!(ax, ticks = false, label = false,ticklabels = false)
    end
    end
end
linkaxes!(fig.content...)

# Regression

# Functions for linear models
regmle(x,y) = optimize(reg(x,y),MLE(), NelderMead()).values[:β]
@model function reg(x, y)
    α ~ Normal(0,4)
    β ~ Normal(0,4)
    σ ~ Exponential()
    for i in eachindex(x)
        μ = α + β*x[i]
        y[i] ~ Normal(μ,σ)
    end
end

# Functions for glm
regmle(x,N,mA,mB,k) = optimize(reg(x,N,mA,mB,k),MLE(), NelderMead()).values[:β]
@model function reg(x,N,mA,mB,k)
    # Priors are assigned to let Turing know that these are parameters to 
    # be estimated. But they are not used in the present analysis, since
    # `regmle` perform maximum likelihood estimation.
    α ~ Normal(0,4)
    β ~ Normal(0,1)
    for i in eachindex(x)
        μ = α + β*x[i]
        k[i] ~ FisherNoncentralHypergeometric(mA[i],N-mA[i],mB[i],exp(μ))
    end
end

# Generate data for figure 2
regr = Float64[]
regt = Float64[]
regb = Float64[]
b = Float64[]
for β in -3:0.15:3
    for rep in 1:1
        x = randn(30)
        y = β .* x
        mA = rand(1:29,30)
        mB = rand(1:29,30)
        N = 30
        k = rand.(FisherNoncentralHypergeometric.(mA,N .-mA,mB,exp.(y)))
        RY =aints.(N,mA,mB,k)
        TY =tmap.(N,mA,mB,k,3.0)
        success = false

        while success == false
            try
            rb = regmle(x,N,mA,mB,k)
            rt = regmle(x,TY)
            rr = regmle(x,RY)
            push!(regr,rr)
            push!(regt,rt)
            push!(b,β)
            push!(regb,rb)
            success = true
            catch e
            end
        end
    end
end

# Plot Figure 2
fig = Figure()
ax1 = Axis(fig[1,1], title = "MLE", xlabel = "β actual", 
                    ylabel = "β estimate", xticks = -3:3, yticks = -5:5)
lines!(ax1,b,b, color = :black, linewidth = 4)
scatter!(ax1,b,regr,color = c1)
hideydecorations!(ax1, ticks = false, label = false, ticklabels=false)
ax2 = Axis(fig[1,2], title = "MAP", xlabel = "β actual", ylabel = "β estimate", xticks = -3:3)
lines!(ax2,b,b, color = :black, linewidth = 4)
scatter!(ax2,b,regt,color = c2)
hideydecorations!(ax2)
ax3 = Axis(fig[1,3], title = "GLM", xlabel = "β actual", ylabel = "β estimate", xticks = -3:3)
lines!(ax3,b,b, color = :black, linewidth = 4)
scatter!(ax3,b,regb,color = c3)
hideydecorations!(ax3)
hidexdecorations!.(fig.content, ticks = false, ticklabels = false, label = false)
linkaxes!(fig.content...)


