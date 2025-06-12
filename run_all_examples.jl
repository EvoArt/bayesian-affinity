run(`quarto render "beetles.qmd"`)
run(`quarto render "data_processing_microbes.qmd"`)
run(`quarto render "microbes.qmd"`)
run(`quarto render "ants_traits.qmd"`)
run(`julia ants_dist.jl --threads = 8`)
                                                                                                         