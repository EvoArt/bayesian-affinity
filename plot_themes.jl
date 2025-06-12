c1 ="#f5c12f"
c2 ="#e13228"
c3 ="#1B6AA5"
c4 =RGB(([0,47,74] ./255)...)
cmap =(c1,c2,c3,c4)

bg = :white
fig_font = "Arial"
basic_theme = Theme(
           backgroundcolor = bg,
           Figure = (backgroundcolor = bg,),
           Scene = (backgroundcolor = bg,),
           Axis = (
               titlealign = :left,
               xgridvisible   = false,
               ygridvisible   = false,
               xlabelsize = 26,
               titlefont = fig_font,
               titlesize = 26,
               xticklabelfont = fig_font,
               xticklabelsize = 16,
               yticklabelfont = fig_font,
               ylabelfont = fig_font,
               xlabelfont = fig_font,
               ylabelsize = 26,
               yticklabelsize = 16,
               backgroundcolor= bg
           ),
           Axis3 = (
               titlealign = :left,
               xgridvisible   = false,
               ygridvisible   = false,
               xlabelsize = 26,
               titlefont = fig_font,
               titlesize = 26,
               xticklabelfont = fig_font,
               xticklabelsize = 20,
               yticklabelfont = fig_font,
               zticklabelfont = fig_font,
               ylabelfont = fig_font,
               xlabelfont = fig_font,
               zlabelfont = fig_font,
               zlabelsize = 26,
               backgroundcolor= bg
           ),
           Legend = (
           backgroundcolor = bg,
             titlefont = fig_font,
             labelfont = fig_font,
             labelsize = 20,
             titlesize = 20,
           ),
           Colorbar = (
             titlefont = fig_font,
             labelfont = fig_font,
             labelsize = 20,
             titlesize = 20,
           )
       )
       