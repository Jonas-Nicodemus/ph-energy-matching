colors = ColorBrewer.colorSchemes["Set1"]["9"]
markershapes = [:square, :circle, :dtriangle, :star, :diamond, :utriangle, :plus]

abstract type AbstractMethod end

struct Reductor <: AbstractMethod 
    name::String
    reduce::Function
    label::String
    color::String
    marker_shape::Symbol
    linestyle::Symbol
end

struct EnergyMatcher <: AbstractMethod 
    name::String
    matchnrg::Function
    label::String
    color::String
    marker_shape::Symbol
    linestyle::Symbol
end

# reductor
phirka_method = Reductor("phirka", (fom, r) -> phirka(fom, r, 5), "\\textsf{pH-IRKA}", colors[2], markershapes[1], :solid) 
prbt_method = Reductor("prbt", (fom, r, Lx, Ly) -> phss(bt(ss(fom), r; Lx=Lx, Ly=Ly)), "\\textsf{PRBT}", colors[3], markershapes[2], :solid)

# energy matcher
em_prbt_br_method = EnergyMatcher("em-prbt-br", (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:BestRicc), "\\textsf{PRBT}(X^\\star)", colors[4], markershapes[4], :dash)
em_prbt_cosmo_method = EnergyMatcher("em-prbt-cosmo", (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:COSMO), "\\textsf{EM-PRBT-COSMO}", colors[4], markershapes[4], :dash)
em_prbt_mosek_method = EnergyMatcher("em-prbt-mosek", (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:Mosek), "\\textsf{EM-PRBT-MOSEK}", colors[7], markershapes[6], :dash)
em_prbt_sedumi_method = EnergyMatcher("em-prbt-sedumi", (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:SeDuMi), "\\textsf{EM-PRBT-SeDuMi}", colors[8], markershapes[7], :dash)
em_prbt_bfgs_method = EnergyMatcher("em-prbt-bfgs", (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:BFGS), "\\textsf{EM-PRBT}", colors[5], markershapes[5], :dash)
em_prbt_bfgschol_method = EnergyMatcher("em-prbt-bfgs-chol", (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:BFGSchol), "\\textsf{EM-PRBT-Chol}", colors[6], markershapes[5], :dash)
em_prbt_bfgs_paul_method = EnergyMatcher("em-prbt-bfgs-paul", (fom, r) -> pauls_roms(fom, r), "EM-PRBT-BFGS-Paul", colors[1], markershapes[1], :dash)


