markershapes = [:none, :auto, :circle, :rect, :star5, :diamond, :cross, :xcross, :utriangle, :dtriangle, :rtriangle, :ltriangle, :pentagon, :star4, :star6, :star7, :star8, :vline, :hline, :+, :x]

abstract type AbstractMethod end

struct Reductor <: AbstractMethod 
    name::String
    reduce::Function
    label::LaTeXString
    color::Int
    marker_shape::Symbol
    linestyle::Symbol
end

struct EnergyMatcher <: AbstractMethod 
    name::String
    rom_name::String
    matchnrg::Function
    label::LaTeXString
    color::Int
    marker_shape::Symbol
    linestyle::Symbol
end

# reductor
phirka_method = Reductor("phirka", (fom, r) -> phirka(fom, r, 5), L"\textsf{pHIRKA}", 1, markershapes[4], :solid) 
prbt_method = Reductor("prbt", prbt, L"\textsf{PRBT}", 3, markershapes[3], :solid)

# energy matcher
# prbt
em_prbt_br_method = EnergyMatcher("em-prbt-br", prbt_method.name, (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:BestRicc), L"\textsf{PRBT}(X^\star)", 5, markershapes[5], :dash)
em_prbt_method = EnergyMatcher("em-prbt", prbt_method.name, (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:BFGS), L"\textsf{EM-PRBT}", 4, markershapes[6], :dash)
em_prbt_hypatia_method = EnergyMatcher("em-prbt-hypatia", prbt_method.name, (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:Hypatia), L"\textsf{EM-PRBT-Hypatia}", 2, markershapes[9], :dash)
em_prbt_cosmo_method = EnergyMatcher("em-prbt-cosmo", prbt_method.name, (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:COSMO), L"\textsf{EM-PRBT-COSMO}", 7, markershapes[10], :dash)
em_prbt_mosek_method = EnergyMatcher("em-prbt-mosek", prbt_method.name, (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:Mosek), L"\textsf{EM-PRBT-MOSEK}", 8, markershapes[5], :dash)

# phirka
em_phirka_method = EnergyMatcher("em-phikra", phirka_method.name, (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:BFGS), L"\textsf{EM-pHIRKA}", 7, markershapes[6], :dash)
em_phirka_hypatia_method = EnergyMatcher("em-phirka-hypatia", phirka_method.name, (fom, rom) -> EnergyMatching.matchnrg(fom, rom; solver=:Hypatia), L"\textsf{EM-pHIRKA-Hypatia}", 9, markershapes[9], :dash)
