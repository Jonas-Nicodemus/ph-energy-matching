var documenterSearchIndex = {"docs":
[{"location":"PortHamiltonianModelReduction/","page":"PortHamiltonianModelReduction","title":"PortHamiltonianModelReduction","text":"CurrentModule = PortHamiltonianModelReduction","category":"page"},{"location":"PortHamiltonianModelReduction/#PortHamiltonianModelReduction","page":"PortHamiltonianModelReduction","title":"PortHamiltonianModelReduction","text":"","category":"section"},{"location":"PortHamiltonianModelReduction/","page":"PortHamiltonianModelReduction","title":"PortHamiltonianModelReduction","text":"Documentation for PortHamiltonianModelReduction.","category":"page"},{"location":"PortHamiltonianModelReduction/","page":"PortHamiltonianModelReduction","title":"PortHamiltonianModelReduction","text":"Modules = [PortHamiltonianModelReduction]","category":"page"},{"location":"PortHamiltonianModelReduction/#PortHamiltonianModelReduction.bt-Tuple{ControlSystemsBase.StateSpace, Any}","page":"PortHamiltonianModelReduction","title":"PortHamiltonianModelReduction.bt","text":"Σr = bt(Σ::StateSpace, r; Lx=grampd(Σ, :o), Ly=grampd(Σ, :c)')\n\nReduces the state dimension of the system Σ to r using standard square root balanced truncation. The cholesky factors of the Gramians can passed as optional arguments Lx and Ly. The default values are the cholesky factors of the observability and controlability Gramians.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianModelReduction/#PortHamiltonianModelReduction.phirka-Tuple{PortHamiltonianStateSpace, Any, Any}","page":"PortHamiltonianModelReduction","title":"PortHamiltonianModelReduction.phirka","text":"Σr = phirka(Σph::PortHamiltonianStateSpace, r, num_runs; tol=1e-3, max_iter=50)\n\nCalls phirka num_runs times and returns the best result.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianModelReduction/#PortHamiltonianModelReduction.phirka-Tuple{PortHamiltonianStateSpace, Any}","page":"PortHamiltonianModelReduction","title":"PortHamiltonianModelReduction.phirka","text":"Σr = phirka(Σph::PortHamiltonianStateSpace, r; tol=1e-3, max_iter=50)\n\nReduces the state dimension of the port-Hamiltonian system Σph to r using the iterative rational Krylov algorithm.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianModelReduction/#PortHamiltonianModelReduction.prbt-Tuple{ControlSystemsBase.StateSpace, Any}","page":"PortHamiltonianModelReduction","title":"PortHamiltonianModelReduction.prbt","text":"Σr = prbt(Σ, r)\n\nReduces the state dimension of the system Σ to r using positve real balanced truncation, which is passivity preserving. Σ can be a StateSpace or a PortHamiltonianStateSpace.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems","text":"CurrentModule = PortHamiltonianSystems","category":"page"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems","text":"","category":"section"},{"location":"PortHamiltonianSystems/","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems","text":"Documentation for PortHamiltonianSystems.","category":"page"},{"location":"PortHamiltonianSystems/","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems","text":"Modules = [PortHamiltonianSystems]","category":"page"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.PortHamiltonianStateSpace","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.PortHamiltonianStateSpace","text":"PortHamiltonianStateSpace{T}\n\nAn object representing a port-Hamiltonian state-space system.\n\ndx(t)/dt = (J-R)Qx(t) + (G-P)u(t)\ny(t)     = (G+P)'Qx(t) + (S-N)u(t)\n\nSee the function phss for a user facing constructor.\n\nFields:\n\nJ::Matrix{T}\nR::Matrix{T}\nQ::Matrix{T}\nG::Matrix{T}\nP::Matrix{T}\nS::Matrix{T}\nN::Matrix{T}\n\n\n\n\n\n","category":"type"},{"location":"PortHamiltonianSystems/#ControlSystemsBase.gram-Tuple{PortHamiltonianStateSpace, Symbol}","page":"PortHamiltonianSystems","title":"ControlSystemsBase.gram","text":"X = gram(Σph, opt; kwargs...)\n\nReturns the Gramian of the system Σph. If opt is :c or :o  it returns the controllability or observability Gramian, respectively by calling the  gram from ControlSystems package.  If opt is :pr_c or :pr_o it returns the  positive-real controllability or observability Gramian, respectively (see prgram).\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#ControlSystemsBase.grampd-Tuple{PortHamiltonianStateSpace, Symbol}","page":"PortHamiltonianSystems","title":"ControlSystemsBase.grampd","text":"L = grampd(Σph, opt; kwargs...)\n\nReturns a Cholesky factor L of the Gramian of the system Σph. If opt is :c or :o  it returns the controllability or observability Gramian, respectively by calling the  grampd from ControlSystems.jl package.  If opt is :pr_c or :pr_o it returns a Cholesky factor L of the  positive-real controllability or observability Gramian, respectively.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#ControlSystemsBase.ss-NTuple{7, Any}","page":"PortHamiltonianSystems","title":"ControlSystemsBase.ss","text":"Σ = ss(J, R, Q, G, P, S, N)\nΣ = ss(Σph)\n\nConverts a PortHamiltonianStateSpace to a standard StateSpace.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#LinearAlgebra.norm","page":"PortHamiltonianSystems","title":"LinearAlgebra.norm","text":"norm(Σph, p=2; kwargs...)\n\nConverts a PortHamiltonianStateSpace to a StateSpace and calls norm on it. For more details see ControlSystems.jl package.\n\n\n\n\n\n","category":"function"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.compose-NTuple{7, Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.compose","text":"A, B, C, D = compose(J, R, Q, G, P, S, N)\n\nComposes the port-Hamiltonian matrices according to standard state-space matrices according to\n\nA = (J - R) * Q\nB = G - P\nC = (G + P)' * Q\nD = S - N.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.decompose-NTuple{5, Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.decompose","text":"J, R, Q, G, P, S, N = decompose(A, B, C, D, X)\n\nDecomposes the standard state-space matrices to port-Hamiltonian matrices according to\n\nQ = X\nJ = skew(A / X)\nR = -sym(A / X)\nG = 0.5 * (X \\ C' + B)\nP = 0.5 * (X \\ C' - B)\nS = sym(D)\nN = skew(D).\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.duplication-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.duplication","text":"D = duplication(n)\n\nReturns the duplication matrix of size n x n.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.ispsd-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.ispsd","text":"ispsd(M)\n\nReturns true if M is positive semi-definite, otherwise false.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.isskew-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.isskew","text":"isskew(M)\n\nReturns true if M is skew-symmetric, otherwise false. \n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.issym-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.issym","text":"issym(M)\n\nReturns true if M is symmetric, otherwise false.    \n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.kyp-Tuple{ControlSystemsBase.StateSpace}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.kyp","text":"X = kyp(Σ; kwargs...)\n\nReturns a solution to the KYP inequality by solving the corresponding linear matrix inequality.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.kypmat-Tuple{ControlSystemsBase.StateSpace, Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.kypmat","text":"W = kypmat(Σ, X)\n\nReturns the KYP matrix of the system Σ for the given matrix X.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.kypmax-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.kypmax","text":"Xmin = kypmax(Σ; kwargs...)\n\nReturns the maximal solution to the KYP inequality by solveing the Riccati equation for the anti-stabilizing solution.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.kypmin-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.kypmin","text":"Xmin = kypmin(Σ; kwargs...)\n\nReturns the minimal solution to the KYP inequality by solveing the Riccati equation for the stabilizing solution.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.lrcf-Tuple{Any, Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.lrcf","text":"lrcf(X, trunc_tol)\n\nComputes an approximate low-rank Cholesky-like  factorization of a symmetric positive semi-definite matrix X s.t. X = Z^T Z (up to a prescribed tolerance trunc_tol).\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.phss-Tuple","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.phss","text":"Σph = phss(J, R, Q, G, P, S, N)\nΣph = phss(J, R, Q, G)\nΣph = phss(Γ, W, Q)\n\nCreates a port-Hamiltonian state-space model Σph::PortHamiltonianStateSpace{T} with matrix element type T.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.phss-Tuple{ControlSystemsBase.StateSpace}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.phss","text":"Σph = phss(Σ)\nΣph = phss(Σ, X)\n\nConverts a StateSpace to a PortHamiltonianStateSpace by executing decompose. If X is not provided, the minimal solution of the KYP inequality is used.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.prare-Tuple{ControlSystemsBase.StateSpace, Symbol, Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.prare","text":"prare(Σ, opt, X; kwargs...)\n\nEvaluates the positive-real algebraic Riccati equation for system Σ and candidat solution X.\n\nIf opt is :o the positive-real controllability algebraic Riccati equation is evaluated,\n\nA'X + XA + (C' - XB) inv(D + D') (C - B'X).\n\nIf opt is :c the positive-real observability algebraic Riccati equation is evaluated,\n\nAX + XA' + (B - XC') inv(D + D') (B' - CX).\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.prgram-Tuple{ControlSystemsBase.StateSpace, Symbol}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.prgram","text":"X = prgram(Σ, opt; kwargs...)\n\nReturns the positive-real Gramian of system Σ. If opt is :c or :o  it returns the positive-real controllability or positive-real observability Gramian, respectively, by solving the corresponding positive-real algebraic Riccati equation (see prare).\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.prgrampd-Tuple{ControlSystemsBase.StateSpace, Symbol}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.prgrampd","text":"L = prgrampd(Σ, opt; kwargs...)\n\nReturns the Cholesky factor of the positive-real Gramian of system Σ (see prgram). If opt is :c or :o  it returns the positive-real controllability or positive-real observability Gramian, respectively.\n\nIn the case that the solution of the positive-real algebraic Riccati equation is not positive definite (due to numerical errors),  it is projected to the set of positive semi-definite matrices calling project_psd.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.project_psd-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.project_psd","text":"Mpsd = project_psd(M; eigtol=1e-8)\n\nReturns the nearest positive semi-definite matrix to M by setting negative eigenvalues to eigtol.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.skew-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.skew","text":"Mskew = skew(M)\n\nReturns the skew-symmetric part of a matrix M.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.sym-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.sym","text":"Msym = sym(M)\n\nReturns the symmetric part of a matrix M.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.truncation-Tuple{Any, Any, Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.truncation","text":"truncation(d, L, trunc_tol) -> (dr, Lr)\n\nComputes a rank revealing factorization for a given LDL-decomposition of S = L * mathrmdiag(d) * L^T (up to a prescribed tolerance trunc_tol) such that L_r * diag(d_r) * L_r^T approx S.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.unvec-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.unvec","text":"M = unvec(v)\n\nReturns the matrix M from the vectorized v, i.e., the inverse of vec.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.unvech-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.unvech","text":"M = unvech(v)\n\nReturns the symmetric matrix M from the half-vectorized v, i.e., the inverse of vech.\n\n\n\n\n\n","category":"method"},{"location":"PortHamiltonianSystems/#PortHamiltonianSystems.vech-Tuple{Any}","page":"PortHamiltonianSystems","title":"PortHamiltonianSystems.vech","text":"v = vech(M)\n\nReturns the lower triangular part of M as a vector. Aka the half-vectorization of M. For the inverse operation see unvech.\n\n\n\n\n\n","category":"method"},{"location":"EnergyMatching/","page":"EnergyMatching","title":"EnergyMatching","text":"CurrentModule = EnergyMatching","category":"page"},{"location":"EnergyMatching/#EnergyMatching","page":"EnergyMatching","title":"EnergyMatching","text":"","category":"section"},{"location":"EnergyMatching/","page":"EnergyMatching","title":"EnergyMatching","text":"Documentation for EnergyMatching.","category":"page"},{"location":"EnergyMatching/","page":"EnergyMatching","title":"EnergyMatching","text":"Modules = [EnergyMatching]","category":"page"},{"location":"EnergyMatching/#EnergyMatching.barrier-Tuple{QuadraticOutputStateSpace, ControlSystemsBase.StateSpace}","page":"EnergyMatching","title":"EnergyMatching.barrier","text":"Σph = barrier(Σ, Σr; Σr0=nothing, kwargs...)\n\nSolves the energy matching problem using the barrier method.\n\n\n\n\n\n","category":"method"},{"location":"EnergyMatching/#EnergyMatching.bestricc-Tuple{QuadraticOutputStateSpace, ControlSystemsBase.StateSpace}","page":"EnergyMatching","title":"EnergyMatching.bestricc","text":"Σphr = bestricc(Σ::QuadraticOutputStateSpace, Σr::StateSpace; kwargs...)\n\nSolves the energy matching problem using the best solution of the riccati equation.\n\n\n\n\n\n","category":"method"},{"location":"EnergyMatching/#EnergyMatching.constraint-Tuple{ControlSystemsBase.StateSpace}","page":"EnergyMatching","title":"EnergyMatching.constraint","text":"fg = constraint(Σ::QuadraticOutputStateSpace, Σr::StateSpace)\n\nReturns the (barrier-)constraint function and its gradient for the energy matching problem, in the standard format for Optim.jl.\n\n\n\n\n\n","category":"method"},{"location":"EnergyMatching/#EnergyMatching.hdss-NTuple{7, Any}","page":"EnergyMatching","title":"EnergyMatching.hdss","text":"Σqo = hdss(J, R, Q, G, P, S, N)\nΣqo = hdss(Σph)\n\nConverts a PortHamiltonianStateSpace to a QuadraticOutputSystem (Hamiltonian dynamic).\n\n\n\n\n\n","category":"method"},{"location":"EnergyMatching/#EnergyMatching.matchnrg-Tuple{PortHamiltonianStateSpace, PortHamiltonianStateSpace}","page":"EnergyMatching","title":"EnergyMatching.matchnrg","text":"Σrem = matchnrg(Σ::PortHamiltonianStateSpace, Σr::PortHamiltonianStateSpace; solver=:BFGS, kwargs...)\n\nSolves the energy matching problem.\n\n\n\n\n\n","category":"method"},{"location":"EnergyMatching/#EnergyMatching.objective-Tuple{QuadraticOutputStateSpace, ControlSystemsBase.StateSpace}","page":"EnergyMatching","title":"EnergyMatching.objective","text":"fg = objective(Σ::QuadraticOutputStateSpace, Σr::StateSpace)\n\nReturns the objective function and its gradient for the energy matching problem, in the standard format for Optim.jl.\n\n\n\n\n\n","category":"method"},{"location":"EnergyMatching/#EnergyMatching.sdp","page":"EnergyMatching","title":"EnergyMatching.sdp","text":"Σphr = sdp(Σ::QuadraticOutputStateSpace, Σr::StateSpace; optimizer=COSMO.Optimizer, ε=1e-8, kwargs...)\n\nSolves the energy matching problem using semidefinite programming.\n\n\n\n\n\n","category":"function"},{"location":"#Energy-matching-in-reduced-passive-and-port-Hamiltonian-systems","page":"Home","title":"Energy matching in reduced passive and port-Hamiltonian systems","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This repository contains the code for the paper Energy matching in reduced passive and port-Hamiltonian systems. The goal is to obtain low-dimensional port-Hamiltonian (pH) models that effectively approximate both the input-output dynamics and the energy (Hamiltonian) of a full-order model (FOM).","category":"page"},{"location":"#Citing","page":"Home","title":"Citing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you use this project for academic work, please consider citing our publication:","category":"page"},{"location":"","page":"Home","title":"Home","text":"T. Holicki, J. Nicodemus, P. Schwerdtner, and B. Unger\nEnergy matching in reduced passive and port-Hamiltonian systems\nArXiv e-print 2309.05778, 2023.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This code base is using the Julia Language and DrWatson to make a reproducible scientific project named","category":"page"},{"location":"","page":"Home","title":"Home","text":"ph-energy-matching","category":"page"},{"location":"","page":"Home","title":"Home","text":"To (locally) reproduce this project, do the following:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Download this code base. Notice that raw data are typically not included in the git-history and may need to be downloaded independently.\nOpen a Julia console and do:\njulia> using Pkg\njulia> Pkg.add(\"DrWatson\") # install globally, for using `quickactivate`\njulia> Pkg.activate(\"path/to/this/project\")\njulia> Pkg.instantiate()","category":"page"},{"location":"","page":"Home","title":"Home","text":"This will install all necessary packages for you to be able to run the scripts and everything should work out of the box, including correctly finding local paths.","category":"page"},{"location":"","page":"Home","title":"Home","text":"You may notice that most scripts start with the commands:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using DrWatson\n@quickactivate \"ph-energy-matching\"","category":"page"},{"location":"","page":"Home","title":"Home","text":"which auto-activate the project and enable local path handling from DrWatson.","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The executable script main.jl is located in the scripts directory.  It performs the following steps:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Set up the experiment.\nApply minimal realization (Kalman-like decomposition).\nDeclare and run the methods (Reductors and EnergyMatcher).\nEvaluate the ROMs.\nAnalyze the results.","category":"page"},{"location":"#Project-structure","page":"Home","title":"Project structure","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This project contains four packages:","category":"page"},{"location":"","page":"Home","title":"Home","text":"EnergyMatching.jl: Contains the methods for solving the energy matching problem.\nPortHamiltonianModelReduction.jl: Contains two structure preserving model reduction algorithms for port-Hamiltonian systems, phirka and prbt.\nPortHamiltonianSystems.jl: Contains methods for the analysis of port-Hamiltonian systems, as well as the PortHamiltonianStateSpace data type.\nQuadraticOutputSystems.jl: Contains methods for the analysis of linear dynamical systems with quadratic output, as well as the QuadraticOutputStateSpace data type.","category":"page"},{"location":"#License","page":"Home","title":"License","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Distributed under the MIT License. See LICENSE for more information.","category":"page"},{"location":"#Contact","page":"Home","title":"Contact","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Jonas Nicodemus - jonas.nicodemus@simtech.uni-stuttgart.de","category":"page"},{"location":"","page":"Home","title":"Home","text":"Tobias Holicki - tobias.holicki@imng.uni-stuttgart.de\nPaul Schwerdtner - paul.schwerdtner@nyu.edu\nBenjamin Unger - benjamin.unger@simtech.uni-stuttgart.de","category":"page"},{"location":"QuadraticOutputSystems/","page":"QuadraticOutputSystems","title":"QuadraticOutputSystems","text":"CurrentModule = QuadraticOutputSystems","category":"page"},{"location":"QuadraticOutputSystems/#QuadraticOutputSystems","page":"QuadraticOutputSystems","title":"QuadraticOutputSystems","text":"","category":"section"},{"location":"QuadraticOutputSystems/","page":"QuadraticOutputSystems","title":"QuadraticOutputSystems","text":"Documentation for QuadraticOutputSystems.","category":"page"},{"location":"QuadraticOutputSystems/","page":"QuadraticOutputSystems","title":"QuadraticOutputSystems","text":"Modules = [QuadraticOutputSystems]","category":"page"},{"location":"QuadraticOutputSystems/#QuadraticOutputSystems.QuadraticOutputStateSpace","page":"QuadraticOutputSystems","title":"QuadraticOutputSystems.QuadraticOutputStateSpace","text":"QuadraticOutputStateSpace{T}\n\nAn object representing a quadratic output state space system.\n\ndx(t)/dt = Ax(t) + Bu(t)\ny(t)     = Cx(t) + M(x(t)⊗x(t))\n\nSee the function qoss for a user facing constructor.\n\nFields:\n\nA::Matrix{T}\nB::Matrix{T}\nC::Matrix{T}\nM::Matrix{T}\n\n\n\n\n\n","category":"type"},{"location":"QuadraticOutputSystems/#ControlSystemsBase.gram-Tuple{QuadraticOutputStateSpace, Symbol}","page":"QuadraticOutputSystems","title":"ControlSystemsBase.gram","text":"X = gram(Σqo::QuadraticOutputStateSpace, opt::Symbol; kwargs...)\n\nReturns the Gramian of system Σqo. If opt is :c  it returns the controllability Gramian by calling the  ControlSystems.gram from ControlSystems package. If opt is :o it returns the observability Gramian by solving the Lyapunov equation\n\nA'X + XA + MPM = 0,\n\nwhere P is the controllability Gramian of Σqo.\n\n\n\n\n\n","category":"method"},{"location":"QuadraticOutputSystems/#ControlSystemsBase.grampd-Tuple{QuadraticOutputStateSpace, Symbol}","page":"QuadraticOutputSystems","title":"ControlSystemsBase.grampd","text":"L = grampd(Σqo::QuadraticOutputStateSpace, opt::Symbol; kwargs...)\n\nReturns a Cholesky factor L of the Gramian of system Σqo. If opt is :c  it returns the controllability Gramian by calling the  ControlSystems.grampd from ControlSystems package. If opt is :o it returns the Cholesky factor of the observability Gramian by solving the Lyapunov equation\n\nA'X + XA + MPM = 0,\n\nwhere P is the controllability Gramian of Σqo.\n\n\n\n\n\n","category":"method"},{"location":"QuadraticOutputSystems/#ControlSystemsBase.lsim-Tuple{QuadraticOutputStateSpace, Function, AbstractVector}","page":"QuadraticOutputSystems","title":"ControlSystemsBase.lsim","text":"result = lsim(Σqo, u, t; kwargs...)\n\nCalculate the time response of the quadratic output state-space model Σqo::QuadraticOutputStateSpace{T} by first treating it as standard state-space model and calling ControlSystems.lsim on it, and then calculating the quadratic output y_2(t) as y_2(t) = M(x(t)⊗x(t)) and add it to the linear output.\n\n\n\n\n\n","category":"method"},{"location":"QuadraticOutputSystems/#LinearAlgebra.norm","page":"QuadraticOutputSystems","title":"LinearAlgebra.norm","text":"norm(Σqo, p=2; kwargs...)\nh2norm(Σqo; kwargs...)\n\nComputes the H2 norm of the quadratic output state-space model Σqo::QuadraticOutputStateSpace.\n\n\n\n\n\n","category":"function"},{"location":"QuadraticOutputSystems/#QuadraticOutputSystems.h2inner","page":"QuadraticOutputSystems","title":"QuadraticOutputSystems.h2inner","text":"h2inner(Σ1, Σ2)\n\nComputes the H2 inner product of the quadratic output state-space models  Σ1::QuadraticOutputStateSpace and Σ2::QuadraticOutputStateSpace.\n\n\n\n\n\n","category":"function"},{"location":"QuadraticOutputSystems/#QuadraticOutputSystems.qoss-Tuple","page":"QuadraticOutputSystems","title":"QuadraticOutputSystems.qoss","text":"Σqo = qoss(A, B, C, M)\nΣqo = qoss(A, B, M)\n\nCreates a quadratic output state-space model Σqo::QuadraticOutputStateSpace{T} with matrix element type T.\n\n\n\n\n\n","category":"method"},{"location":"API/","page":"API","title":"API","text":"","category":"page"}]
}
