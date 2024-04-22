# FibonacciFit

The following example session shows how to run the `FibonacciFit` muon reconstruction.

```julia
julia> using NeRCA, KM3io

julia> det = Detector("aux/KM3NeT_00000133_20221025.detx")
Detector 133 (v5) with 21 strings and 399 modules.

julia> ffparams = FibonacciFitParameters(;tmaxlocal=18.0, roadwidth=200.0, nfits=12, α₁=7.0, α₂=0.4, θ=7.0)
FibonacciFitParameters(18.0, 200.0, 50, 12, 10, 5.0, 7.0, 0.4, 7.0)

julia> ffit = FibonacciFit(ffparams, det)
FibonacciFit with a coarse scan of 7.0ᵒ and a fine scan of 0.4ᵒ.

julia> f = ROOTFile("/Volumes/Ultraspeed/HPSS_Lyon/mc/atm_muon/KM3NeT_00000133/v8.1/trigger/mcv8.1.mupage_tuned_100G.sirene.jterbr00013288.1.root")
ROOTFile{OnlineTree (11412 events, 11412 summaryslices), OfflineTree (11412 events)}

julia> ffit(f.online.events[23].snapshot_hits)
12-element Vector{FibonacciFitCandidate}:
 FibonacciFitCandidate([88.93184137346843, 339.09585624318083, 418.6361293065551], [0.17810965449491756, -0.10848711561882224, -0.9780120125644721], 1.5834566902250916e7, 47.940130047899885, 0.5426162591010238, 0.044859575569164145, 2, 8)
 FibonacciFitCandidate([89.08091581215234, 340.34586759852624, 418.5547927780048], [0.17663799377277853, -0.11617474934909236, -0.9773957472639233], 1.5834565772132114e7, 47.936745453848154, 0.5426162591010238, 0.044859575569164145, 2, 8)
 FibonacciFitCandidate([88.71738426137938, 339.51767701352196, 418.6530740654206], [0.17263474610610405, -0.11889742891861323, -0.9777835373094723], 1.5834566628640838e7, 47.93555523513469, 0.5426162591010238, 0.044859575569164145, 2, 8)
 FibonacciFitCandidate([88.68018117411658, 338.55166370495215, 418.7048867808974], [0.17503881844497432, -0.11116961496780249, -0.9782651627985623], 1.583456748187883e7, 47.89273846705427, 0.5426162591010238, 0.044859575569164145, 2, 8)
 FibonacciFitCandidate([88.92348748131037, 340.8201543426498, 418.5087870382899], [0.17431113629732115, -0.12319422516471291, -0.9769538426398674], 1.5834565063683022e7, 47.88922147604587, 0.5426162591010238, 0.044859575569164145, 2, 8)
 FibonacciFitCandidate([88.49903181162912, 339.9858817112237, 418.6188045359057], [0.16938556379379227, -0.1268537689226887, -0.9773518568501152], 1.583456596506952e7, 47.88315627435971, 0.5426162591010238, 0.044859575569164145, 2, 8)
 FibonacciFitCandidate([88.33500463892133, 336.6297736412863, 418.6465039846342], [0.18351420416559014, -0.09413314536360323, -0.9784996105330982], 1.5834567917471658e7, 47.8771730128378, 0.5426162591010238, 0.044859575569164145, 2, 8)
 FibonacciFitCandidate([88.6112184882073, 340.88428643669346, 418.4948246456611], [0.1712486202541263, -0.1296335434685557, -0.9766622007986405], 1.5834564676049871e7, 47.841189208589036, 0.5426162591010238, 0.044859575569164145, 2, 8)
 FibonacciFitCandidate([89.06821069843838, 338.7434590937738, 418.5691401666228], [0.18478520379540103, -0.09669122983426234, -0.9780108560396608], 1.5834566795501318e7, 47.83543683516487, 0.5426162591010238, 0.044859575569164145, 2, 8)
 FibonacciFitCandidate([88.21610397084605, 336.1470720660109, 418.57912176091025], [0.18714443443904788, -0.08772505956646029, -0.9784075196882687], 1.5834567757131675e7, 47.810986334354645, 0.5426162591010238, 0.044859575569164145, 2, 8)
 FibonacciFitCandidate([96.7134749624555, 359.1848068570001, 418.4664433214492], [0.21393497522746657, -0.059264166957793095, -0.9750485038649214], 1.5834558094776593e7, 47.788389553327626, 0.5426162591010238, 0.044859575569164145, 2, 8)
 FibonacciFitCandidate([89.54414464834538, 340.52834817554526, 418.50586547857824], [0.1836347574587258, -0.1039257041807162, -0.977485408516985], 1.5834565731202858e7, 47.779737246769365, 0.5426162591010238, 0.044859575569164145, 2, 8)
```
