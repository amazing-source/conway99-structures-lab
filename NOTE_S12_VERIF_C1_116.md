# S12 — Vérification indépendante de notes/116 de C1 (second ordre, critère lointain, Fano tendu)

25/09/2026, labo « structures », à la demande de C1. Preuves refaites à la main ; contrôle numérique sur de vrais
graphes de la famille srg(v,k,1,2), à savoir srg(9,4,1,2) et BvLS srg(243,22,1,2) (`scripts/verify_c1_p2_far.py`,
`logs/verify_c1_p2_far.log`). **Verdict : les trois points TIENNENT.** Il y a une précision à inscrire en (1).

## 1. Second sous-constituant et (P2)

**Équations de fentes.** A_Δχ_u = 2·1 − χ_u − χ_ū, en trois cas : y ∋ u (λ = 1), y ∋ ū, et y ∌ u, ū (μ = 2).

**Identité.**
- Pour y ≠ z, (A_Δ²)_{yz} = 2 − [y ~ z] − |y ∩ z|. Pour y = z, elle vaut k − 2. D'où
  **A_Δ² + A_Δ − (k − 2)I = 2J − ∂ᵀ∂**, soit 12I pour 99.
- Le sous-espace im ∂ᵀ est invariant, avec les valeurs propres 12, 0⁷, (−2)⁶. Sur ker ∂, qui est aussi invariant,
  on a θ² + θ = 12.
- Par la trace, spec Δ_x = 12¹ 3⁴⁰ 0⁷ (−2)⁶ (−4)³⁰. Contrôle : tr A_Δ² = 1008.

**(P2).** On encadre l'identité par 1_R et 1_R′, pour des ensembles R, R′ quelconques de lointains :
Σ_z cnt(z,R)cnt(z,R′) = (k−2)|R ∩ R′| − N(R,R′) + Σ_{y∈R,w∈R′}(2 − |y∩w|). Les termes y = w sont nuls.

**Diagonale sur un bloc (4 paires)** : **Σ cnt = 4(k − 2) et Σ cnt(cnt − 1) = 16 − 2e(R).**
- Le « 16 » de C1 suppose donc **e(R) = 0**, c'est-à-dire la torsion maximale, que note/E3 énonce bien.
- Contrôle : sur Paley 9 et sur les 55 blocs de BvLS, e(R) = 4, et la diagonale vaut 8. La forme générale est exacte.

## 2. Critère lointain
Soient L une droite, W = L ∪ N(L), F les lointains, M = A[N(L),F] et G = G(L), cubique.
- **A_F Mᵀ = Mᵀ(J/3 − A_G − I)** : on lit A² = (k−2)I + 2J − A sur F × N(L), et chaque lointain a 6 voisins dans N(L).
- **MMᵀ = (k−2)I + 2J − A_G − A_G² − K**, où K vaut J sur chaque part. On a A_G 1_{X_p} = 1_N, donc A_GK = KA_G = J, et
  tout commute.
- Valeurs propres de MMᵀ :
  - 6k − 24 = 60 sur 1 ;
  - 0 sur U (indicatrices de parts de somme nulle, où A_G = 0) ;
  - (r − A_G)(A_G − s) = (3 − A_G)(A_G + 4) sur le reste.
- ker M est invariant par A_F et vaut K_r ⊕ K_s (vecteurs propres de A portés par F).
- On écrit f = P_M f + f₀. Alors fᵀA_F f = nᵀ(MMᵀ)⁺(J/3 − A_G − I)n + r‖f₀‖² − (r − s)‖Π_s f‖², avec n = Mf. Le
  profil n n'a pas de composante sur U (2 voisins par part), et n_1 = (2|S|/(k−2))1 = (|S|/6)1.
- On obtient les **formes exactes** :
  - (r − s)‖Π_s 1_S‖² = r|S| + nᵀYn − 2e(S), avec Y = (MMᵀ)⁺(J/3 − A_G − (1 + r)I) ;
  - (r − s)‖Π_r 1_S‖² = 2e(S) − s|S| − nᵀY′n, avec Y′ = (MMᵀ)⁺(J/3 − A_G − (1 + s)I).
- **Pour 99 :**
  - sur 1, Y = 5/60 = 1/12 et Y′ = 12/60 = 1/5 ;
  - sur le complément, Y = −(3 − A_G)⁻¹ et Y′ = (A_G + 4)⁻¹ = Y + 7P ;
  - cela donne exactement (U) et (L), avec égalité pour S = F (480).
  Les versions par blocs 𝒰, 𝒯 et 𝒰 + 𝒯 = 28I − 7N̂ᵀPN̂ en découlent.
- Il faut G(L) connexe pour l'inverse sur le complément, comme C1 l'indique.
- Contrôle BvLS (r = 4, s = −5), sur 4 droites :
  - les identités de MMᵀ et de A_F Mᵀ sont exactes ;
  - les deux formes exactes tiennent sur 60 ensembles S aléatoires par droite ;
  - l'égalité vaut pour S = F.

## 3. Fano tendu : table des g, κ et (F5)
- **Table des g, re-dérivée** à partir de deux ingrédients :
  - les comptes des familles tendues : 2 dans la sienne, 1 si une classe est commune, 3 sinon ;
  - les équations de fente par classe : 4 voisins contenant la classe c, 4 contenant m.
  - Les près-voisins de y = {e,f} sont les deux τ, dans {ℓ,m} et {ℓ′,m}, et les deux σ, dans {ℓ,a} et {ℓ′,b}.
  - On en tire : g = 1 − [a ∈ ·] sur M_ℓ, 1 − [b ∈ ·] sur M_ℓ′, **3 − s_c − s_d sur M_m**, et g_{m,c} comme chez C1.
    Il faut pour cela que {c, M_ℓ(c)} ∪ {M_m(c), M_ℓ(M_m(c))} = C₄. ✓
- **κ ⟺ b ≠ M_m(a).**
  - Si b = M_m(a), alors par Klein M_ℓ′M_m = M_ℓ. Le 4ᵉ point c ∉ {a, b, M_ℓ(a)} a g_{m,c} = −1.
  - Réciproquement, les cas b ∈ {a, M_ℓ(a), M_ℓ′(a)} couvrent C₄. ✓
- **(F5).** Soit R un bloc de droite {c,d,p}, et R′₁, R′₂ ses deux blocs de type A (les droites par p).
  - **Σ = 16 sur les 8 paires de type A.** On a N(R, R′₁) + N(R, R′₂) = Σ_{y∈R}(6 − Σ_{C₄} s) = 4·4. En effet, par
    tension, les 8 atterrissages des paires de R tombent dans C₄(R) = R′₁ ∪ R′₂. La symétrie des arêtes
    N(R′,R) = N(R,R′) conclut.
  - **Toutes les autres paires ont cnt ≤ 1** :
    - R lui-même donne 0 (torsion maximale) ;
    - {c,p} et {d,p} donnent 1 (τ) ;
    - {c,x}, {d,x} et {p,x} donnent au plus 1 (famille tendue : un seul voisin dans Q_ℓ) ;
    - le type B donne 1 − [·].
  - La diagonale de (P2) donne alors Σ_A cnt(cnt − 1) = 16, d'où Σ_A (cnt − 2)² = 32 − 64 + 32 = 0. Donc **g = 2
    partout**.
  - Autrement dit, a et b tombent dans les deux blocs M_m différents, soit b ∈ {M_ℓ(a), M_ℓ′(a)}. ✓
