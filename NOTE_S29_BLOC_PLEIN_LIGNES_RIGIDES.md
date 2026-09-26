# NOTE S29 : le bloc plein (bloc A), ses ombres lointaines, les lignes rigides et le bloc impair lointain

26/09, nuit. Cadre : involution σ de srg(99,14,1,2), Fix σ = {x} ; bloc c **plein** (toute orbite de l'étoile de c pivote en c,
c'est le « bloc A », d_c = 6 en AD). Code : `involution/fullblock_rigid.py`, `involution/fullblock_cpsat.py`,
`involution/tff_check.py`, `ad/step1_fullblock.py`. DÉMONTRÉ = preuve écrite ci-dessous.

**Notations.**
- Étiquettes v = (b, ±) des 6 autres blocs ; π(v) est l'étiquette opposée.
- m : compagnons dans N(c+), c'est-à-dire (c+, v) ~ (c+, m v).
- h : (c+, v) ~ (c−, h v).
- Racine lointaine z = {a₁, a₂} (deux blocs ≠ c) ; Φ(z) = {v : z ~ (c+, v)}.
- Ombre s_f = 1_{Φ(z)} − 1_{Φ(σz)} ; d_f = 1_z − 1_{σz} ; B_f = signes de z par bloc.

## 1. Faits de base au bloc plein (DÉMONTRÉ)
- **h = π∘m**, c'est-à-dire τ_c = η_c.
  - o_v pivote en c. Ses voisins par c sont (c+, m v) et (c−, h v), et ils forment D(o_v) = o_{m v}. D'où h v = π m v.
  - En particulier S_{o_w, o_v} = [m w = v] − [h w = π v] = 0 : le bloc étoile × étoile de S est nul.
- **D préserve étoile(c).** Le partenaire D d'une orbite lointaine est donc lointain.
- **Pas d'annulation.** Si z ~ (c+, l) et z ~ (c−, πl), alors n(z, o_l) = 2, donc o_l = D(f) serait dans l'étoile : exclu.
  - Donc Φ(z) ∩ Φ(σz) = ∅, |supp s_f| = 4, et chaque ligne lointaine de T a exactement 6 entrées ±1 parmi les 29 autres
    orbites lointaines (10 − 4).
- **(i)–(iii) valent en tout bloc**, sans hypothèse d'étoile libre. On compte les voisins communs, avec des sommets de
  correction explicites et deux à deux distincts :
  - Φ est une bijection des racines lointaines sur les paires non-m ;
  - K₂ = 2J − I − P_π − P_m − P_{τπ} et K₃ = 2J − P_π − P_τ − P_{τm} − P_{mτ}, avec τ = π∘h.
  - Au bloc plein (τ = m, sans frère) :
    - K₂[a,b] = 2 − [a=b] − [b=πa] − [b=ma] − [b=mπa] ;
    - K₃[b,d] = 2 − [d=mb] − [d=πb] − 2[d=b].
- **Unicité du voisin par a** (a ∈ z) : on n'a jamais à la fois a ∈ Φ(z) et πa ∈ Φ(σz), ni à la fois πa ∈ Φ(z) et
  a ∈ Φ(σz).

## 2. Orbites lointaines ancrées ; critère de rigidité (DÉMONTRÉ)
- **Zéros d'ancrage.** Soit f dans e₁e₂, ancrée en e₁. Ses voisins par e₁ sont les racines de D(f), qui est lointain. Donc
  Φ(z) et Φ(σz) évitent lab(e₁).
- **Formule (C1, relue).** Au bloc plein, Z_ff ≡ 4x + y (mod 7), avec :
  - x = s·d = Σ_i (α_i − γ_i), où α_i = [le voisin de z par a_i passe par c] et γ_i = [le voisin par πa_i passe par c] ;
  - y = Σ_b s(b+) s(b−).
- f est rigide ⟺ (x, y) ∈ {(0,0), (2,−1), (−2,1)}.
- **Si f est ancrée en e₁**, alors x = α₂ − γ₂ ∈ {−1, 0, 1}.
  - Le cas α₂ = γ₂ = 1 force deux étiquettes de e₂ de même signe, donc y ∈ {1, 2}.
  - Il reste : **f ancrée est rigide ⟺ son ombre a exactement une étiquette dans chacun des 4 blocs hors de {c, e₁, e₂}**,
    c'est-à-dire φ_f(c e₂) = 0 et φ_f(c b) = 1 pour ces 4 blocs.
  - C'est la forme des ombres frères du 3K₂ de C1.
- **Rep 5827**, bloc A = 6 en 6-cycle (0,1,5,2,4,3). Une solution du modèle par paires a 3 orbites ancrées rigides :
  (0,5)₀, (2,3)₀ et (2,5)₀. Le niveau paires ne force PAS la rigidité en général (la plupart des ancrées ont φ(c e₂) = 1).

## 3. Le bloc impair lointain T = T_FF au bloc plein (DÉMONTRÉ)
S² + S = 12I, lue sur les lignes lointaines et les colonnes de la boule, donne des équations **linéaires exactes** en T :
- (E1) Σ_g T_fg s_g = −(s_f + d_f) ;
- (E2) Σ_g T_fg B_g = −C s_f, où (C s)(e) = s(e+) − s(e−) ;
- (E3) diagonale : 6 entrées non nulles par ligne ; hors diagonale : B_f·B_g + s_f·s_g + (T²)_fg + T_fg = 0.
- Autrement dit X_F T = N X_F, avec X_F = (B ; s) de taille 18 × 30 et N = −(S_II′ + I) de spectre (−1)⁶ (−2)⁶ 1⁶.

**Spectre de T (exact).** On a Σ_f s_f s_fᵀ = 11I − J + P_π (entrées étoile × étoile), donc les lignes-étiquettes de X_F sont
de somme nulle.
- T agit sur col(X_Fᵀ) comme N sur ℝ¹⁸/ker X_Fᵀ, et sur ker X_F avec T² + T = 12I.
- tr T = 0 et tr T² = 180 imposent 12j = Σ μ(μ+1) ≤ 2j pour les j valeurs propres de N perdues au-delà de 1_labels, donc
  j = 0.
- D'où **rang X_F = 17** et **spec T = 3⁹ (−4)⁴ (−1)⁵ (−2)⁶ 1⁶**.
- C'est un graphe signé 6-régulier sur les 30 orbites lointaines.

**Résidu réel au bloc plein (forme close).** T = A_Φ − 7Π, où A_Φ = X_Fᵀ(X_F X_Fᵀ)⁺ N X_F + 3Π_ker est déterminé par
Φ, et Π est un projecteur de rang 4 dans ker X_F.
- Avec la formule de C1 (4|r_f|² = 10 − w_fᵀM_c⁺w_f, où M_c = 11I + P_π − J au bloc plein et w_f = 4s_f + d_f), on
  obtient **|r_f|² = (14 + y − 3x)/15**.
- Contrôles : Σ_f |r_f|² = 28 (car Σx = Σy = 0), et |r_f|² ≡ 4x + y ≡ Z_ff (mod 7).
- ⚠ Une orbite rigide a |r_f|² = 14/15 (7-adiquement divisible par 7), PAS zéro. La rigidité est un énoncé purement 7-adique
  (T_{f·} ≡ Z_{f·} mod 7). Il n'existe pas de version réelle « ligne exacte » : une première rédaction de cette note le
  prétendait à tort.

## 4. Essais de calcul (portable, ≤ 5 min, un seul processus)
- Modèle local Φ au bloc plein du rep 5827 (6-cycle + ancrages) :
  - HiGHS (45 696 variables σ-paires) : inconnu à 200 s ;
  - CP-SAT (2 cœurs) : inconnu à 240 s.
  - Même difficulté que l'ÉTL : la bijection Φ est dure pour les solveurs génériques.
- Modèle joint Φ + T (E1, E2, 6 par ligne) en CP-SAT : en cours.

## 5. Au passage (relectures de cette nuit)
- **Sol, A111111 (3,3)** : le signe de Witt est juste.
  - Rang 14 ; les 4389 mineurs principaux 14 × 14 non singuliers sont carrés, donc la boule est non scindée ; g ≥ 24 − 2r₀.
  - Énergie 3‖ψ‖²/2 re-dérivée ; certificats revérifiés par C2.
  - **A33 est mort.**
- **Lemme de Witt général**, en tout bloc : g ≥ 18 + t − 2r₀ + 2·[boule réduite non scindée].
  - C1 en tire **g⁻ ≠ 16** (notes/132). J'ai relu ses trois points (décomposition de rang, Witt, (i)–(iii) en tout bloc) :
    ils tiennent.
  - **g⁻ ∈ {18, 20}.**
- Un bloc plein force g = 20. À g = 18, il n'y a donc aucun bloc plein, et les résidus AD ne vivent qu'à g = 20.
