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

## 6. Le système local (i)+(ii)+(iii) semble infaisable à m = 7 (27/09, nuit) — piste pour tuer C2 d'un coup
**Comptage heuristique** (indépendance naïve, CLT local).
- Il y a 60! ≈ 10^81,9 bijections Φ.
- (ii) : 121 contraintes indépendantes, de probabilité ≈ 0,33 chacune.
- (iii) : 54 à 66 contraintes indépendantes, de probabilité ≈ 0,3.
- Espérance : **≈ 10^−10 solution par type de bloc à m = 7**, contre 10^48 à m = 11 (BvLS existe) et 10^24 pour (i)+(ii) seuls
  (que l'on résout en 14 s).
- Aucune solution de (i)–(iii) n'a jamais été trouvée à m = 7, par aucun de nous.
- Tests d'ordre 2 et d'ordre 3 :
  - les moments M± = V(I ± R)Vᵀ sont PSD sur les 110 359 types ;
  - le tenseur T3[a,b,d] avec ses 3 marges est PL-faisable sur 3000 types.
  - Il n'y a donc pas de contradiction « agrégée ».

**Deux étapes** (`involution/iota_design.py`).
- (iii) ⟺ ι = ΦσΦ⁻¹ est un couplage des 60 paires non-μ dont les comptes croisés valent K₃.
  - Au bloc plein, c'est un **quasi-design** : 30 quadruplets, chaque paire d'étiquettes couverte 3 fois, les paires π 2 fois,
    les paires m 1 fois.
- Étape 2 : Φ = affectation σ-paires → ι-paires orientées, sous (ii).
- **Chaque ι trouvé meurt en étape 2, au niveau PL** :
  - bloc plein 6-cycle : 9/9, dont 5 avec le lemme de C1 imposé ;
  - étoile libre : 1/1.
- Pour le 6-cycle, le (ii) agrégé par blocs suffit déjà. Le mécanisme est « propre » ou « compagnon » (couverture 2) combiné
  avec le « reste » (couverture 4). Le niveau des profils de blocs est faisable.

**Lemme mono-bloc (C1, notes/134).** En un bloc p d'un D-cycle de longueur ≥ 3, notons Sp_p = {p±, m(p±)}. Au plus UN
quadruplet contient 3 étiquettes de Sp_p, et aucun ne les contient toutes les 4.
- Il ne tue qu'une partie des ι (1 sur 4).
- La mort générique est globale : l'IIS s'étend sur plusieurs blocs et utilise les orientations au niveau des étiquettes.
- C1 ne trouve pas d'argument uniforme. En ÉTL, les certificats analogues ont un surplus nul.

**Solveurs sur le modèle monolithique** (portable, ≤ 290 s) : tous inconnus.
- HiGHS sur les σ-paires (65 000 variables) ;
- CP-SAT : Φ seul, et Φ + T_FF ;
- étape 1′ (quasi-design + (ii) par blocs, 19 680 variables) en CP-SAT et en HiGHS, avec et sans les coupes du lemme ;
- CaDiCaL (55 000 variables, 124 000 clauses).
En branch and bound sur ι, le PL de l'étape 1′ meurt après ~12 paires fixées, à 0,2 s par PL.

**Reformulation (C1).** (i)+(ii) ⟺ une **double couverture orthogonale** du graphe cocktail CP₆. Les pages H_a (a = étiquette
de N(x)) sont des quasi-2-facteurs de K₁₂ − μ (deux chemins plus des cycles) d'extrémités {a, πa, μa, τπa}. On a
|H_a ∩ H_{a′}| = 1 si a′ ∉ {a, πa} et 0 si a′ = πa. Le vrai verrou est (iii) par-dessus.

**À faire (cloud).**
- Preuve exhaustive par type : branch and bound sur ι avec PL (npz exportés : s1p_fullblock_6cyc, s1p_fullblock_24), ou SAT
  (`sat_local.py`) sur plusieurs cœurs.
- Commencer par les 4 types de bloc plein.
- Tous les types (~10³ à 10⁴ classes) demanderont une méthode plus efficace.

## 7. Contrôle positif sur BvLS (27/09, nuit)
BvLS = srg(243, 22, 1, 2), avec les mêmes λ = 1 et μ = 2 et m = 11 blocs. C'est un graphe de Cayley sur F₃⁵, et σ : v ↦ −v n'a
que 0 comme point fixe (`involution/bvls_local_control.py`).
- Aux **11 blocs**, la vraie structure locale vérifie EXACTEMENT :
  - (i) Φ bijective, compatible avec σ ;
  - (ii) et (iii), avec les formules GÉNÉRALES K₂ = 2J − I − P_π − P_μ − P_{τπ} et K₃ = 2J − P_π − P_τ − P_{τμ} − P_{μτ} ;
  - les conditions ponctuelles L1, S1 et S2 de C1.
- Tous les blocs y sont pleins (μ = τ).
- L'étape 2 accepte le VRAI ι de BvLS au niveau PL (`bvls_step2_control.py`).
- Les formules et le mécanisme de mort des modèles locaux sont donc validés sur un vrai graphe. Les morts à m = 7 ne sont pas
  des artefacts de modélisation.
- ⚠ **Nuance (C1, 27/09)** : sur BvLS avec σ = −1, on a μ = τ = π en TOUT bloc, donc tous les blocs gauches sont frères. Alors
  P_μ = P_τ = P_π et P_{τπ} = P_{τμ} = P_{μτ} = I. Le contrôle ne distingue donc PAS τπ de πτ, ni les autres termes entre eux :
  c'est un cas dégénéré. Les formules générales reposent sur la dérivation à la main (voisins communs avec sommets de
  correction), refaite indépendamment par C1 et par nous, et non sur ce contrôle.

## 8. Calibrage honnête de la piste locale (27/09, ~2 h)
- **Types pleins et mixtes** : chaque ι trouvé meurt en étape 2 AU PL : 9/9 au bloc plein, et 2/2 sur des types mixtes
  aléatoires.
- **Étoiles libres (ÉTL, μ ∩ τ = ∅)** : l'étape 1 est plus lente (une paire sécante par étiquette, car K₃[b,b] = 2). Sur le type
  aléatoire rand_type_7_3, l'étape 2 est **PL-faisable**, et la PLNE (150 s) comme CP-SAT (90 s) restent inconnues.
  - Le comptage heuristique y est plus proche de 1 : ≈ 10^−4,5 sans la diagonale nulle de K₃.
  - C1 a noté un surplus nul de ses certificats en ÉTL.
- **Portée réaliste** : la piste locale tuerait probablement toutes les configurations contenant un bloc NON libre, et
  réduirait C2 au cas « tous les blocs sont des étoiles libres ». Ce n'est pas une mort de C2 d'un coup.
- **Cas ÉTL** :
  - chercher une vraie solution locale ((i)+(ii)+(iii)) sur le cloud, avec sat_local.py et un cfg {"mu": …, "tau": …} ;
  - si on en trouve une, c'est la première solution locale à m = 7, et C1 pourra y tester le q:7 (lignes rigides à det carré,
    résidu K anisotrope) ;
  - sinon, il faut un UNSAT par type.

## 9. Pour les 79 types clés (Σ = ∅, « tout un bloc ») — 27/09, ~2 h 45
- **Relecture de la réduction large-T et de l'échelle des pivots de C1** : correctes (relecture envoyée à C1 et C2).
  - Recompte indépendant (`recount_pivots.py`) : 0 : 1787, 2 : 958, 4 : 407, 6 : 106, 8 : 39, 12 : 11. Pertinentes à ≥ 6 : 79.
  - Profil (2,4) : 15 classes, dont 9 pertinentes.
- **Le lemme mono-bloc de C1 s'étend à tout bloc gauche p ANCRÉ en c et NON FRÈRE** (p+ et p− pivots, μp+ ≠ p−), quel que
  soit le type. Relu par C1 : correct. Pour un bloc frère, Sp = {p±} a une couverture 0 et la preuve ne s'applique pas ; le code
  n'y émet d'ailleurs aucune clause.
  - sur Sp_p = {p±, μp±}, la couverture est 2 et Σλ(Sp_p) = 13 − [digone] ;
  - aucun centre sécant n'est dans Sp_p, car K₃ est nul en diagonale aux pivots ;
  - donc a₃(p) + 3a₄(p) ≤ 1 (ou 0 pour un digone), ce qui donne des coupes valides sur les ι-paires pour les types à blocs
    ancrés.
- **Forme égalité (C1)** : a₃ + 3a₄ + Σ_{Q∈S_p}(k_Q − 1)⁺ + #{Q ∉ S_p : k_Q = 0} = 1 − [digone], où S_p désigne les 10
  σ-paires passant par p et k_Q = |Q ∩ Sp|.
  - **Digone ancré {p, q}, dans tout type** : chaque σ-paire Z vérifie k_Z = 2 − |blocs(z) ∩ {p, q}|. C'est l'équilibre de
    D-cycle (v), désormais valable hors bloc plein.
    - Preuve : k = 1 sur les 8 σ-paires de S_p ∖ S_q (et symétriquement), k = 0 sur les 2 de {p, q}, k = 2 sur les 12
      autres.
    - Codé : `--digon`, dans sat_local.py et sat_local2.py. Il touche 21 des 79 types, dont tier8_0.
  - **Bloc chemin** (ancré, non frère, compagnons dans deux blocs) : toute σ-paire passant par p a k ≤ 2, et au plus une a
    k = 2.
- **Calibrage SAT (portable, ≤ 285 s)** : tier8_0 est inconnu avec l'encodage 1 (CaDiCaL) comme avec l'encodage 2
  (ι-paires explicites, `sat_local2.py`, Glucose). La difficulté est intrinsèque : elle se tranchera sur le cloud. Lancer les
  deux encodages en parallèle et garder le premier UNSAT (même contenu logique).
