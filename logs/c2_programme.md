# Branche C2 (involution) : programme strictement mathématique

Démarré le 2026-09-25 à la demande de l'opérateur : **aucun solveur (SAT, CP-SAT, recherche exhaustive)**. Les calculs machine ne servent
qu'à *tester* un lemme sur les deux graphes réels de la même famille qui possèdent une involution à un point fixe :
srg(9,4,1,2) (m = 2) et Berlekamp–van Lint–Seidel srg(243,22,1,2) (m = 11, `code/z2/control_bvls.py`).
Règle de cohérence : tout lemme général doit être vrai pour m = 2 et m = 11 ; une contradiction à m = 7 doit utiliser une
arithmétique propre à m = 7 (ici θ − τ = 7 = m).
Consigne opérateur : quand un raisonnement bute, inventer l'outil plutôt que chercher une solution existante.
Étiquettes : PROUVÉ / VÉRIFIÉ (sur m = 2, 11) / CONJECTURÉ / IMPASSE.

## 0. Notations (voir `NOTE_involution_case.md` §2)
Γ = srg(99,14,1,2), σ involution, Fix(σ) = {x}. Blocs B_i = {(i,0),(i,1)} (i = 1..7) de N(x) ; σ échange les deux sommets de
chaque bloc. Sommets extérieurs = étiquettes {a,b} (a, b dans des blocs différents), 84 ; orbites U = (ij,p), 42.
Représentant de U = (ij,p), i < j : u = {(i,0),(j,p)}. a_UW = A[u,w], b_UW = A[u,σw].
S (49×49, orbites intérieures puis extérieures) = [[−I7, Cᵀ],[C, S_out]], S² + S = 12I, spectre 3²⁷, (−4)²².
C (42×7) : C_Ui = +1 si u ∋ (i,0), −1 si u ∋ (i,1), 0 sinon. CᵀC = 12I, CᵀS_out = 0, S_out² + S_out = 12I − B, B = CCᵀ.
Q_out (42×42, côté invariant) : Q_UW = nombre de voisins de u dans l'orbite W ; Q_out² + Q_out = 12I + 4J − Sh (Sh = NbNbᵀ).
**G := Q_out − 2D** (ajouté le 25/09 après relecture) : graphe des orbites SANS les arêtes D, 10-régulier, et D est la matrice du couplage
des carrés. D'où G² + G = 8I + 4J − Sh − 2(GD + DG) − 2D ; les termes en D sont ceux qui produisent les « − 2e » du §15.

## 1. Structure locale (PROUVÉ)
**Lemme 1 (21 carrés).** Pour u extérieur, N(u) ∩ N(σu) = {v, σv} et (u, v, σu, σv) est un 4-cycle induit, σ-invariant, σ = rotation
d'un demi-tour. Les 84 sommets extérieurs sont partitionnés en 21 tels carrés (D est un couplage parfait des 42 orbites).
*Preuve.* μ = 2 et σ-invariance donnent {v,σv} ; v ≁ σv (NOTE §2(b)) ; D symétrique. ∎

**Lemme 2 (triangles extérieurs).** Les 3 étiquettes d'un triangle entièrement extérieur sont deux à deux disjointes.
*Preuve.* Si u ∩ w ∋ c, alors c est le sommet commun de u et w, donc le triangle unique de l'arête uw est {c,u,w}. ∎
Comptage : 7 triangles par x, 84 triangles {c,u,w} (un sommet intérieur), 140 triangles extérieurs (70 σ-orbites).

**Lemme 3 (deux couplages par bloc).** Soit O_i l'ensemble des 12 orbites passant par le bloc i, a = (i,0). Chaque u ∋ a a
exactement un voisin contenant a (couplage m_i) et exactement un voisin contenant a′ (appariement π_i : u ~ σw avec w ∋ a).
m_i et π_i sont deux couplages parfaits de O_i ; R_i = m_i ∪ π_i est 2-régulier (cycles pairs alternés + digones).
Les digones de R_i sont exactement les paires D contenues dans O_i. Pour tout U : Σ_{W∈O_i} Q_UW = 2.
En holonomie signée : une arête U–W de R_i a h = s_UW C_Ui C_Wi = +1 si elle est dans m_i, −1 si elle est dans π_i.

## 2. Un invariant 7-adique nouveau (PROUVÉ)
**Théorème A (CORRIGÉ le 25/09 après l'audit adverse d'Astra).** Pour toute srg(99,14,1,2), α := rank_F7((A − 3I)|_{1⊥}) est pair,
donc e′ := rank_F7(A − 3I) = α + 1 est **impair** (la version initiale disait « e′ pair » : faux, voir la preuve ; concorde avec L2 de C1).
**Théorème B.** Dans le cas involution, r := rank_F7(S − 3I) est pair ; de façon équivalente
e := rank_F7((S_out − 3I)|_{C⊥}) = r − 7 est impair (1 ≤ e ≤ 15).

*Outil (réseaux propres + formule d'oddité).* Soit X une matrice symétrique entière n×n, X² + X = 12I (valeurs propres 3, −4),
V₋ = ker(X + 4I), P₋ = (3I − X)/7 le projecteur orthogonal sur V₋, K₋ = Zⁿ ∩ V₋.
(i) 7P₋ est entière, donc 7K₋* = 7P₋(Zⁿ) ⊂ K₋ : le groupe discriminant est (Z/7)^d, et K₋*/K₋ ≅ Fⁿ₇ / ker(X̄ − 3) ≅ Im(X̄ − 3),
    donc d = rank_F7(X − 3I).
(ii) P₋ est 2-entière : K₋ ⊗ Z₂ est facteur direct orthogonal de Z₂ⁿ, donc unimodulaire.
(iii) w = P₋(1) est un vecteur caractéristique de K₋ ⊗ Z₂ (x·w = x·1 ≡ x·x mod 2), donc oddity(K₋) ≡ w·w = (3n − σ_X)/7 (mod 8),
    où σ_X = 1ᵀX1 = somme des coefficients.
(iv) Formule d'oddité (réseau défini positif de rang f = dim V₋) : oddity ≡ f + Σ_{p impair} p-excess (mod 8). Seul p = 7 intervient :
    K₋ ⊗ Z₇ = 1^{±(f−d)} ⊥ 7^{ε d}, 7-excess = 6d + 4k (k = 1 ssi ε = −1).
Donc (3n − σ_X)·7⁻¹ ≡ f + 6d + 4k (mod 8).

*Preuve de B.* X = S, n = 49, f = 22. Au-dessus de la diagonale, S a 7·12 + 210 = 294 coefficients ±1, donc
σ_S = tr S + 2Σ_{X<Y} S_XY ≡ −7 + 2·294 ≡ 1 (mod 4). Modulo 4 : w·w ≡ 3(3·49 − 1) ≡ 2, et 22 + 6d ≡ 2 + 2d. Donc d = r est pair.
Le bloc intérieur de S − 3I est −4I ≡ 3I (inversible mod 7), d'où r = 7 + rank(complément de Schur) ; ce complément est
Y = S_out − 3I + 2B (mod 7), nul sur span(C), égal à S_out − 3I sur C⊥, d'où e = r − 7. ∎
*Preuve de A.* X = A − (terme J) : on applique (i)–(iv) à K₋₄ = Z⁹⁹ ∩ E₋₄ (E₋₄ = (3I − A)/7 + J/63, 2-entière). **Erreur corrigée :**
63E₋₄ = −9(A − 3I) + J ≡ 2(3I − A) + J (mod 7), et non 2(3I − A) (le J manquait ; sur 1 l'ancienne formule donnait 0 d'un côté,
6·1 de l'autre). Plus proprement, 7E₋₄ ≡ (3I − A) + 4J (mod 7) : nul sur 1 (−11 + 396 = 385 ≡ 0), égal à 3I − A sur 1⊥ ; donc
d = α, alors que (A − 3I)1 = 11·1 ≢ 0 donne e′ = α + 1.
1 ∈ E₁₄ donc w = E₋₄(1) = 0 : K₋₄ est pair, oddity 0. La partie 3-adique est 9^{±1} (recollement avec Z·1, norme 99), 3-excess ≡ 8 ≡ 0.
Donc 0 ≡ 44 + 6α + 4k (mod 8) : α pair, **e′ impair**. ∎

*Contrôles.* m = 2 (Paley 9, p = 3, θ − τ = 3) : même argument, prédiction rank_F3(S − I) pair ; calcul à la main : rang 2. VÉRIFIÉ.
*Cohérence (IMPASSE partielle).* Les décompositions alternatives (K_C ⊥ K₃′ ⊥ K₋₄′ dans l'espace extérieur, partie invariante Q̃,
graphe extérieur H) redonnent les mêmes relations : toutes les formules d'oddité sont liées par l'invariance de l'oddité par
sous-réseau d'indice 2^a. Seule information nouvelle : r pair, et une relation (r/2 + k ≡ n₋ mod 2) entre le signe 7-adique k et
la parité du nombre n₋ d'arêtes négatives de Γ₁ dans la jauge standard.
Pour conclure par cette voie il faut une *seconde* détermination, indépendante, de r mod 2 ou de k.

## 3. Comptages de points fixes (IMPASSE, consigné pour ne pas refaire)
Pour une involution sur un ensemble fini X : |X| ≡ |X^σ| (mod 2). Essayé avec X = sommets, arêtes, triangles, non-arêtes,
4-cycles (2079 ; 21 fixes), 5-cycles (41 580 ; 0 fixe), 6-cycles (429 660 par Ihara–Bass ; 224 fixes : 42 réflexions par deux arêtes
de blocs + 42 hexagones antipodaux à une paire intérieure + 140 hexagones antipodaux extérieurs, ce dernier nombre étant forcé :
20 chemins u–q–r–σu extérieurs par u). Tout est pair et forcé : aucune information.

## 4. Raffinements 2-adiques et reformulation (PROUVÉ ; aucune contradiction)
- Module libre : K₋₄ ⊗ Z₂ ≅ Z₂[C₂]²² (Broué/Krull–Schmidt, χ₋₄(σ) = 0). Pour la moitié invariante β (forme x·(1+σ)y), un vecteur
  caractéristique w vérifie (1+σ)w = c := E₋₄(e_x), d'où β(w,w) = (E₋₄)_xx / 2 = 2/9 ≡ 2 (mod 8) ; det β est un carré.
  Avec oddity(K₋₄) = 0 on obtient **σ_S ≡ 1 (mod 8)**, puis **r/2 + k ≡ 0 (mod 2)**. Même résultat par E₃ (c·c/2 = 3/11).
- MAIS σ_S ≡ 1 (mod 8) se redémontre combinatoirement : pour chaque paire {U,W}, la chaîne P_UW des 2-chemins U–Z–W a une holonomie
  fixée par S_out² + S_out = 12I − B ; chaque arête apparaît 18 fois dans Σ P_UW, d'où Σ_{U<W}((G²)_UW + s_UW + B_UW)/2 ≡ 0,
  soit n₋ pair. Les réseaux ne donnent donc, au-delà de la combinatoire, que l'information 7-adique.
- **Reformulation (Théorème B′).** N := S − 3I mod 7 est symétrique, N² = 0, rang r pair, et la forme x ↦ xᵀNx sur F₇⁴⁹/ker N est
  *hyperbolique* (discriminant (−1)^{r/2}). Pour le graphe entier (Théorème A′, CORRIGÉ) : on travaille sur W = 1⊥, en dimension α = e′ − 1
  (l'exposant e′/2 initial n'était pas entier) ; avec α = 2m l'oddité donne m + k ≡ 1 (mod 2), et la forme induite par A − 3I sur W
  mod 7 est *non déployée* dans les deux cas (m pair : dét non carré ; m impair : dét carré, alors que (−1)^m = −1 n'est pas un carré mod 7). La partie σ-invariante est donc non déployée elle aussi.
  Une contradiction viendrait d'un sous-espace totalement isotrope explicite trop grand, ou d'un calcul indépendant du type de Witt.

## 5. Heuristique « 7 | k » (CONJECTURÉ comme piste)
Famille srg(1+2m², 2m, 1, 2), t = θ − τ = √(8m−7). Paramètres admissibles : m = 2, 7, 11, 56, 497 (t = 3, 7, 9, 21, 63).
Un premier p divise à la fois t et k = 2m ssi p = 7 (t² = 8m − 7) ssi 7 | m : c'est le cas exactement pour m = 7, 56, 497, les trois
paramètres dont l'existence est inconnue ; m = 2 et 11 (graphes existants, avec involution) ont t = 3, 9 premiers à 7.
Un argument 7-adique où « k ≡ 0 (mod 7) » intervient est donc le bon type d'argument (il ne peut pas s'appliquer à m = 2, 11).
Là où 7 | m intervient déjà : sur la partie invariante, Nbᵀ Nb = 10I + 2J₇ a les valeurs propres 24 et 10 congrues mod 7 (J₇ nilpotent mod 7).

## 6. Observations sur BvLS (m = 11) — à ne pas prendre pour des lemmes
- Les 55 carrés y sont tous du type « frère » (D(U) = U′). Les types « un bloc commun » et « disjoints » ne sont exclus par aucun
  argument local que j'ai trouvé (les identités Q², S² et le comptage des apex sont automatiquement satisfaits).
- G ∘ (DGD) = 0 (si U ~ W dans Γ₁ alors D(U) ≁ D(W)) ; [S_out, J] ≠ 0 et {S_out, J} ≠ 0 pour la structure complexe J des carrés.

## 7. Théorème D : le cas « tous les carrés frères » est impossible (PROUVÉ, 2026-09-25)
Notation : P = ij une paire de blocs ; ses 4 étiquettes {a,b},{a,b′},{a′,b′},{a′,b} forment deux orbites sœurs P₀, P₁.
Type de P : n_s(P) = Q_{P₀P₁} ∈ {0,1,2} ; n_s = 2 ⟺ D(P₀) = P₁ (carré « frère »).
**Théorème D.** Il n'existe pas de srg(99,14,1,2) avec involution σ telle que D(U) = U′ (orbite sœur) pour toute orbite U.
*Preuve.*
(1) Si D(U) = U′ et u = {a,b} : les deux étiquettes sœurs {a,b′} ∋ a et {a′,b} ∋ a′ sont voisines de u, donc m_a(u), P_a(u),
    m_b(u), P_b(u) sont ces deux étiquettes (condition (I) : un seul voisin contenant a, a′, b, b′). Les 10 autres voisins
    extérieurs de u ne rencontrent ni le bloc i ni le bloc j : ils sont dans des paires de blocs disjointes de P.
(2) Deux voisins y₁ ≠ y₂ de u dans une même paire P′ sont soit σ-conjugués (alors y₁ ∈ N(u) ∩ N(σu) = D(u) ⊂ P, absurde), soit
    partagent un sommet c ; dans un carré frère ils sont adjacents, et c, u seraient deux voisins communs d'une arête (λ = 1). Absurde.
(3) Il y a exactement C(5,2) = 10 paires disjointes de P et 10 voisins lointains : **u a exactement un voisin dans chaque paire
    disjointe** (c'est ici que m = 7 sert : 2m − 4 = C(m−2,2) ⟺ m = 7 ; pour m = 11, BvLS a 18 voisins pour 36 paires).
    Par (3) appliqué aux deux côtés, H[P,P′] est un couplage parfait σ-invariant, qui induit une bijection des orbites P → P′.
(4) Donc, en indexant les orbites par (P, ε) : Q_out = 2·(échange frère) + Σ_{P∩P′=∅} E_{PP′} ⊗ Π_{PP′}, Π ∈ {I₂, X}, et 0 entre
    paires partageant un bloc. Les espaces Sym (constants sur P) et Alt (h(P,ε) = (−1)^ε g(P)) sont invariants :
    Q_out|Sym = A(KG) + 2I et Q_out|Alt = A_χ − 2I, où KG = graphe de Kneser K(7,2) = srg(21,10,3,6), de spectre 10, 1¹⁴, (−4)⁶,
    et A_χ est une signature de KG (χ_PP′ = +1 si Π = I₂, −1 si Π = X).
(5) Spec(Q_out) = {12, (−2)⁶, 3²⁰, (−4)¹⁵} (NOTE §2(e) / DEEP_THOUGHT §1). Sym en prend {12, (−2)⁶, 3¹⁴}, donc Alt porte {3⁶, (−4)¹⁵} :
    A_χ a le spectre {5⁶, (−2)¹⁵}, soit A_χ² = 3A_χ + 10I.
(6) P, P′ disjointes ont 3 voisins communs dans KG : (A_χ²)_{PP′} = 3χ_PP′ force χ_PR χ_RP′ = χ_PP′ pour chacun : **tous les
    triangles de KG sont équilibrés**. P = 12, P′ = 13 (non adjacentes) ont 6 voisins communs, les paires de {4,5,6,7}, groupées en
    3 couples complémentaires {R, R̄} ; R ~ R̄ et les triangles (P,R,R̄), (P′,R,R̄) équilibrés donnent χ_PR χ_RP′ = χ_PR̄ χ_R̄P′.
    Donc (A_χ²)_{PP′} = 2(t₁ + t₂ + t₃) avec tᵢ = ±1 : c'est ±2 ou ±6, jamais 0. Or A_χ² = 3A_χ + 10I exige 0. Contradiction. ∎
*Cohérence.* Pour m = 11 l'étape (3) échoue (18 < 36) : aucun conflit avec BvLS, dont les 55 carrés sont frères.
**Conséquence.** Toute involution devrait avoir au moins un carré de type « un bloc commun » ou « disjoint ».
Types : sib ⟺ n_s = 2 ; « un bloc » ⟹ n_s = 0 ; « disjoint » ⟹ n_s ∈ {0,1} ; n_s = 1 ⟹ D de type disjoint.

## 8. Cas mixte : outils (PROUVÉ) et premières contraintes
Compression sur les paires de blocs. Sym = vecteurs constants sur chaque paire P (21), Alt = h(P,ε) = (−1)^ε g(P) (21) ;
Sym = span(Nb) ⊕ Sym′. T = compression de Q_out sur Sym (T_PP = n_s(P), T_PP′ = e(P,P′)/4, e = nb d'arêtes de H entre les 4+4
étiquettes), Q_A = compression sur Alt ((Q_A)_PP = −n_s(P), (Q_A)_PP′ = (q₀₀ + q₁₁ − q₀₁ − q₁₀)/2 avec q_ab = Q_{P_a P′_b}),
K : Alt → Sym′ (lignes indexées par Sym′, colonnes par Alt ; l'ancienne mention « Sym′→Alt » était inversée, audit d'Astra), K_PP′ = (q_P(P′₀) − q_P(P′₁))/2, q_P(Z) = nb de voisins de z dans P. Comme Nb⊥ est Q-invariant avec Q² + Q = 12 :
  T² + T + KKᵀ = (12I + 4J − Sh)|_Sym  et  Q_A² + Q_A + KᵀK = 12I   (identités exactes).
Sommes de lignes : Σ_{P′∩P=∅} T_PP′ = 8 + n_s(P), Σ_{|P′∩P|=1} T_PP′ = 4 − 2n_s(P).
**Lemme 8.1 (paires frères, PROUVÉ).** Si P est frère : toute étiquette disjointe de P a exactement un voisin dans P, et aucune
étiquette partageant un bloc avec P n'y a de voisin (c'est la perfection de la grille G_P, Lemme 9.1). Donc e(P,R) = 4 pour R
disjointe, T_PR = 1, et la ligne P de K est nulle. Si P et P′ sont frères et disjointes, H[P,P′] est un couplage parfait.
**CORRECTION (2026-09-25, même séance).** Une première version affirmait aussi que la *colonne* P de K est nulle, donc que chaque
étiquette de P a exactement un voisin dans *toute* paire disjointe R. C'est FAUX en général : pour R disjointe non frère, les deux
orbites de P peuvent avoir (2,0) voisins dans R (l'orbite à 2 voisins est alors le μ-partenaire d'une paire partageante non
adjacente de R). L'identité diagonale donne μ_U + μ_U′ = ‖K_{·P}‖², qui est *automatiquement* vraie dans ce cas. Les anciens
lemmes 8.2–8.4 (équilibre des triangles, contrainte sur les couplages de K₄) ne valent donc que si les paires R concernées sont
elles aussi frères ; l'argument de parité exige localement « tout frère ». RETIRÉS en l'état.
**Lemme 8.2 (diagonale, PROUVÉ).** Pour toute paire P : μ_U + μ_U′ + n_s(P) + ‖K_{P·}‖² − ‖K_{·P}‖² = 2, avec
μ_U = Σ_R Q_{UR₀}Q_{UR₁} ; sommé sur P : Σ_U μ_U = 42 − N_s (c'est l'entrée (Q²)_{UU′} = 2 − n_s, rien de neuf).
**Lemme 8.5 (graphe non frère, PROUVÉ).** Soit N ⊂ E(K₇) l'ensemble des paires non frères. Chaque bloc d'une paire X ∈ N porte
une autre paire de N (les voisins proches d'une étiquette de X ne peuvent pas être dans une paire frère partageant un bloc,
par 8.1) : N n'a pas de sommet de degré 1.
Reste à traiter : tout cas où N ≠ ∅.

## 9. Grilles et quasi-perfection (PROUVÉ)
**Lemme 9.1 (grilles parfaites, sans involution).** Si S est un ensemble de 9 sommets induisant 18 arêtes (p. ex. une grille
K₃□K₃ = srg(9,4,1,2)), alors ‖E₋₄1_S‖² = (3·9 − 36)/7 + 81/63 = 0, donc 1_S ∈ E₁₄ ⊕ E₃ et A1_S = 3·1_S + 1 : **tout sommet hors
de S a exactement un voisin dans S**. (Critère général : E₋₄f = 0 ⟺ Af = 3f + (1ᵀf/9)·1.)
Au point fixe x, P = ij est frère ⟺ G_P = {x, a, a′, b, b′} ∪ (4 étiquettes de P) est une grille ; le lemme redonne 8.1 et 8.4.
Pour P quelconque : ‖E₋₄1_{G_P}‖² = (8 − 4n_s(P))/7 et ⟨E₋₄1_{G_P}, E₋₄1_{G_P′}⟩ = (4/7)(T₀ − T)_{PP′} (T₀ = A(KG) + 2I) :
c'est exactement T ≼ T₀ (rien de neuf côté invariant).
**Lemme 9.2 (carré « un bloc commun » presque parfait).** Soit Q = (u,v,σu,σv), u = {a,b} ∈ 12, v = {a,d} ∈ 13, β(Q) = {1,2,3},
F_β = les 24 étiquettes des 6 paires de {4,5,6,7}, w l'apex extérieur de l'arête v–σu.
(i) une étiquette y ∉ Q a au plus un voisin dans Q, sauf les apex (a, a′, w, σw) qui en ont deux (voisins adjacents) ;
(ii) Σ_{y∈F_β} |N(y) ∩ Q| = 24 + 2(ω(u) + ω(v)), où ω(u) = nb de {m_b(u), P_b(u)} dont l'autre élément est dans le bloc 3
     (idem ω(v) avec d et le bloc 2). D'où ω(u) + ω(v) ≤ [w ∈ F_β] ≤ 1 : presque toute étiquette lointaine a exactement un voisin dans Q.
(iii) si ω(u) = 0 : u a 6 voisins lointains dans les 6 paires de {4,5,6,7} ; si ce multigraphe est simple c'est K₄ et les quatre
     « demi-arêtes » (autres éléments de m_b, P_b, y_d, y_d′) tombent dans les quatre blocs lointains, une chacune.
Les 6 paires lointaines d'un carré « un bloc » forment 3 couples complémentaires : c'est la configuration exacte de l'argument
de parité du Théorème D (P = 12, P′ = 13, R ⊂ {4,5,6,7}). Piste : un argument de parité local sur chaque carré non frère.

## 10. Théorème D local (PROUVÉ) et conséquence
Pour P frère et R disjointe : chaque orbite de R touche exactement une orbite de P, donc le bloc 2×2 est une permutation
((Q_A)_PR = ±1, K_RP = 0) ou du type (2,0) ((Q_A)_PR = 0, K_RP = ±1) ; (2,0) exige R non frère et une étiquette de P
μ-partenaire d'une paire partageante non adjacente de R. Si toutes les paires qui interviennent sont frères, les identités
exactes (Q_A² + Q_A + KᵀK)_{PX} = 12δ n'ont plus de correction et l'argument du Théorème D s'applique tel quel.
**Théorème D local.** S'il existe trois blocs t, a, b tels que ta et tb soient frères et que toutes les paires de [7]∖{t},
sauf éventuellement ab, soient frères, alors contradiction.
*Preuve.* P = ta, P′ = tb. Pour R ⊂ [7]∖{t,a,b} (6 paires, 3 couples complémentaires) : (Q_A²)_{PP′} = Σ_R χ_PR χ_RP′ car les
colonnes P, P′ de K sont nulles (toutes les paires disjointes de P ou de P′ sont dans [7]∖{t}, frères sauf ab qui partage
un bloc avec les deux) ; l'équilibre des triangles (P,Y,R̄) s'obtient comme au Théorème D par la ligne exacte de P en R̄
(tous les Y concernés sont dans [7]∖{t}, frères). D'où 2(t₁+t₂+t₃) = 0, absurde. ∎
**Remarque (IMPASSE partielle).** Sous l'hypothèse du Théorème D local, toute arête tc de N (c ∉ {a,b}) laisse c de degré 1
dans N, ce que 8.5 interdit ; donc N ⊂ {ab}, puis N = ∅ par 8.5. La version locale n'exclut donc rien de plus que 8.5 + Théorème D.
Ce que 8.5 + Théorème D excluent : N = ∅ et |N| ≤ 2 ; plus généralement tout N ayant un sommet de degré 1.
Plus petit N non exclu : un triangle de K₇ (p. ex. {12, 13, 23}). Là, pour deux paires frères partageant un bloc, chaque couple
complémentaire mélange frère et non frère : il faut comprendre l'interaction frère / non frère (motifs (2,0)).

## 11. Le graphe N des paires non frères (PROUVÉ)
Notations : N ⊂ E(K₇) (paires non frères, vues comme arêtes sur les 7 blocs), d_i = degré du bloc i dans N,
n_far(P) = nombre d'arêtes de N disjointes de P, φ_u(R) = nombre de voisins de u dans la paire R.
**Lemme 11.1 (comptage lointain).** Pour P ∈ N et u ∈ P : u a exactement un voisin dans chaque paire frère disjointe de P
(perfection de la grille, 9.1 : son seul voisin dans G_S est une étiquette de S), et 8 + n_s(P) voisins lointains en tout. Donc
  Σ_{R ∈ N, R∩P=∅} φ_u(R) = n_far(P) − 2 + n_s(P),   en particulier n_far(P) ≥ 2 − n_s(P) ≥ 1.
**Lemme 11.2.** Si D(U) est de type disjoint (dans R), φ_u(R) = 2, donc n_far(P) ≥ 4 − n_s(P) ≥ 3.
Les partenaires D des orbites non frères sont non frères : type « un bloc » ⟺ arêtes de N adjacentes, type « disjoint » ⟺ arêtes
de N disjointes.
**Lemme 11.3 (liaison).** Si d_i = 2 (arêtes ij, ik de N), la partie non frère de R_i vit sur les 4 orbites de ij ∪ ik ; si une
orbite de ij est D-appariée à ik à travers i, l'orbite sœur l'est aussi (R_i 2-régulier, digones saturés). De plus, sans digone
en i ni arête sœur, les 4 orbites forment un 4-cycle de R_i.
Rappel : une arête sœur simple U–U′ (n_s = 1) appartient à R_i ET à R_j (m dans l'un, π dans l'autre).
**Conséquences.**
- (C1) aucun bloc de degré 1 dans N (8.5) ; (C2) toute arête de N a au moins 2 arêtes de N disjointes, 3 si elle porte une orbite
  de type disjoint ; (C3) n_s = 1 ⟹ type disjoint.
- |N| ≤ 4 impossible. Avec au plus 3 arêtes et degré minimal ≥ 2 (C1), le seul graphe possible est le triangle, où aucune arête n'a
  d'arête disjointe (contre C2). Avec 4 arêtes et degré minimal ≥ 2, c'est C₄, où 12 n'a qu'une arête disjointe (contre C2).
  *(Corrigé le 25/09 après relecture de la session « structures » : l'ancienne rédaction « 3 arêtes ⟹ couplage » oubliait le triangle.)*
- 5 arêtes, degré minimal ≥ 2 : C₅ ou le losange K₄ − e. Le losange est exclu par (C2) : l'arête joignant ses deux sommets de degré 3
  n'a aucune arête disjointe. *(Cas ajouté le 25/09, même relecture.)*
- **N = C₅ impossible.** Chaque arête n'a que 2 arêtes disjointes, donc pas de type disjoint, donc tout est « un bloc » ; tous
  les degrés valent 2, donc (11.3) D induit un couplage parfait des 5 arêtes du cycle en arêtes adjacentes : impossible, 5 est impair.
- Triangle composante {ij, ik, jk} de N (d_i = d_j = d_k = 2) : un appariement « un bloc » dans le triangle force la troisième
  arête à être de type disjoint avec n_s = 0 (R_j et R_k saturés), donc n_far ≥ 4.
- Deux triangles disjoints {12,13,23} ∪ {45,46,56} : n_far = 3 partout, donc (point précédent) aucun « un bloc », tout est
  « disjoint » avec n_s = 1 sur les 6 paires, et chaque étiquette n'a comme voisins dans l'autre triangle que sa paire D.
  Survit aux comptages actuels : À TRAITER.
**Lemme 11.4 (parité d'orientation, PROUVÉ).** Une arête sœur simple de P = ij (n_s = 1) est de type m en un bout et π en
l'autre (u = {a,b} ~ {a,b′} : m pour M_a au bloc i, π au bloc j). Si le bloc i est de degré 2 dans N avec deux arêtes ij, ik
de n_s = 1, R_i est le 4-cycle formé des deux arêtes sœurs (opposées) et de deux arêtes croisées ; R_i alternant m/π, les deux
arêtes sœurs ont le même type en i. Donc le long d'un cycle de N dont tous les blocs sont de degré 2 et toutes les arêtes de
n_s = 1, le type s'inverse sur chaque arête et se conserve à chaque bloc : **un tel cycle est de longueur paire**.
**Deux triangles disjoints : impossible.** D'après §11, tout y est de type disjoint avec n_s = 1 ; chaque triangle est un cycle
impair de blocs de degré 2 : contradiction avec 11.4.
**Liste à 6 arêtes.** Les graphes à 6 arêtes, de degré minimal 2, dont chaque arête a au moins 2 arêtes disjointes sont C₆, 2C₃ et
K₂,₃ (K₄, papillon, C₅ + corde échouent au critère des arêtes disjointes). 2C₃ est exclu ; restent C₆ et K₂,₃.
**Lemme 11.5 (règle des 2-chemins entre sœurs, PROUVÉ).** Soit P = ij avec n_s(P) = 0, orbites sœurs U, U′. L'entrée (U,U′) de
S_out² + S_out = 12I − B donne Σ_Z s_UZ s_ZU′ = 0 (B_UU′ = 0 pour des sœurs), et côté Q, (G²)_UU′ = 2 − 2(G_{U,D(U′)} + G_{D(U),U′})
∈ {0, 2}. Si W est voisin commun de U et U′ par R_i (chemin U–W–U′ dans R_i : une arête m, une π, car R_i alterne), son 2-chemin
a le signe h₁h₂ C_Ui C_U′i = −C_Ui C_U′i ; par R_j c'est −C_Uj C_U′j = +C_Ui C_U′i (car C_Ui C_U′i C_Uj C_U′j = −1 pour des sœurs).
Donc : **U, U′ ont au plus un voisin commun par R_i et au plus un par R_j** ; deux voisins communs par le même R_i sont impossibles.
**Corollaire 11.6 (blocs de degré 2).** Si d_i = 2 (arêtes e, e′), alors n_s(e) = n_s(e′) (une arête sœur de e est dans R_i et
sature le côté e), et R_i est soit deux digones (paire liée « un bloc » à travers i, n_s = 0), soit le cas S (deux arêtes sœurs,
n_s = 1). Le 4-cycle sans arête sœur mettrait U, U′ à distance 2 deux fois dans R_i : exclu par 11.5.
**Corollaire 11.7 (composantes cycles).** Si une composante de N est un cycle (tous ses blocs de degré 2), n_s est constant sur
elle ; « tout lié » est impossible (une arête serait liée à ses deux voisines, or D(U) est unique), donc n_s = 1, cas S partout,
tout de type disjoint (n_far ≥ 3), et longueur paire (11.4). D'où à nouveau C₅ exclu, et **C₆ ⟹ n_s = 1, cas S partout, tout disjoint**.
**K₂,₃** (N = {13,14,15,23,24,25}) : n_far = 2 ⟹ pas de disjoint ⟹ n_s = 0 ; les blocs 3, 4, 5 sont de degré 2 ⟹ paires liées
13–23, 14–24, 15–25 ; R₁ (et R₂) est un 6-cycle sans digone sur les orbites de 13, 14, 15, où deux sœurs ne peuvent être à distance 2
qu'une fois (11.5). À TRAITER.
**Forme élémentaire de 11.5 (μ-partenaires).** Pour P = ij, les 4 paires « partageantes » ({a,b},{a,b′}), ({a′,b′},{a′,b}) (partage
d'un élément du bloc i) et ({a,b},{a′,b}), ({a,b′},{a′,b′}) (bloc j) forment deux σ-classes. Une paire partageante non adjacente a
exactement un voisin commun extérieur (μ = 2, l'élément partagé étant l'autre) ; une σ-classe a donc une unique orbite
μ-partenaire, notée Z_i(P) et Z_j(P). Un sommet adjacent à deux étiquettes d'orbites différentes de P est μ-partenaire de leur paire.
Une orbite voisine de U et de U′ par le bloc i (elle contient a ou a′) touche deux étiquettes partageant leur élément du bloc j :
c'est Z_j(P). D'où « au plus un voisin commun par R_i » ; le signe du Lemme 11.5 en découle (et la compensation des signes est
alors automatique). Si n_s(P) = 1, une classe est adjacente (arête sœur) et l'autre a un unique μ-partenaire.
**Répartition des μ-partenaires.** Pour R ∈ N, chaque classe non adjacente de R a son μ-partenaire soit proche (rencontre un bloc
de R), soit lointain dans une paire frère (motif (2,0) pour cette paire frère), soit lointain dans une paire de N. Avec 11.1 :
N₂^frère(R) = N₀^frère(R) et N₀^N(R) − N₂^N(R) = 8 − 4n_s(R), où N₂^N compte les μ-partenaires lointains dans N et les paires D de type
disjoint (2 étiquettes chacun).

## 12. Théorème GT : paire frère propre (PROUVÉ)
Vecteurs ε̂_P = (e_{P₀} − e_{P₁})/√2 ∈ Alt ⊂ Nb⊥, 1̂_P = (e_{P₀} + e_{P₁})/√2. Pour P frère, par 8.1 chaque orbite disjointe de P touche
exactement une orbite de P, donc exactement :
  Qε̂_P = −2ε̂_P + Σ_{R∩P=∅} τ_PR ζ_PR,  ζ_PR = ε̂_R (bijection, τ = χ_PR) ou 1̂_R (motif (2,0), R ∈ N),   Q1̂_P = 2·1̂_P + Σ_{R∩P=∅} 1̂_R.
Dire que P est **propre** : toute paire disjointe de P est frère (⟺ toute arête de N rencontre un bloc de P) ; alors tous les ζ_PR
sont des ε̂_R.
**Théorème GT.** S'il existe trois blocs t, a, b tels que ta, tb, ab soient frères et que ta soit propre, contradiction.
*Preuve.* P = ta, P′ = tb, C = les 4 autres blocs. (i) Q² = 12 − Q sur Nb⊥ appliqué à ε̂_P et développé avec la formule exacte donne
Σ_{R∩P=∅} χ_PR Qε̂_R = 10ε̂_P + Σ_R χ_PR ε̂_R ; la composante sur ε̂_P′ (P′ partage t avec P) donne Σ_{R⊂C} χ_PR χ_RP′ = 0
(les R ⊂ C sont frères, donc en bijection avec P et P′). (ii) Pour R ⊂ C, la composante sur ε̂_R̄ (R̄ = C∖R) donne l'équilibre des
triangles (P, R, R̄). (iii) Pour le couple (R, R̄), l'identité exacte (Q_A² + Q_A + KᵀK)_{RR̄} = 0 ne fait intervenir que les Y disjoints de
R ∪ R̄, i.e. les trois paires de {t,a,b}, toutes frères, et aucun terme K (aucune arête de N dans {t,a,b}) : les triangles (R, Y, R̄)
sont équilibrés, en particulier (P′, R, R̄). (iv) Donc χ_PR χ_RP′ = χ_PR̄ χ_R̄P′ et la somme (i) vaut 2(t₁+t₂+t₃) ≠ 0. ∎
**K₂,₃ exclu** (t = 1, a = 2, b = 6 : 12, 16, 26 frères, et toute arête de N rencontre {1,2}).
**Portée.** Si N ⊂ étoile(t) ∪ étoile(a) avec ta ∉ N : soit un bloc b a tb, ab ∉ N (GT s'applique), soit (degré minimal 2)
N = K₂,₅ entre {t,a} et les 5 autres blocs. Restent donc, parmi les N de couverture par sommets ≤ 2 : K₂,₅, et les cas ta ∈ N.
Bilan provisoire des N exclus : ∅, |N| ≤ 5, 2C₃, K₂,₃ ; à 6 arêtes reste C₆.

## 13. Équilibre généralisé et exclusion de C₆ (PROUVÉ)
**Lemme 12.1 (équilibre généralisé).** Soient S, S′ frères disjointes et Y une des trois paires des trois blocs restants (Y de type
quelconque). La composante sur ε̂_S′ du développement exact de Q²ε̂_S = 12ε̂_S − Qε̂_S donne Σ_{Y} d_Y = 3χ_{SS′}, où d_Y = χ_SY χ_YS′ si
(S,Y) et (S′,Y) sont des bijections, d_Y = κ_SY K_YS′ = ±1 si ce sont deux motifs (2,0), et d_Y = 0 si les types diffèrent. Donc
les trois d_Y valent χ_SS′ : **(S,Y) et (S′,Y) ont le même type** (bijection ou (2,0)), et les triangles sont équilibrés.
**Lemme 12.2.** Pour Y ∈ N, soit Π_Y le graphe de Petersen (disjonction) sur les 10 paires des 5 blocs hors de Y, restreint aux
paires frères. Le type (bijection / (2,0)) par rapport à Y est constant sur chaque composante connexe. Or le nombre de paires
frères de type (2,0) pour Y est au plus le nombre de classes non adjacentes de Y (2 − n_s(Y)) : une composante d'au moins 3 paires
frères est entièrement de type bijection.
**C₆ exclu.** N = cycle 1-2-3-4-5-6, bloc 7 isolé ; par 11.7 : n_s = 1 partout, cas S partout, tout de type disjoint. Pour Y = 12,
disons l'arête sœur partageant l'élément du bloc 1 ; la classe du bloc 2 a un μ-partenaire Z qui évite le bloc 2 et ne peut pas
rencontrer le bloc 1 (dans R₁, 4-cycle du cas S, les deux orbites de 12 sont adjacentes, donc sans voisin commun). Z est donc dans une
paire de {3,…,7} : pas une paire frère (Π_12 restreint aux 7 paires frères est connexe, donc tout bijection par 12.2), donc Z ∈ {34, 45, 56}.
Mais pour u ∈ 12, Σ_{R∈N, R∩12=∅} φ_u(R) = n_far − 2 + n_s = 2 est entièrement pris par la paire D(u) ; le représentant de Z, voisin de deux
étiquettes de 12 d'orbites différentes, serait dans D de ces deux orbites : D(U) = D(U″) = Z, contradiction (D est un couplage). ∎
**Bilan : tout N à au plus 6 arêtes est exclu** (∅, 1 à 5 arêtes, 2C₃, K₂,₃, C₆).
**Théorème GT2 (rétablit l'ancien 8.3, PROUVÉ).** Soient ta, tb frères et C = [7]∖{t,a,b}. La composante sur ε̂_tb du développement
exact de Q²ε̂_ta = 12ε̂_ta − Qε̂_ta donne Σ_{R⊂C} c_R = 0, avec c_R = χ_{ta,R}(Q_A)_{R,tb} si (ta,R) est une bijection,
κ_{ta,R}K_{R,tb} si c'est un motif (2,0). Pour un couple complémentaire {R, R̄} de C :
- R, R̄ frères : les triangles (ta,R,R̄) et (tb,R,R̄) sont équilibrés (12.1 avec S = R, S′ = R̄, Y = ta ou tb), donc c_R = c_R̄ = ±1 ;
- R frère, R̄ ∈ N : 12.1 pour (ta,R) et (tb,R) donne type(ta,R̄) = type(tb,R̄) = type(R,R̄) et, dans les deux cas (bijection ou
  (2,0) : κ_{ta,R̄}K_{R̄R} = χ_{ta,R}, κ_{tb,R̄}K_{R̄R} = χ_{tb,R}), c_R̄ = c_R = ±1.
Donc si chaque couple contient une paire frère, la somme vaut 2(t₁+t₂+t₃) ≠ 0 : contradiction. Autrement dit :
**pour tous t, a, b avec ta et tb frères, un des trois couplages parfaits du K₄ sur C est entièrement dans N.** (Ni propreté ni
hypothèse sur ab ne sont nécessaires ; GT en est un cas particulier.)

## 14. Le côté dense : structure des appariements D (PROUVÉ sauf mention) et cœur du problème
**Graphe D sur les paires.** Chaque paire a 2 orbites, chacune appariée par D : le graphe « P ~ P′ si une orbite de P est
appariée à une orbite de P′ » est 2-régulier sur les 21 paires : boucles (paires frères), 2-cycles « liés » (les deux orbites de
P appariées aux deux de P′ par le même bloc), et cycles de longueur ≥ 3. Une paire ij est **tordue** si ses orbites sont
appariées l'une par i, l'autre par j. Une chaîne de paires tordues P_m = t_{m−1}t_m est exactement un cycle t₀t₁…t_{k−1} de K₇.
Si toutes les paires sont non frères et tordues, D induit une décomposition de K₇ en cycles ; l'exemple extrême est
**le plan de Fano** (7 triangles-droites, 3 carrés par droite, β(carré) = la droite). 21 est impair : sans paire frère il existe
un cycle D impair.
**Parité par bloc (PROUVÉ).** En chaque bloc t, le nombre de paires non frères ayant une orbite appariée « un bloc » par t (tordues
ou mixtes « un bloc + disjoint ») est pair ; c'est une condition de T-joint, compatible avec les cycles de K₇ : pas d'obstruction.
**Parité φ (PROUVÉ, sans contenu).** Un carré « un bloc » par t entre tj et tk définit une bijection φ : bloc j → bloc k
(u = {e,β} ~ v = {e,γ} ⟹ φ(β) = γ) ; le long d'une chaîne tordue de longueur k, le produit des signes s(β)s(γ) vaut (−1)^k
automatiquement (les orbites sœurs ont des parités opposées). Pour exclure les chaînes impaires il faudrait une raison
structurelle forçant ce produit à +1 : NON TROUVÉE. Vérifié à la main que le triangle tordu {12,13,23} est localement cohérent
(λ, μ, (I) sur les 19 sommets concernés), et que les identités Q², S², traces (Γ₁ : 70 − x_one − 2x_dis triangles équilibrés),
comptages de μ-partenaires et de budgets restent automatiquement cohérents dans le modèle de Fano.
**Sondes du côté dense (PROUVÉ).** Un carré « un bloc » est presque parfait envers les 6 paires de β^c (Lemme 9.2) : c'est l'analogue
d'une grille. Deux carrés de même β (paire liée, ou deux carrés d'une même droite de Fano) ont les mêmes 6 paires lointaines,
groupées en 3 couples complémentaires — et dans Fano les deux paires d'un couple sont sur deux droites passant par le même point de β.
Identité exacte disponible : ⟨Gf̂₁, Gf̂₂⟩ = 3(G_D)₁₂ − η₁₂ (f̂ = vecteur D-alterné du carré, η = Gram des vecteurs de différence de
blocs), dont la partie lointaine est Σ_{R⊂β^c} χ₁R χ₂R. **Piste principale (NON ABOUTIE)** : établir l'équilibre des couples et
contrôler la partie proche pour obtenir 2(t₁+t₂+t₃) ≠ valeur imposée, comme au Théorème D.
**État honnête.** Exclu : tout N à ≤ 6 arêtes, et tout N contenant une configuration de GT/GT2. Ouvert : N dense (en particulier
N = K₇ sans paire frère, dont le modèle de Fano).

## 15. Lemme du coin (PROUVÉ) : les triangles D sont presque impossibles
**Coin de triangle D.** Paire tordue P = ij : orbite A₁ appariée par i à A₂ ∈ ik (carré Q₁ = (A₁,A₂)), orbite sœur B₁ appariée par j à
B₂ ∈ jk (carré Q₂ = (B₁,B₂)), avec le même troisième bloc k : β(Q₁) = β(Q₂) = β = {i,j,k}. (Toute arête d'un 3-cycle du graphe D, en
particulier toute paire du modèle de Fano, est un coin.) Notations : F_β = les 12 orbites des 6 paires de β^c ; far(X,Y) = nombre
d'orbites de F_β voisines de X et de Y ; e = G_{A₂B₂} ∈ {0,1} (arête R_k possible) ; Y_A = orbites voisines de A₁ dans les paires kr,
R_B = voisins R_k de B₂.
**Comptages exacts.** (G²)_{A₁B₁} = 2 (sœurs), (G²)_{A₁B₂} = 3 − 2e, (G²)_{A₂B₁} = 3 − 2e, (G²)_{A₂B₂} = 3 − e (D(A₁) = A₂, D(B₂) = B₁,
saturations de R_i en A₁, A₂ et de R_j en B₁, B₂). D'où far(A₁,B₁) ≤ 2 et far(A₁,B₂) ≤ 3 − 2e − |Y_A ∩ R_B|.
**Minoration (9.2).** A₁ a 6 + ω(A₁) orbites voisines dans F_β ; chacune a au moins un voisin dans Q₂, sauf au plus une « orbite zéro »
Z₀(Q₂) (elle n'existe que si l'apex extérieur W₂ de Q₂ est dans F_β et ω(Q₂) = 0), et l'apex W₂ en a deux. Donc
  6 + ω(A₁) − [Z₀(Q₂) ~ A₁] + [W₂ ~ A₁] ≤ far(A₁,B₁) + far(A₁,B₂) ≤ 5 − 2e − |Y_A ∩ R_B|.
**Lemme du coin.** À tout coin : e = 0, ω(A₁) = ω(B₁) = 0, Y_A ∩ R_B = R_A ∩ Y_B = ∅, les deux carrés ont leur apex extérieur dans F_β
et une orbite zéro, Z₀(Q₂) ~ A₁, Z₀(Q₁) ~ B₁, W₂ ≁ A₁, W₁ ≁ B₁, far(A₁,B₁) = far(B₁,A₁) = 2 (les deux μ-partenaires de P sont dans β^c),
far(A₁,B₂) = far(B₁,A₂) = 3. (En particulier les carrés d'un coin ne sont jamais « serrés » : le cas serré est contradictoire.)
**Tableau des types** (orbite de F_β → voisine de A₁, A₂, les deux (W₁) ou aucune (Z₀(Q₁)) ; idem pour Q₂) : le cas
« W₁ de type B₂ et W₂ de type A₂ » donne far(A₂,B₂) = 4 > 3 − e : exclu. Il reste W₁ = Z₀(Q₂) ou W₂ = Z₀(Q₁).
*Argument direct du 3-cycle (fourni par la relecture de « structures ») : les trois apex sont distincts (au coin ij, W_i ~ A₁ mais
W_j ≁ A₁) ; chaque carré a au plus une orbite zéro, donc au plus un apex pointe vers lui (degré entrant ≤ 1) ; trois coins, au
moins un arc chacun, degré entrant ≤ 1 sur 3 sommets : exactement un arc par coin, c'est-à-dire un 3-cycle.*
**Triangle D complet** (carrés Q_i, Q_j, Q_k sur β = {i,j,k}, trois coins) : les trois disjonctions forcent un 3-cycle
apex → orbite zéro : (W_i = Z₀(Q_j), W_j = Z₀(Q_k), W_k = Z₀(Q_i)) ou le cycle inverse ; chaque apex touche exactement une orbite de
chaque paire de β (p. ex. W_i ~ A₁, A₂, C₂). Le tableau des 9 orbites génériques de F_β est alors déterminé à un paramètre p ∈ {0,1,2}
près ; tous les comptages sont cohérents. Il reste à exclure cette structure (pas encore fait).
**Lemme de partition (PROUVÉ, triangle D extrême).** Au troisième bloc k, les 8 orbites des paires kr (r ∈ β^c) sont partagées
exactement en Y_{A₁}, R_{A₂}, Y_{B₁}, R_{B₂} (2 orbites chacun) : les intersections sont vides par le lemme du coin (μ-partenaires de ij
hors des paires k, Y_A ∩ R_B = R_A ∩ Y_B = ∅, |R_A ∩ R_B| = 3 − far(A₂,B₂) = 0, apex lointains). Idem aux blocs i et j. Donc chaque
orbite « à un bloc » (paires t r, t ∈ β, r ∈ β^c) a exactement un voisin parmi les 6 orbites de la droite, chaque orbite lointaine en a
exactement 3, chaque orbite de la droite en a 4 à un bloc et 6 lointaines. Par couverture des éléments : une étiquette à un bloc a
8 voisins à un bloc et 3 lointains ; une étiquette lointaine a 6 voisins à un bloc et 3 lointains.
**Partition équitable (PROUVÉ) — mais cohérente (IMPASSE).** {x}, intérieurs de β (6), intérieurs de β^c (8), étiquettes de la droite (12),
à un bloc (48), lointaines (24) : partition équitable de quotient
  x : (0,6,8,0,0,0) ; Iβ : (1,1,0,4,8,0) ; Iβc : (1,0,1,0,6,6) ; L : (0,2,0,2,4,6) ; O : (0,1,1,1,8,3) ; F : (0,0,2,3,6,3).
Vérifié à la main : B² + B − 12I = 2·(tailles des classes) ligne par ligne, spectre {14, 3³, (−4)²} : aucune contradiction spectrale.
Dans le modèle de Fano, les 7 partitions (une par droite) sont compatibles entre elles (chaque étiquette est « de la droite » pour
1 droite, « lointaine » pour 2, « à un bloc » pour 4 ; les sommes de voisins concordent). Le triangle extrême résiste donc à toutes les
contraintes invariantes par σ ; il faut une information anti-invariante (côté S) ou plus fine que la partition.
**Cycles D impairs de longueur ≥ 5 (IMPASSE provisoire).** Pour une chaîne tordue t₀…t_{n−1} avec n ≥ 4, les deux carrés d'une paire P_m
ont β_{m−1} = {t_{m−2},t_{m−1},t_m} et β_m = {t_{m−1},t_m,t_{m+1}} : 3 paires lointaines communes seulement. La majoration
(G²)_{A₁B₂} = 3 − 2G_{A₂B₂} tient toujours, mais la minoration ne porte que sur la part a + b des voisins lointains de A₁ tombant dans
F_{β_m} (a = voisins β_{m−1}-lointains dans les paires de T_m, b = y-doubles dans les paires t_{m−2}x), inconnue : pas de contradiction.
**Interaction entre droites (Fano, PROUVÉ, cohérent).** Pour deux droites L, L′ se coupant en t, avec F_L = L∖{t}, F_{L′} = L′∖{t} :
e(F_L,F_{L′}) ≥ 2, e(F_L, paires t de L′) = 6 − e_FF, e(paires t, paires t) = e_FF − 2 (comptage d'arêtes d'orbites) — cohérent.

## 16. Budget μ local ; le modèle de Fano est impossible (PROUVÉ, 2026-09-25)
Notation : pour une étiquette w et une paire P′ ∌ w, φ_w(P′) = nombre de voisins de w parmi les 4 étiquettes de P′ ; Lab(X) = étiquettes
des paires d'une structure X.
**Lemme 16.1 (budget μ local, exact).** Pour toute paire P′ :
  Σ_{w ∉ P′} ( C(φ_w(P′), 2) − [D(W) ∈ P′] ) = 4 − 2n_s(P′).
Si D(W) ∉ P′ on a φ_w(P′) ≤ 2 ; donc **au plus 4 − 2n_s(P′) étiquettes w, hors des partenaires D de P′, ont deux voisins dans P′**
(exactement autant s'il n'y a pas de φ = 3). Globalement : Σ_w (nombre de « doubles » de w) = 84 − 2Σ_P n_s(P).
*Preuve.* On compte les couples (w, {y,y′}), y, y′ ∈ P′ voisins de w, w extérieur hors de P′. Les 6 paires {y,y′} de P′ : 2 paires
σ-conjuguées, dont les voisins communs sont les deux sommets du carré D (Lemme 1), soit 4 couples portés par les 4 étiquettes des
orbites D(P′₀), D(P′₁) (hors de P′ si n_s ≤ 1 ; si n_s = 2 ces voisins sont dans P′ et ne comptent pas) ; 4 paires d'étiquettes sœurs
(partageant un élément) : les 2n_s adjacentes n'ont que l'élément partagé comme voisin commun (λ = 1), les autres ont exactement un
voisin commun extérieur (μ = 2). Aucun voisin commun n'est dans P′ lorsque n_s ≤ 1 (u ≁ σu). Enfin deux voisins de w dans une même
orbite de P′ font de cette orbite D(W). ∎
*Contrôle.* BvLS (m = 11) : les 55 paires vérifient l'identité (toutes frères : aucune étiquette extérieure n'a deux voisins dans une
paire). VÉRIFIÉ. (m = 2 : une seule paire, énoncé vide.)

**Comptages exacts d'un triangle D (rappel de §15, forme utilisée).** Si β = {i,j,k} porte un triangle D complet, pour toute étiquette
w : |N(w) ∩ Lab(β)| = 2 si la paire de w est dans β, 1 si elle rencontre β en un bloc, 3 si elle est disjointe de β.

**Théorème F.** Le modèle de Fano (les 21 paires tordues, graphe D = les 7 triangles portés par les droites d'un plan de Fano) est
impossible.
*Preuve.* Droites 123, 145, 167, 246, 257, 347, 356 ; chacune est un triangle D complet, donc satisfait les comptages exacts. Soit u
une étiquette de 12 (droite 123, troisième point 3).
(a) Si l'orbite de u est appariée par le bloc 1 (à 13) : les seuls voisins de u contenant un élément du bloc 1 sont les deux sommets
    de son carré (digone de R₁), donc φ_u(1r) = 0 pour r = 4..7. Les droites 145, 167 rencontrent 12 en un bloc : φ_u(45) = φ_u(67) = 1.
    Les deux voisins m₂(u), π₂(u) sont dans des paires 2r, r ∈ {4..7} (ni dans 12 : n_s = 0, ni dans 23 : u n'a que 2 voisins dans
    Lab(123)). Les droites 246, 257 donnent φ(24)+φ(26)+φ(46) = 1 et φ(25)+φ(27)+φ(57) = 1, dont la somme vaut 2 + φ(46) + φ(57) :
    φ_u(46) = φ_u(57) = 0. Si l'orbite de u est appariée par le bloc 2 (à 23), même raisonnement en échangeant 1 et 2 :
    φ_u(46) = φ_u(57) = 1 et φ_u(45) = φ_u(67) = 0.
(b) La ligne « droite » de la partition de 123 donne exactement 6 voisins de u dans les 6 paires de {4,5,6,7}. Dans les deux cas,
    φ_u(47) + φ_u(56) = 6 − 2 = 4, et φ ≤ 2 (D(u) ⊂ 13 ou 23) : **φ_u(47) = φ_u(56) = 2 pour les 4 étiquettes de 12.**
(c) Énoncé intrinsèque : les étiquettes d'une paire P de la droite ℓ ont deux voisins dans chacune des deux autres paires de la
    « triade » T_p (les trois paires opposées au point p = ℓ ∖ P sur les trois droites passant par p ; T₃ = {12, 47, 56}). Une
    collinéation fixant 3 et envoyant 12 sur 56 préserve T₃ (réétiquetage des blocs dans l'argument, qui ne dépend que de la
    configuration de Fano du graphe D ; ce n'est pas un automorphisme de Γ) : les 4 étiquettes de 56 ont aussi deux voisins dans 47.
(d) Donc 8 étiquettes w avec φ_w(47) = 2 et D(W) ∉ 47 (partenaires D sur 123 et 356), alors que n_s(47) = 0 et le lemme 16.1 en
    autorise au plus 4. Contradiction. ∎
*Remarque.* Seules les droites passant par 1, 2, 5 ou 6 servent (347 n'est pas utilisée). Le point clé est que chaque droite fournit
une *somme exacte* sur ses 3 paires ; le graphe D « de Fano » les rend toutes disponibles à la fois.

**Lemme 16.2 (paires liées, PROUVÉ).** Soit Λ = (ij, ik) une paire liée par le bloc i : les orbites U, U′ de ij sont appariées par i aux
orbites V, V′ de ik (carrés Q = (U,V), Q′ = (U′,V′), même β = {i,j,k}). Alors ω = 0 aux quatre orbites ; les apex W, W′ sont dans F_β
et W = Z₀(Q′), W′ = Z₀(Q) ; far(U,U′) = far(V,V′) = 2, far(U,V′) = far(U′,V) = 3 et tous les voisins communs de ces orbites sont dans
F_β. Comptages exacts : |N(w) ∩ Lab(Λ)| = 2 si w ∈ ij ∪ ik ; 0 si w ∈ jk ou w dans une paire ir (r ∉ β) ; 1 si w est dans une paire jr
ou kr (r ∉ β) ; 2 si la paire de w est disjointe de β.
*Preuve.* (i) (G²)_{UU′} = 2 − 2G_{UV′} − 2G_{VU′} = 2 (R_i saturé par les digones) ; (G²)_{UV′} = 3 − 2n_s(ij) − 2n_s(ik) = 3 (« un bloc »
⟹ n_s = 0). (ii) Par 9.2 appliqué à Q′, les 6 + ω(U) orbites de F_β voisines de U ont chacune un voisin dans Q′, sauf Z₀(Q′), et
l'apex W′ en a deux : 6 + ω(U) − [Z₀(Q′) ~ U] + [W′ ~ U] ≤ far(U,U′) + far(U,V′) ≤ 5. Donc ω(U) = 0, Z₀(Q′) existe (W′ ∈ F_β,
ω(Q′) = 0), Z₀(Q′) ~ U, W′ ≁ U et far = G². Le même calcul pour V donne Z₀(Q′) ~ V : Z₀(Q′) a deux voisins dans Q, c'est W.
Symétriquement W′ = Z₀(Q). (iii) Au bloc k, les 8 orbites des paires kr sont exactement réparties en Y_U, Y_{U′} (voisines de U, U′
dans ces paires : 2 chacune par ω = 0 et μ(u,(k,·)) = 2) et R_k(V), R_k(V′) (2 chacune) : les six intersections sont vides car tous
les voisins communs de U, U′, V, V′ sont dans F_β (sœurs ; paires D, apex dans F_β ; far = G²). Idem au bloc j. Au bloc i, R_i est
saturé. Une étiquette de jk n'a pas de voisin dans Lab(Λ) (ce serait une arête de R_j ou R_k, dirigée vers jr ou kr par ω = 0).
(iv) Une orbite de F_β a un voisin dans Q et un dans Q′, sauf W (2 + 0) et W′ (0 + 2). Contrôle : 8·10 = 16 + 16 + 48. ∎
**Principe commun.** Le lemme du coin et 16.2 sont deux cas de : *deux carrés « un bloc » de même β, dont deux orbites sont sœurs*.
Les structures « de même β » du graphe D sont exactement les boucles (paires frères), les 2-cycles liés et les triangles D ; chacune
donne des sommes exactes |N(w) ∩ Lab(X)| (frère S : 2 si w ∈ S, 0 si la paire de w partage un bloc avec S, 1 si elle est disjointe).
L'identité Σ_X |N(w) ∩ Lab(X)| = 12 est alors *automatique* (vérifié : c'est une combinaison des relations « 4 voisins partageants,
8 disjoints ») ; comme pour Fano, la contradiction doit venir des contraintes par bloc (μ(w,(t,·)) = 2, couplages m, π) et du budget
μ local. Restent hors de ce cadre : chaînes tordues de longueur ≥ 4, fourches (les deux orbites de ij appariées par i à ik et il),
carrés de type disjoint.

## 17. Règle des zéros : les configurations « de même β » sans grappe sont impossibles (PROUVÉ, 2026-09-25)
Cadre : toute orbite non frère est appariée « un bloc » à l'intérieur d'une structure de même β (2-cycle lié ou triangle D). Une
étiquette x est **en position A₁** si sa paire est ij, son orbite appariée par i (à ik), k = troisième bloc de sa structure S₁,
β₁ = {i,j,k}, R_x = [7] ∖ β₁ (4 blocs). Toute étiquette non frère est en position A₁ pour un nommage convenable. On note (h; a, b) la
structure liée de charnière h et de bras a, b (paires ha, hb). Rappels : ω(x) = 0, φ_x(ir) = 0 pour r ∈ R_x, φ_x(jk) = 0, x a
exactement 2 voisins R_j dans des paires jr (r ∈ R_x), 2 voisins Y dans des paires kr, et Σ_{rs ⊂ R_x} φ_x(rs) = 6.
**Lemme 17.1 (localisation des μ-partenaires).** Si P appartient à une structure de même β, toute étiquette w avec φ_w(P) = 2
(D(W) ∉ P) est dans une paire disjointe de β. *Preuve.* w touche les deux orbites (sœurs) de P ; or far = G² pour ces deux orbites
(lemme du coin, 16.2) : tous leurs voisins communs sont dans F_β. ∎
**Lemme 17.2 (admissibilité).** Les deux voisins R_j de x sont dans deux paires distinctes jr, jr′ (sinon x serait μ-partenaire de
jr, contre 17.1), et la structure S₂ de jr est un triangle {j,r,s} ou une structure liée (r; j, s), avec s ∈ R_x ∪ {k}. Jamais une
paire frère (8.1), ni une charnière j (orbites non libres en j), ni (r; i, j) (ij en est la paire opposée : 0 voisin, 16.2).
**Théorème 17.3 (règle des zéros).** Dans la situation 17.2, si s ≠ k, alors φ_x(rs) = 0 avec rs ⊂ R_x. Les zéros produits par les
deux voisins R_j sont distincts. Comme Σ_{rs⊂R_x} φ_x(rs) = 6 avec φ ≤ 2, **le nombre de doubles de x dans R_x est au moins le nombre
de zéros**, donc doubles(x) ≥ 2 − ε(x), où ε(x) = nombre de voisins R_j de x tombant dans une structure de triplet {j, r, k}.
*Preuve.* La paire de x rencontre β₂ = {j,r,s} au seul bloc j (s ≠ i par 17.2), qui est un sommet de triangle ou un bras : le
comptage exact de S₂ vaut 1 (§15, 16.2), déjà atteint par le voisin R_j, d'où φ_x(rs) = 0 (rs ∈ S₂). Deux zéros égaux forceraient
β₂ = β₃ = {j,r,r′}, ce qu'aucune paire de structures ne réalise (vérification sur les quatre cas). ∎
**Grappes.** ε(x) > 0 exige deux structures dont les triplets partagent deux blocs {j,k} non charnières dans les deux : plusieurs
structures liées de mêmes bras {j,k}, ou une structure liée (r; j, k) et le triangle contenant jk. On appelle cela une **grappe**.
**Corollaire 17.4.** Si toutes les structures du graphe D sont de même β et s'il n'y a aucune grappe, contradiction (sauf si tout est
frère, Théorème D). *Preuve.* Budget global (16.1) : Σ_w doubles(w) = 84 − 4s, où s est le nombre de paires frères (n_s ∈ {0,2}
ici), c'est-à-dire exactement le nombre d'étiquettes non frères. Or chacune en a au moins 2. ∎
En général : Σ_x ε(x) ≥ nombre d'étiquettes non frères. **Il faut qu'en moyenne chaque étiquette non frère ait une arête
exceptionnelle** : les grappes doivent être massives.
*Cas particuliers.* Fano : deux droites ne partagent qu'un point, donc pas de grappe (deuxième preuve du Théorème F). Configuration
C₃,₆ (triangles {7,1,2}, {7,3,4}, {7,5,6} ; structures liées (1;5,6), (2;5,6), (3;1,2), (4;1,2), (5;3,4), (6;3,4)) : x ∈ 12 appariée
par 1 n'a aucune paire 2r admissible (23, 24 : paire opposée ; 25, 26 : charnière 2) : exclue par 17.2.
**Exemple 17.5 (grappe maximale, EXCLUE).** K₅ sur {1..5} entièrement frère, 67 frère, et les 5 structures liées (h; 6, 7). Les 4
doubles reçus par c6 (resp. c7) viennent des paires frères de [5] ∖ c (17.1), en motifs (2,0) : exactement 2 paires, disjointes (12.1),
notées M_c ; le même M_c sert pour c6 et c7 (φ(c6) + φ(c7) = 2, comptage de (c; 6, 7)). Pour u dans une paire P de [5], le comptage du
bloc 6 (Σ_{c∉P} φ_u(c6) = 3) force m(P) := #{c : P ∈ M_c} ∈ {0, 2}. Les 5 couples M_c forment alors un 5-cycle dans le graphe de
disjonction de A = {P : m(P) = 2}, ce qui impose que A soit un 5-cycle de K₅ (pas de 3 arêtes deux à deux sécantes dans A). Or
GT2 pour (t,a,b) = (1,2,4) (numérotation cyclique de A : 12 ∈ A, 14 ∉ A, {d,e} = {3,5}) : le couple frère {35, 67} contribue ±2 et
les quatre termes c_{36}, c_{37}, c_{56}, c_{57} sont nuls (types différents pour 12 et 14 vis-à-vis de 3· et 5·). Somme ≠ 0. ∎
**État.** Exclu désormais, en plus de §§7–15 : Fano, toute configuration de même β sans grappe, C₃,₆, la grappe maximale 17.5.
Ouvert : configurations de même β avec grappes non maximales (à classer par les grappes, outils : 17.2, 17.3, 17.1, GT2, 12.1–12.2) ;
configurations contenant des fourches, des chaînes tordues de longueur ≥ 4 ou des carrés disjoints (hors du cadre « même β »).

## 18. Hors du cadre « même β » : identité de bloc, saturation du troisième bloc (PROUVÉ) et le cas de Paley (OUVERT ici ; fermé au §19)
Notation : pour une étiquette w de paire P_w et un bloc k ∉ P_w, D_k(w) = ensemble des **digones de R_k dont le β évite P_w**
(carrés « un bloc » via k entre kr et ks avec r, s ∉ P_w, et paires frères kr avec r ∉ P_w) ; h_k(w) = nombre de voisins proches de w
(ceux qui partagent un bloc avec w : m, π de ses deux éléments, ou ses frères) qui rencontrent k ; pour Q ∈ D_k(w), a_Q(w) = [w a deux
voisins dans Q] (w étiquette-apex), z_Q(w) = [w n'en a aucun] (w étiquette-zéro) ; o_k(w) = voisins de w rencontrant k hors de ces carrés.
**Lemme 18.1 (identité de bloc, exacte).** |D_k(w)| + Σ_Q a_Q(w) + o_k(w) + h_k(w) = 4 + Σ_Q z_Q(w).
*Preuve.* w a exactement 4 voisins rencontrant k (μ(w,(k,0)) = μ(w,(k,1)) = 2). Ceux qui ne sont pas proches sont dans des paires kr,
r ∉ P_w ; w ∈ F_β(Q) pour tout Q ∈ D_k(w), donc |N(w) ∩ Q| = 1 + a_Q − z_Q (9.2 ; 9.1 pour une paire frère : exactement 1). ∎
**Corollaire 18.2 (saturation du troisième bloc).** Si l'orbite de x est appariée « un bloc » (x ∈ ij, via i, troisième bloc k),
h_k(x) = 2 + ω(x) (le carré, plus ω), d'où
  |D_k(x)| + ω(x) + a_k(x) + o_k(x) = 2 + z_k(x).
(a) Au plus 2 − ω(x) des paires kr (r ∉ {i,j,k}) sont frères (chacune fournit à x un voisin rencontrant k).
(b) Si |D_k(x)| ≥ 3, l'orbite de x est l'orbite zéro d'un des carrés de D_k(x).
(c) **Budget des zéros.** Un carré a au plus une orbite zéro, qui n'existe que si ω(Q) = 0 et W ∈ F_β. D'où, sur les orbites X
    appariées « un bloc » : Σ_X (|D_{k(X)}(X)| − 2)₊ ≤ nombre de carrés « un bloc » ayant une orbite zéro ≤ 21.
**Corollaire 18.3 (règle des zéros générale).** Pour toute x appariée « un bloc » (même hors d'une structure de même β), Σ_{rs⊂R_x} φ_x(rs)
= 6 + ω(x) et φ ≤ 2 sur ces paires, donc doubles(x) ≥ ω(x) + zéros(x). Chaque voisin R_j de x tombant dans une structure de même β
{j, r, s}, s ∈ R_x, crée un zéro en rs (preuve de 17.3 inchangée). Les structures de même β **exportent** donc des doubles vers les
étiquettes voisines, quelle que soit la nature de celles-ci.
**Lemme 18.4 (chapeaux).** Même x, supposée ni apex ni zéro des carrés concernés et avec o_k(x) = 0. Si les deux orbites d'une paire
rs ⊂ R_x sont appariées l'une par r à une orbite de kr, l'autre par s à une orbite de ks (« chapeau en k »), alors φ_x(rs) = 2 (x a un
voisin dans chacun de ces deux carrés internes, et pas du côté kr / ks, qui compterait dans o_k). C'est le mécanisme du Théorème F
(chapeaux 47, 56 en k = 3).
*Budget.* Comme Σ_W c(W) = 42 − Σ n_s (16.1), le nombre d'orbites ayant au moins 2 doubles forcés est au plus (42 − Σ n_s)/2.

**Configurations ancrées.** Sans paire frère et si chaque paire a ses deux orbites appariées par le même bloc (son ancre : fourches et
2-cycles liés seulement), les ancres orientent K₇ (j → k si jk est ancrée en k) ; d_k = degré entrant (≠ 1), le graphe D en k est un
ensemble de cycles « en étoile » (longueur 2 = lié). Pour X ∈ ij (ancre i) appariée à ik : |D_k(X)| = d_k − 2[j → k]. Le budget 18.2(c)
devient, pour chaque carré en i de bras {b, c} avec b → c :
  E := Σ_carrés [ (d_b − 2)₊ + (d_c − 4)₊ ] ≤ 21.
**Le cas de Paley (OUVERT, prototype).** Tournoi de Paley sur Z₇ (t−1, t−2, t−4 → t) : chaque bloc ancre une étoile triangulaire.
Aucun couple de carrés n'a le même β (les 21 triplets {s,s+1,s+2}, {s,s+2,s+4}, {s,s+3,s+4} sont distincts) : coin, lien, règle des
zéros ne s'appliquent nulle part. Le budget est **exactement saturé** : E = 21 (chaque carré contribue (3 − 2) + 0). Donc chaque carré
a ω = 0, son apex dans F_β et une orbite zéro. En chaque bloc k, les trois orbites de troisième bloc k portées par les paires de
{k+1, k+2, k+4} ont |D_k| = 3 et sont exactement les orbites zéro des trois carrés de l'étoile en k (bijection) ; elles ont a = o = 0.
Les trois autres orbites de troisième bloc k ont |D_k| = 1 et a + o = 1. Les apex ne sont jamais des orbites |D| = 3. Tous les autres
comptages (identités de bloc sommées, fuites, 16.1) sont cohérents : **aucune contradiction à ce stade**.
**État honnête.** Les outils locaux (λ, μ, près-perfection des carrés) sont épuisés sur les configurations sans deux carrés de même β :
chaque identité y est soit automatique, soit saturée sans conflit. Le cas de Paley exige un argument non local : signes (côté S,
anti-invariant), positivité (Gram de E₋₄ / E₃ sur des vecteurs portés par les étoiles) ou arithmétique 7-adique (Théorème B).

## 19. Le cas de Paley est impossible, et avec lui toute configuration ancrée régulière (PROUVÉ, 2026-09-25)
(L'argument n'utilise pas les signes S : il combine la saturation du budget des zéros (18.2) et la 2-régularité de R_s.)
**Lemme 19.1 (orbites libres en un bloc).** En un bloc s, appelons *libres* les orbites passant par s qui ne sont pas appariées par s.
Sous R_s elles forment un graphe 2-régulier (cycles pairs, alternés m/π, de longueur ≥ 4 : les digones sont exactement les paires D
via s, Lemme 3), sans arête entre deux orbites sœurs d'une paire à n_s = 0, et deux orbites sœurs (n_s = 0) y ont au plus un voisin
commun (un tel voisin est le μ-partenaire de la classe partageant l'autre bloc, forme élémentaire de 11.5). En particulier les orbites
libres en s ne peuvent pas être exactement les 2 orbites d'une paire, ni exactement les 4 orbites de deux paires (K₂,₂ interdit).
**Corollaire 19.2 (configurations ancrées).** Dans une configuration ancrée (18 ; chaque paire a ses deux orbites libres en son bloc non
ancre), le bloc s porte 2(6 − d_s) orbites libres, donc d_s ∉ {4, 5} ; d_s ≠ 1 par définition. D'où d_s ∈ {0, 2, 3, 6} et la suite des
degrés entrants est 3⁷ (tournoi régulier), (6,3,3,3,2,2,2) ou (6,0,3,3,3,3,3).
**Théorème P.** Aucune configuration ancrée sur un tournoi régulier n'existe ; en particulier le cas de Paley est impossible.
*Preuve.* (i) *Saturation.* d ≡ 3 : chaque étoile est un triangle D. Pour un carré en t de bras j → k, l'orbite de {t,j} (troisième bloc
k) a |D_k| = 3 − 2 = 1 et l'orbite de {t,k} (troisième bloc j) a |D_j| = 3 (t ∉ In(k) car k → t). Les 21 orbites avec |D| = 3 ont z ≥ 1
(18.2) ; chaque carré a au plus une orbite zéro : injection de 21 orbites dans 21 carrés, donc **chaque carré a une orbite zéro, donc
ω(Q) = 0 pour les 21 carrés** (une orbite zéro exige ω(Q) = 0).
(ii) *Cibles interdites.* En un bloc s, Out(s) = {t₁,t₂,t₃}, paires libres Π_m = {s,t_m} (ancrées en t_m). L'orbite de Π_m appariée à
{t_m, c} (c ∈ In(t_m) ∖ {s}) a pour troisième bloc c ; si c ∈ Out(s), ses arêtes R_s vers Π_c compteraient dans ω : par (i) elles sont
interdites, et ses deux voisins R_s sont alors les deux orbites de la paire restante. Soit f(s,t) = |In(t) ∩ Out(s) ∖ {s}| (nombre de
chemins s → c → t) = nombre d'orbites de {s,t} ayant une cible interdite.
(iii) *f ≤ 1.* Si f(s,t₁) = 2 (t₂ → t₁, t₃ → t₁), l'orbite de Π₁ de troisième bloc t₂ a ses deux voisins dans Π₃ et l'autre les a
dans Π₂. Si t₂ → t₃, l'orbite de Π₃ de troisième bloc t₂ a sa cible Π₂ interdite, donc ses deux voisins dans Π₁ : impossible, l'une des
orbites de Π₁ est déjà saturée par Π₂. Si t₃ → t₂, symétrique. Donc f ≤ 1 partout.
(iv) *f ≡ 1 impossible.* Σ_{arcs} f = nombre de triples transitifs = Σ_v C(3,2) = 21 = nombre d'arcs, donc f ≡ 1. En un bloc s, chaque
Π_m a alors une orbite X_m forcée vers Π_{φ(m)} (φ(m) ≠ m). X_m est voisine de X_{φ(m)}, dont les deux voisins sont dans Π_{φ²(m)} ;
il faut donc φ²(m) = m pour tout m : φ serait une involution sans point fixe de {1,2,3}. Impossible. ∎
*Remarque.* Le seul tournoi doublement régulier sur 7 sommets est celui de Paley (f ≡ 1) ; les autres tournois réguliers tombent déjà
en (iii). Le Théorème P ferme le prototype ouvert au §18.
**Reste parmi les configurations ancrées :** (6,3,3,3,2,2,2) (étoile totale en un bloc, trois étoiles triangulaires, trois paires liées)
et (6,0,3,3,3,3,3) (un bloc ancre tout, un bloc n'ancre rien, les cinq autres forment le tournoi régulier sur 5 sommets). Le budget des
zéros n'y est pas saturé (au plus 9 orbites |D| = 3 dans le second cas), donc ω n'est pas forcé à 0 : il faut y raffiner (ii)–(iv).
Hors configurations ancrées : chaînes tordues de longueur ≥ 4, mélanges, carrés disjoints ; 19.1 s'y applique tel quel.

## 20. Arrêt de la piste « configurations », retour au cœur du graphe (2026-09-25)
**Suite des pistes (IMPASSE honnête).** Cas ancré (6,3,3,3,2,2,2) : une structure liée (u; s, w) avec s → w impose out(s) ≥ 4 (sinon
K₂,₂ en s, 19.1), d'où un triangle cyclique L₁ → L₂ → L₃ → L₁ et In(L_j) = {L_{j−1}, q_j}, q_j ∈ C, L_{j−1} → q_j ; de même un bras C d'une
composante liée de l'étoile en A doit être battu par son partenaire. Cas (6,0,3⁵) : 9 orbites |D| = 3 seulement, budget non saturé,
aucune contradiction. Le mécanisme du Théorème P (ω ≡ 0 forcé) ne s'y applique pas. Piste arrêtée ici, sur consigne de l'opérateur.
**Le cœur, sans involution.** Autour d'un sommet x : N(x) = 7K₂ (blocs) ; les 84 sommets à distance 2 sont en bijection avec les paires
de voisins de x non adjacents, c'est-à-dire avec les **arêtes du graphe cocktail CP(7) = K₇ₓ₂** sur les 14 voisins. Γ existe si et
seulement s'il existe un graphe H (84 sommets, 12-régulier) sur E(CP(7)) avec, N étant l'incidence arête–sommet de CP(7) et P l'échange
des partenaires de bloc :
  (E1) HN = 2J − N(I + P)   (u = {a,b} a 1 voisin contenant a, 1 contenant a′, 2 contenant tout autre élément) ;
  (E2) H² + H = 12I + 2J − NNᵀ   ((H² + H)_{uw} = 2 − |u ∩ w| pour u ≠ w : c'est tout λ = 1, μ = 2).
Spectre forcé : {12, (−2)⁶, 0⁷, 3⁴⁰, (−4)³⁰} (les −2 sont les vecteurs d'incidence de bloc, les 0 les vecteurs c_t, 7 dimensions).
**Plan.** Pour une paire de blocs P = ij, les 4 arêtes {a,b}, {a,b′}, {a′,b}, {a′,b′} et x, a, a′, b, b′ forment une grille 3×3 ou non ;
« x est plat » si les 21 paires sont des grilles et si chaque couplage entre deux grilles disjointes (parfait, par 9.1) préserve les
droites (isomorphisme de C₄).
**Théorème D′ (PROUVÉ, aucun automorphisme supposé).** Dans une srg(99,14,1,2), aucun sommet n'est plat.
*Preuve.* Si x est plat, les vecteurs χ_P(u) = (−1)^{ε_i(u)+ε_j(u)} (21, orthogonaux à 1, n_t, c_t) engendrent un espace H-invariant sur
lequel H = A_s − 2I, A_s signature de Kneser KG(7,2) (les couplages préservant les droites envoient χ_P sur ±χ_{P′}). Le spectre (E2) y
impose {3, −4}, soit A_s² = 3A_s + 10I, et l'étape (6) du Théorème D (équilibre des triangles, 2(t₁+t₂+t₃) ≠ 0) conclut. ∎
(Le Théorème D en est le cas « involution » : un couplage commutant avec σ préserve les antipodes, donc les droites.)
**Lecture.** Les deux membres connus de la famille (m = 2, 11) sont des graphes de Cayley sur F₃ⁿ issus de codes ternaires : tous leurs
sommets sont plats. Les cas ouverts sont exactement 7 | m (v = 99, 6273, 494019, aucun n'est une puissance de 3). Un graphe à 99 sommets
devrait donc être **non plat en chaque sommet** : un objet sans aucun modèle affine connu. Piste de cœur : montrer que λ = 1, μ = 2
force la platitude (au moins en un sommet), ou mesurer globalement le défaut de platitude (grilles comptées : 9 par grille).
**Crédit et correction (2026-09-25, après lecture en seule lecture du dépôt amazing-source/conway99-c1, commit 9afb4d2).** Le
Théorème D′ est un cas particulier du **L8 de C1** (« aucun point de Hamming », notes/105 de C1) : C1 exclut tout sommet dont les 21
paires sont des grilles, sans supposer que les couplages préservent les droites ; C1 a aussi L10 (au plus 18 grilles par sommet) et
L16 (Krein : P ≤ 1151, soit Σ|W| ≥ 1410). La question « demande de réparation de la parité contre offre » coïncide avec le verrou
déclaré de C1 : il manque un **plafond de torsion** (P ≥ 1152). Les « 1 522 solutions locales » de C1 sont des modèles en un point de
Hamming, désormais exclus par L8 : ils ne servent pas de cartes tordues. Seul élément possiblement transférable d'ici : dans le
Théorème P, le plafond venait d'une ressource bornée (au plus une orbite zéro par carré presque parfait). Analogue général : un
9-ensemble à e arêtes a ‖E₋₄1_S‖² = (36 − 2e)/7, donc une déviation quadratique totale ≤ 7(36 − 2e) chez les sommets extérieurs.
Réserve de C1 : un Γ tordu « générique » aurait P ≈ 126, donc aucune moyenne ne donnera le plafond ; il devrait venir d'une rigidité.

## 21. Tentative « le désordre sature le plafond 7-adique » (2026-09-25) — IMPASSE, déjà cartographiée par C1
Idée : k = rank₇(A − 3I) − 15 a un plafond (Brouwer–van Eijl : k ≤ 30 ; L22 de C1, par le genre de Λ₋₄ et la formule d'oddité :
**k ≤ 28**) et la torsion pousse k vers le haut ; il suffirait de montrer rank₇(M_x) ≥ 29 en un sommet (M_x = complément de Schur
84 × 84 de la boule de x, C1). Constat après lecture de C1 (notes/97, 99, 101) :
- l'idée elle-même est déjà chez C1 (L22 §5 : « la structure générique est interdite ») ;
- les minorations locales prouvées s'arrêtent à 22 (L21, une droite générique) ; un bloc entièrement déterminé par la structure
  locale plafonne à 24 (deux étoiles couplées, 24 × 24) : **aucune contradiction possible avec des blocs connus** ;
- toute preuve de k ≥ 29 doit donc passer par les arêtes génériques inconnues (technique des relations de C1 sur quatre étoiles :
  44 lignes, au moins 16 relations de lignes complètes forcées à k = 28), c'est-à-dire un gros calcul exact, ou une idée neuve.
Reformulation qui pourrait servir d'idée neuve (non exploitée) : à k = 28, ker₇(M_x) dépasse le noyau rationnel de **exactement
2 dimensions** en chaque sommet (« accidents 7-adiques »), et ces accidents portent la partie U de Λ₋₄, **anisotrope** (L22 : U non
hyperbolique ; en dimension 2, anisotrope). Montrer qu'un accident est toujours isotrope (z·z ≡ 0) exclurait k = 28 ; montrer qu'ils
contiennent un sous-espace totalement isotrope de dimension (30 − k)/2 exclurait tout k. Aucun argument en ce sens pour l'instant.

## Statut de relecture externe (mis à jour le 2026-09-25)
Il n'existe pas de contrôle positif en N = K₇ (BvLS n'a aucune torsion). La relecture indépendante est donc le seul garde-fou.
- **Théorème P (§18–19) : VÉRIFIÉ EN ENTIER**, maillon par maillon, par la session « structures » (sur b56e836). Maillons
  relus : 9.2(ii), 18.1, 18.2 (h_k = 2 + ω), |D_k| = d_k − 2[j → k], la saturation (i), les cibles interdites (ii), f ≤ 1 (iii),
  et f ≡ 1 impossible (iv).
- **Théorème F (§15–16) : VÉRIFIÉ EN ENTIER** par la même session : 16.1, étapes (a)–(d), lemme du coin, tableau des types,
  3-cycle du triangle complet (argument direct ajouté au §15), lemme de partition et comptages exacts 2 / 1 / 3. Signalement
  de notation : « (G²)_{A₁B₂} = 3 − 2e » est correct avec G = Q_out − 2D (arêtes D exclues) ; la relecture appliquait à G l'identité de Q_out
  ((G² + G)_{UV} = 4 − Sh_UV), d'où 3. Définition de G désormais explicite au §0.
- Borne g(x) ≤ 14 au point fixe (§11–13) : conclusion vérifiée ; deux trous de rédaction corrigés (b56e836). Preuves internes de
  GT, 11.4–11.7, 12.1–12.2 et du Théorème D : logique suivie par la relecture, sans re-dérivation.
- **Version sans σ du Lemme 16.1** (établie le 25/09 avec la session « structures »). Pour une paire de classes R en x (4 étiquettes) :
  Σ_{w extérieur ∉ R} C(n_R(w), 2) = 8 − s_R − a_R − i_R.
  - s_R : paires sœurs adjacentes.
  - a_R : paires antipodales adjacentes.
  - i_R : couples (paire d'étiquettes de R, voisin commun situé dans R).
  Contrôles : grille (4, 0, 4) → 0 ; srg(9,4,1,2) et les 55 paires de BvLS → 0 échec (« structures »). Le terme a_R n'a aucun contrôle
  positif : les deux graphes connus ont a = 0. En torsion maximale de C1 (s = a = i = 0), le total vaut 8.

## Audit adverse d'Astra (autre modèle, 25/09, relayé par « structures ») — bilan
- **Erreur réelle** : Théorème A (§2) et A′ (§4), corrigés ci-dessus (e′ impair, α = e′ − 1 pair). Le Théorème B n'est pas touché ;
  GT, F, P, D n'utilisent pas A.
- Notation corrigée : K : Alt → Sym′ (§8), sans changement des identités.
- **Tient** : 12.1, 12.2 ; g(x) ≤ 14 ; 9.2(ii) ; lemme du coin ; 11.4–11.7 ; F (collinéation = réétiquetage de l'argument) ;
  P ; D. 17.3–17.4 tiennent dans le cadre « même β » : ne pas retirer l'hypothèse sur les grappes.
- Bibliographie : le texte de Paduchikh identifié porte sur (85,14,3,2) et ne se transpose pas ; citer Cesarz–Woldar dans sa
  version publiée (Algebraic Combinatorics 8(2), 2025) ; ajouter Wilbrink 1984 (pas d'ordre 11) ; Thakkar–Severini (arXiv
  2608.11211 v2) décrit encore le cas Z₇ à un point fixe comme ouvert, à tort (Behbahani–Lam).

## 22. Réduire Γ à des problèmes et résultats déjà connus (consigne du 25/09 : pas de moteur, des leviers existants)
Principe : exhiber dans Γ une sous-structure qui serait un objet classique dont l'inexistence ou la classification est acquise.
**Déjà exploité dans le corpus.**
- Lam (pas de plan projectif d'ordre 10) : tissu de 11 grilles à couplages anti-isomorphes (L25 de C1).
- Lloyd, Tietäväinen, Makhnev : cas plat (quotient de H(7,3)).
- Unicité du STS(7) (Fano) : Théorème F, L8 de C1.
- Théorème de Hall (injection de demandes dans des carrés) : Théorème P.
- Classification des tournois réguliers à 7 sommets : Théorème P.
**Candidats examinés ici.**
1. **Ensembles intrigants** (théorie des tight sets et m-ovoïdes, géométrie finie). [RÉSULTAT] Une grille S est 3-intrigante (4 voisins
   dedans, 1 dehors, par 9.1) et un ensemble c-serré T est (−4)-intrigant (|T| = 11c). Leurs vecteurs centrés vivent dans des espaces
   propres différents, donc sont orthogonaux, d'où **|S ∩ T| = |S||T|/99 = c**. En particulier toute grille coupe toute coclique de
   Hoffman en exactement 2 sommets ; un tissu de 11 grilles la découpe en 11 × 2 = 22. Mord seulement si les deux objets existent.
2. **Coclique de Hoffman C (22 sommets).** [RÉSULTAT]
   - Les ensembles N(u) ∩ C (u ∉ C) forment un design 2-(22, 4, 2) à 77 blocs, tous distincts, deux à deux d'intersection ≤ 2
     (μ = 2 et C indépendant).
   - Les 77 sommets extérieurs portent exactement 77 droites entièrement extérieures (chaque sommet extérieur en a 3 ; les 4 autres
     droites passent par ses 4 voisins dans C) : une configuration (77₃).
   - Reste à trouver une classification connue de ces designs (question bibliographique pour Astra). Ici non plus, rien ne force
     l'existence de C.
3. **Plans projectifs, Bruck–Ryser–Chowla, designs symétriques** : aucune autre source de design symétrique que les tissus. Rien.
4. **Théorèmes min-max** (Hall, König, Tutte, flot maximal/coupe minimale). Les budgets du §16–19 sont des problèmes de b-couplage :
   une coupe minimale violée est une contradiction. [PISTE] Écrire les deux cas ancrés restants, (6,3,3,3,2,2,2) et (6,0,3⁵), comme
   problèmes de transport, et chercher une coupe violée (à la main : ce sont des réseaux de taille ≤ 42).
**Verdict provisoire.** Les réductions qui mordent exigent toutes l'existence d'une sous-structure (tissu, grille, coclique de
Hoffman) que rien ne force. On retombe sur le même verrou : forcer l'existence d'un objet serré. Seule la piste 4 évite ce verrou,
mais elle reste dans le cadre C2 des configurations.
**Cas ancré (6,3,3,3,2,2,2) par la voie du transport (§22, piste 4), 25/09.**
- [RÉSULTAT] **Le sous-cas « triangle des C transitif » est impossible.** Si C₁ → C₂ → C₃ et C₁ → C₃, les arcs des C vers les L
  sont exactement les q_j → L_j, donc C₃ serait q_j pour deux valeurs de j. Or le degré entrant de C₃ depuis les L vaut 3 − 2 = 1, alors
  que chaque q_j = C₃ exige L_{j−1} → C₃. Contradiction.
- [RÉSULTAT] **Le tournoi restant est unique** (indices mod 3) : C_j → C_{j+1}, C_j → L_j, L_j → L_{j+1}, L_j → C_{j+1},
  L_j → C_{j−1}, et tout → A. Les étoiles : liée (L_j; L_{j−1}, C_j), de queue L_{j−1} ; triangle en C_j sur {C_{j−1}, L_{j−1}, L_{j+1}}.
  Dans l'étoile en A, les composantes liées ont une queue L (out(C) = 3 interdit une queue C).
- [RÉSULTAT] **Demandes.** Les orbites |D| = 3 sont les X_j ∈ {C_j, L_{j−1}} (3 au total) et les orbites de l'étoile en A dont la queue
  est un C. Il y en a au plus 4 : avec e_CC = e_LL (égalité des degrés), t_C ≤ e_CC + min(6 − 2e_CC, 3) ≤ 4. Une double arête C_x–L_x
  aurait une queue C, ce qui est interdit. Total des demandes ≤ 7, pour une capacité de 9 carrés (3 étoiles en C).
- [RÉSULTAT] **Blocages.** Dans chaque carré d'étoile en C, une seule orbite peut porter ω : celle de {C_j, C_{j−1}} en C_{j−1} pour S^a,
  et celle de {C_j, L_{j+1}} en L_{j+1} pour S^b et S^c. Au bloc L_j, le graphe libre est forcé biparti entre
  {U, U′, V₁, V₂} et {A₁, A₂, Z^b, Z^c}, mais une configuration avec ω = 0 partout existe localement :
  U–Z^b–U′–A₂–V₂–Z^c–V₁–A₁–U.
  Au bloc C_j, ω(P^a) = 0 force l'unique 6-cycle P^a–A₁′–Y–P^b–Y′–A₂′ ; ce n'est pas une contradiction.
- **Verdict (IMPASSE honnête).** Aucun blocage n'est forcé, donc la capacité (9) couvre la demande (≤ 7) : Hall ne mord pas. Pour
  conclure, il faudrait forcer au moins 3 blocages (ω = 1, ou un apex hors de F_β), et rien ne les force.
**Cas ancré (6,0,3⁵) par la même voie** (travail de la session « structures », vérifié par échange le 25/09). [RÉSULTAT], négatif :
- L'étoile en A est un 6-cycle, 3 + 3, ou (A; Z, B_x) plus un 4-cycle : 2 + 2 + 2 est exclu, car une composante liée a une queue de
  degré sortant ≥ 4, et seul Z en a un.
- Il y a 9 demandes : une par étoile en B, dirigée vers le bloc t − 2, et 4 venues de l'étoile en A. Chacune n'est compatible qu'avec
  les 3 carrés d'un seul bloc B.
- Hall par bloc donne demande(b) ≤ 1 + 2 = 3 = offre(b) : aucune violation.
- Aux blocs saturés (« tendus »), la seule contrainte forcée est ω(P₂b) = 0. Elle est réalisable dans 19.1.
- Les égalités des orbites |D| = 0 sont automatiques. Plus généralement, l'identité 18.2 n'a de contenu que par ses plafonds : au
  plus une orbite zéro par carré ; un zéro exige ω(Q) = 0 et W ∈ F_β ; au plus une orbite apex par carré.
- Cascade apex → zéro (seul couplage non vide trouvé) : une étiquette apex w d'un carré Q via k, avec k(w) = k, a ω = o = 0 et
  est l'orbite zéro des |D_k(w)| − 1 autres carrés de D_k(w).
**Bilan de la piste 4 (IMPASSE).** Aucun des deux cas ancrés restants ne se ferme avec les ingrédients actuels : il reste du jeu
(9 demandes pour 15 carrés, et au plus 7 pour 9), et aucun blocage n'est forcé. Pour que Hall morde, il faudrait forcer des blocages
(ω = 1, ou un apex hors de F_β), c'est-à-dire encore forcer une structure : c'est le même verrou sous une autre forme.

## 23. La démarche de C3 transposée à C2 (consigne relayée par C1, 25/09) : cohérente, sans contradiction
Démarche (d'après la description relayée par C1 ; le dépôt C3 n'a pas été ouvert, faute d'accord direct de l'utilisateur) : on
compare deux spécialisations d'un même module hermitien libre sur l'anneau local Z_p[G], dont les discriminants sont congrus, puis on
exploite l'asymétrie globale (seul le caractère trivial voit j = 1, de norme 99 = 9·11).
[RÉSULTAT]
1. **Seul p = 2 compte.** Pour p impair, Z_p[C₂] ≅ Z_p × Z_p, les deux spécialisations sont indépendantes, et il n'y a aucune congruence.
2. **Invariant local à p = 2.** Pour θ ∈ {3, −4}, F = E_θZ₂⁹⁹ est libre sur Z₂[C₂] (quotient de Brauer nul : le numérateur de
   (E_θ)_xx est pair) et pair (car ⊥ 1). Posons R = det(Λ_θ^+)/det(Λ_θ^−) (parties unités 2-adiques). Alors
   **R ≡ 1 + 2t (mod 4)**, avec t = (E_θ)_xx/2 (mod 2).
   *Preuve.*
   - det(H₋ + 2S′) ≡ det H₋·(1 + 2 tr(H₋⁻¹S′)) (mod 4).
   - tr(H₋⁻¹S′) ≡ h̄(u, w̄) = h̄(w̄, w̄), où u est le vecteur de Wu de h̄ et w̄ représente q(x) = b(x, σx) (mod 2), qui est linéaire.
   - Sur Z⁹⁹, q(x) ≡ x_x, d'où (1 + σ)w̃ = E_θe_x et h̄(w̄, w̄) = (E_θ)_xx/2.
   *Tests* (code/z2/c3_transpose_test.py), 4 sur 4 conformes : Paley 9 (θ = 1 et −2 : t = 0, R = 1) ; BvLS (θ = 4 : t = 0, R = 1 ;
   θ = −5 : t = 1, R ≡ 3 mod 8).
3. **Global pour 99.**
   - det Λ₃^− = 2²⁷·7^r (réseau de S sur Z⁴⁹, r pair par le Théorème B).
   - det Λ₃ = 7^α·11 (recollement d'indice 9 avec Z1, α pair par le Théorème A corrigé), et [Λ : Λ⁺ ⊕ Λ⁻] = 2²⁷.
   - D'où R₃ = 7^{α−2r}·11 ≡ 3 (mod 8), alors que le local donne t = 3/11, impair, donc R ≡ 3 (mod 4) : **cohérent**.
   - Pour θ = −4, R = 9·7^{α−2r} ≡ 1, et t = 2/9 : **cohérent**.
   - Mod 8, le local ne tranche pas entre 3 et 7 pour t = 1 (blocs de rang 1 ; BvLS réalise 3). Les invariants plus fins sont ceux du
     genre 2-adique, déjà épuisés au §4.
**Verdict (IMPASSE).** L'asymétrie « j voit 11 » existe à p = 2, mais elle tombe exactement sur la valeur locale forcée.
**Complément après lecture de la preuve de C3 (dépôt conway99-c3-experimental, lu pour la démarche seulement, avec l'accord de
l'utilisateur).**
- La démarche de C3 : les deux spécialisations de E₃ (quotient Q et matrice d'Eisenstein M) donnent **le même espace quadratique
  sur F₃**, d'où l'égalité des symboles de Legendre de leurs discriminants. Or seul le côté trivial voit j (norme 33 = 3·11), et 11
  n'est pas un carré modulo 3.
- **Pourquoi l'analogue échoue en C2 :**
  1. à p = 2, F₂^× est trivial, donc l'étape « même espace mod p » ne dit rien ;
  2. au premier ordre (mod 4), le point fixe unique fournit exactement le terme correctif.
- [RÉSULTAT] **Forme générale.** Soit une srg de v impair munie d'une involution à exactement un point fixe, et θ une valeur propre de
  multiplicité f avec E_θ 2-entier. Comme u₊u₋ = det Λ_θ, on a R ≡ det(Λ_θ) modulo les carrés, et donc
  **det(Λ_θ)_impair ≡ 1 + f/v (mod 4)**. Contrôles : Paley 9 et BvLS (4 cas sur 4).
- Pour 99 (θ = 3) : det Λ₃ = 7^α·11 et 1 + 6/11 ≡ 3 (mod 4), ce qui force α pair. **C'est exactement le Théorème A corrigé** : la
  démarche C3 transposée redonne une parité déjà connue, sans contradiction.
- [RÉSULTAT] Au passage, l'analogue du « point 1 » de C3 (§4) : les vecteurs invariants de Λ₃ ont une coordonnée paire au point fixe,
  par la ligne x de N = 11(A + 4I) − 2J, qui donne v_x ≡ Σ_{y∼x} v_y = 2Σ_blocs ≡ 0 (mod 2).
  (Plus simplement : 3v_x = Σ_{y∼x} v_y = 2Σ, et θ = 3 est impair. Le fait est trivial.)
- **CORRECTION (25/09, laboratoire `Desktop/c2_experimental`, note 24 et N24).** La « forme générale » ci-dessus est une
  **tautologie** et n'exclut jamais rien.
  - Λ_θ ⊗ Z₂ est unimodulaire **pair** (E_θ1 = 0 est un vecteur caractéristique), donc det ≡ (−1)^{f/2} (mod 4) sans involution.
  - Dès que f est pair, 1 + f/v ≡ (−1)^{f/2} (mod 4).
  - Le cran suivant, l'invariant d'Arf, a une formule locale : Arf = w̃·w̃/2, contrôlée dans 50 cas sur 50 sur de vrais graphes.
    Pour 99, il est automatique pour Λ₃ ; pour Λ₋₄, il équivaut à σ_S ≡ 1 (mod 8), déjà établi au §4.
  - **Antériorité** : le laboratoire C2 du 23/09 (`c2_experimental`, N09 et notes 09–10) avait déjà conclu que le transfert de C3 est
    vide pour une involution et qu'aucune obstruction de genre ne tuera C2. Le §23 le redémontre sans l'avoir su.

**No-go de « structures » (reçu le 25/09, qui explique pourquoi la transposition de C3 ne mord pas).**
- La démarche C3 demande un opérateur entier M ≡ Q (mod 𝔭) vérifiant (M + 4)(M − 3) = 0, sans terme en J.
- Or (A + 4I)(A − 3I) = A² + A − 12I = 2J (vérifié : A² = 12I − A + 2J pour (14, 1, 2)). Sur les orbites de σ, c = 4. Réduire
  mod 𝔭 force donc 𝔭 | 2 : seul p = 2 est possible.
- En 2, on ne voit que la forme bilinéaire, donc au plus la parité, c'est-à-dire le Théorème A. L'invariant d'Arf demanderait une
  congruence mod 4, qui échoue car 2J ≢ 0 (mod 4).
- Il faudrait une symétrie d'ordre impair pour que le symbole de Legendre morde : c'est le cas C3.
- « structures » a aussi exclu le candidat sans symétrie H = (−A + i√3·T)/2, où T oriente cycliquement les 231 droites (triangles).
  **Vérifié ici, à la main :** H n'a pas de spectre entier à au plus 3 valeurs propres.
  - H est hermitienne à coefficients dans Z[ω] (ω sur un arc orienté, ω̄ en sens inverse) et H ≡ A (mod π). Donc χ_H ∈ Z[x] et
    χ_H ≡ χ_A ≡ x⁵⁴(x + 1)⁴⁵ (mod 3). Avec un spectre entier : 54 valeurs propres ≡ 0, 45 valeurs propres ≡ 2, aucune ≡ 1 (mod 3).
  - tr H = 0 ; tr H² = 99·14 = 1386 ; tr H³ = 6·231 = 1386 (holonomie ω³ = 1 sur chaque droite).
  - H1 = −7·1 (T1 = 0), et 7I + H = N_ω N_ω* ⪰ 0, donc −7 est la plus petite valeur propre.
  - *Deux valeurs propres {−7, a}.* tr[(H + 7)(H − a)] = 1386 − 693a = 0 donne a = 2, et −7 ≡ 2 ≡ 2 (mod 3) : il n'y a pas
    54 valeurs propres ≡ 0. Exclu.
  - *Trois valeurs propres {−7, a, b}.* tr[(H + 7)(H − a)(H − b)] = 1386(8 − a − b) + 693ab = 0, soit **(a − 2)(b − 2) = −12**.
    Avec a, b ≥ −6, cela laisse {1, 14}, {4, −4}, {0, 8}, {5, −2}, {−1, 6}.
    - La règle mod 3 (aucune valeur ≡ 1) élimine les trois premiers et {5, −2}, puisque 1, 4, −2 ≡ 1. Il reste {0, 8} et {−1, 6}.
    - {0, 8} : m₀ = 54, m₋₇ + m₈ = 45, tr H = 0 donne m₋₇ = 24 et m₈ = 21. Mais tr H² = 1176 + 1344 = 2520 ≠ 1386.
    - {−1, 6} : m₆ = 54, m₋₇ + m₋₁ = 45, tr H = 0 donne 6m₋₇ = 279. Pas entier.
  - Donc aucun spectre, pas même à trois valeurs propres, sans recherche.
- **Conclusion commune :** la famille genre/réseaux est épuisée pour C2 ; le verrou reste combinatoire.

## 24. Lemme de la queue : aucune configuration ancrée n'existe (25/09, croisement avec le labo `c2_experimental`)
**Origine.** Le **lemme d'équilibre** du labo (note 23 §8, preuve humaine, sans involution), appliqué aux configurations ancrées des
§18–§22. Rédaction détaillée : labo, `notes/25_LEMME_DE_LA_QUEUE.md`.

**Lemme 24.0 (équilibre, labo).** Soit ℓ une droite, et C une réunion de composantes du graphe 3-régulier induit sur N(ℓ). Tout
sommet à distance 2 de ℓ (60 sommets) a exactement |C|/6 voisins dans C.
*Preuve.* BBᵀ = (12I − A_NN + 2J) − Cls − A_NN², et C coupe chaque classe en tiers. Donc Σ_z n_C(z)² = 5|C|²/3 et
Σ_z n_C(z) = 10|C| : la variance est nulle. ∎ (C'est propre à k = 14 : 1_C − (|C|/36)1_N est un vecteur propre de A_NN de valeur
propre 3.)

**Composantes en un bloc k.**
- N(ℓ_k) contient les intérieurs b ∈ X_k = {±e_m : m ≠ k} et les racines ±e_k + b.
- Les trois sommets d'étiquette b sont dans la même composante, et b, −b aussi (la droite b – (−b) passe par x).
- Les composantes sont donc les orbites de ⟨σ, μ_k, τ_k⟩ sur X_k. Ce sont des réunions de blocs : c'est un graphe G_k sur les 6
  autres blocs.

**Lemme 24.1.** Soit X ⊂ ij une orbite appariée par i à D(X) ⊂ ik (troisième bloc k). Alors, dans G_k, la composante de i contient
au moins 3 + [j dans la même composante] blocs.
*Preuve.* x évite k. Dans cette composante, il voit forcément :
- l'intérieur de son signe en i ;
- les deux sommets v, σv de son carré (des racines d'étiquettes ±e_i, dans les classes ±e_k) ;
- l'intérieur de son signe en j, si j est dans la composante.

Or il y a exactement (nombre de blocs) voisins dans la composante, par 24.0. ∎

**Corollaire 24.2 (queue).** Dans une configuration ancrée (§18), X_In(k) est stable par μ_k = τ_k. En effet, une paire ancrée en k
est appariée par k à une paire elle aussi ancrée en k. Donc X_Out(k) est une réunion de composantes, de taille 2(6 − d_k). Pour un
carré de pivot i et de bras j, k : **d_k ≤ 3 − [k → j] et d_j ≤ 3 − [j → k]**. La queue de l'arc entre bras a donc d ≤ 2, et la tête
d ≤ 3.
Pour une structure liée, on retrouve la règle « s → w ⇒ out(s) ≥ 4 » du §20.

**Théorème 24.3.** Aucune configuration ancrée n'existe.
*Preuve.* Par 19.2, les degrés sont dans {0, 2, 3, 6}. En un bloc de degré 3, le triangle D porte des carrés sur les trois paires
de bras. Tout arc entre éléments de In(i) a donc une queue de degré ≤ 2, et au moins deux éléments de In(i) sont de degré ≤ 2.
- **3⁷** : impossible. C'est une nouvelle preuve du Théorème P.
- **(6,0,3⁵)** : In(B) = {Z, B₁, B₂}, et l'arc B₁–B₂ a une queue de degré 3. Impossible.
- **(6,3,3,3,2,2,2)** :
  - un In(C) contient au plus un C, donc les C forment un cycle ;
  - les deux L de In(C_j) pointent vers C_{j−1}, donc les trois C ont le même couple {L, L′} ;
  - les trois C pointent vers le troisième L, qui est alors de degré ≥ 3 : impossible.
  (Pour contrôle, dans le tournoi unique du §22, l'arc C_{j−1} → L_{j−1} dans In(C_j) a une queue C.) ∎

**Statut.** PROUVÉ et **VÉRIFIÉ deux fois** :
- par « structures » : relecture adverse étape par étape, avec contrôle numérique de l'identité d'équilibre générale
  BBᵀ1_C = (k − 14)1_C + (5/3)|C|1_N sur BvLS (891 droites, 8910 composantes) et contrôle négatif ;
- par une relecture interne indépendante.
Les deux cas ancrés qui restaient ouverts au §22 sont fermés.

*Simplification (relecture interne).* 19.2 est inutile :
- un degré 1 est impossible (paire sœur) ;
- tout bloc ayant un arc sortant est un bras, donc de degré ≤ 3 par 24.2 ;
- d'où d ∈ {0, 2, 3, 6} sans l'argument K₂,₂.

**Reste**, en Σ = ∅ :
- les paires tordues (§14) et les chaînes tordues de longueur ≥ 4 ;
- les mélanges ;
- les carrés disjoints, y compris dans les paires mixtes.

Les fourches (§16) sont ancrées, donc fermées.

24.1 s'y applique tel quel. Mais les composantes de G_k y sont plus grandes : les orbites libres forment des cycles de R_k (19.1),
et σ recolle les paires mixtes. Il faut donc contrôler G_k.
