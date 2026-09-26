# Trois problèmes ouverts précis (cas « involution » du graphe de Conway 99)

Pour un modèle extérieur (Sol, Astra) : à traiter en recherche indépendante. Chaque énoncé se suffit à lui-même.
- Distinguer « démontré » et « vérifié par calcul ».
- Donner tout calcul sous une forme revérifiable.
- **Ne plus travailler sur les vecteurs courts (normes 8 à 10) ni sur g⁻ = 12 : ces pistes sont closes** (voir le point 0).

## 0. Cadre et acquis
**Graphe et involution.**
- Γ = srg(99, 14, 1, 2), et σ est une involution de Γ. σ fixe exactement un sommet x (Makhnev–Minakova).
- **Blocs** : N(x) est formé de 7 arêtes {(c,+), (c,−)}, c = 0…6, et σ échange (c,+) et (c,−).
- **Racines** : les 84 sommets à distance 2 de x. Chaque racine a exactement deux voisins intérieurs, dans deux blocs
  différents. On écrit r = (i,α)(j,β), et r ↦ {(i,α), (j,β)} est une bijection vers les paires d'intérieurs de blocs
  différents. σr = (i,−α)(j,−β).
- 42 **orbites** {r, σr}, 21 **paires** de blocs, 2 orbites par paire (dites **sœurs**).

**Partie impaire.**
- L⁻ = {f ∈ Z⁹⁹ : σf = −f}, de base d_o = e_o − e_{σo} pour les 7 blocs et les 42 orbites. Avec la forme divisée par 2, c'est
  Z⁴⁹ standard.
- A y agit par une matrice entière symétrique S (49 × 49), avec **S² + S = 12I** et spectre 3²⁷, (−4)²².
- **Entrées**, en prenant pour représentants (c,+) pour les blocs et (a,+)(b,ε) pour les orbites :
  - S_cc = −1 et S_cc′ = 0 ;
  - S_{c,o} = ±1 (signe de l'intérieur de o en c) si c ∈ supp o, et 0 sinon ;
  - S_oo = 0, et S_{oo′} ∈ {0, ±1} ;
  - S_{oo′} = A(o, o′) − A(o, σo′) sur les représentants.
- **Rang** : g⁻ := rang₇(S − 3I).

**Démontré (et recoupé par deux codes indépendants) : g⁻ ∈ {16, 18, 20}.**
- Parité et genre 7-adique : g⁻ est pair et ≤ 20.
- Droite-bloc : pour chaque bloc c, la sous-matrice principale de S − 3I sur les 19 coordonnées {d_c, les 6 autres blocs,
  les 12 orbites contenant c} a un rang₇ égal à 18 − t_c, avec **t_c ≤ 4** (énumération exhaustive).
- Les cas t_c = 4 forcent une incompatibilité entre deux blocs voisins, d'où g⁻ ≠ 14.
- **Bloc « plein »** (les 6 paires en c appariées « via c », c'est-à-dire que les deux voisins par c de chaque racine de
  l'étoile sont dans une même orbite) : sa droite-bloc est de rang 18 et de signe local +1, **donc g⁻ = 20**.

## Problème 1 — exclure g⁻ = 16
**Formule.** Pour tout bloc c : g⁻ = 18 + t_c − 2r₀_c + ρ_c, avec ρ_c ≥ 0.
- t_c = dim ker Q_c, où Q_c est la matrice 12 × 12 (mod 7) de la boule impaire de c.
- r₀_c est la dimension des φ ∈ ker Q_c **équilibrés**, c'est-à-dire orthogonaux mod 7 aux 30 « ombres réduites » des orbites
  qui évitent c.

**Contrainte à g⁻ = 16.** Dans CHACUN des 7 blocs : t_c ≥ 2 et r₀_c ≥ 1 + (t_c + ρ_c)/2. Il faut donc de la torsion
équilibrée partout.

**Question.** Montrer que c'est impossible, par exemple par une incompatibilité entre les torsions équilibrées de deux blocs,
ou par un argument de réseau. Le réseau L_T (rang 15, det 2¹²·7^{g⁻−7}, min ≥ 8) donne seulement min L_T ≤ 16,9 à g⁻ = 16
par Blichfeldt.

## Problème 2 — le réseau de déterminant 49 (g⁻ = 20)
**Énoncé.** À g⁻ = 20, le réseau N de Gram 3I − S (49 générateurs) est :
- entier, défini positif, de rang 22, de **déterminant 49** ;
- de groupe discriminant (Z/7)² anisotrope ;
- de minimum ≥ 2, avec min N* ≥ 8/7.

Il est engendré par 42 vecteurs u_o de norme 3, et contient 7 vecteurs u_c de norme 4, avec :
- (u_a, u_b) = −S_ab ∈ {0, ±1} pour a ≠ b ;
- **Σ_a u_a u_aᵀ = 7 I₂₂** (repère serré, car (3I − S)² = 7(3I − S)) ;
- 3u_c + Σ_{o∋c} (±u_o) = 0 ;
- chaque u_o a un produit ±1 avec exactement 10 autres u_o.

**Question.** Montrer qu'un tel réseau muni d'un tel repère n'existe pas. Il suffirait déjà de l'exclure quand un bloc est
« plein » : alors les 12 u_o de ce bloc ont des produits explicites.

Pistes :
- N ⊕ diag(7,7) se recolle en un réseau unimodulaire impair de rang 24 (classification de Borcherds) ;
- séries thêta et formes harmoniques ;
- vecteurs de norme 2 : seulement deux profils possibles (14 × ±1, ou 2 sur un bloc + 10 × ±1).

## Problème 3 — l'étoile totalement libre (le cœur commun à tous les cas restants)
**Formulation compacte** (dérivée de λ = 1, μ = 2 autour de la droite {x, (c,+), (c,−)} ; validée sur le graphe de
Berlekamp–van Lint–Seidel).
- X est l'ensemble des 12 intérieurs des 6 autres blocs, écrits ±e_m, et σ(b) = −b.
- μ est un couplage parfait de X : μ(b) = b′ ⟺ (c,+)b ~ (c,+)b′.
- τ est un couplage parfait de X : τ(b) = d ⟺ (c,+)b ~ (c,−)(−d). τ n'a pas de point fixe.
- Φ associe à chaque racine lointaine z (les 60 racines qui évitent le bloc c, soit les paires {a, a′} ⊂ X de blocs
  différents) la paire Φ(z) = {b ∈ X : (c,+)b ~ z}.

**Contraintes.**
- (i) Φ est une bijection des 60 racines lointaines sur les 60 paires de X qui ne sont pas des paires de μ.
- (ii) Pour a ≠ b : #{z : a ∈ z, b ∈ Φ(z)} = 2 − [b = σa] − [b = μa] − [b = τσa]. Pour a = b : 1 − [τa = σa].
- (iii) #{z : b ∈ Φ(z), d ∈ Φ(σz)} = 2 − [d = σb] − [d = τμb] − [d = μτb] − [d = τb].

**Étoile totalement libre** : μ et τ n'ont **aucune arête commune**. Aucune orbite n'est appariée via c. Elle apparaît dans
tous les résidus actuels du cas involution.

**Question.** Montrer qu'aucun (μ, τ, Φ) totalement libre ne satisfait (i)–(iii). Ou, à défaut, en exhiber un : ce serait
une information précieuse.

**Faits utiles (démontrés).**
- Chaque étiquette b possède exactement une orbite lointaine {z_b, σz_b} avec b ∈ Φ(z_b) ∩ Φ(σz_b) : c'est son partenaire D.
  Les 12 orbites ainsi obtenues sont distinctes.
- ρ := Φ ∘ σ ∘ Φ⁻¹ est une involution sans point fixe des 60 arêtes de K₁₂ − μ, avec exactement 12 paires sécantes (une par
  étiquette) et 18 paires disjointes.

**Déjà tenté, sans succès.**
- Les solveurs SAT et CP génériques sur (i)–(iii) : délais dépassés. Même avec une solution plantée, ils en ratent un tiers.
- (iii) seul : souvent faisable.
- Conditions de moments d'ordre 2 et positivité des Gram : automatiquement satisfaites.

**Renfort possible, non encore exploité.** Les résidus réels. Le Gram 3I − S est de rang 22, et la droite-bloc de rang 18
sur ℝ. Les 30 orbites lointaines ont donc des résidus dans ℝ⁴ qui forment un repère serré Σ r rᵀ = 7I₄. Même chose côté pair,
avec 𝒫 = 6I − 2Q − Sh + J de rang 15 et l'étoile de rang 11 : Σ r rᵀ = 14I₄.
