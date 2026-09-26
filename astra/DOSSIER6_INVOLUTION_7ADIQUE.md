# Pour Astra : dossier n° 6 : tuer l'involution par l'arithmétique 7-adique (nouveau sujet, autonome)

**Contexte.** Le graphe de Conway Γ = srg(99,14,1,2) a Aut(Γ) ∈ {1, Z₂, Z₃} (littérature). Notre branche C3 a exclu
l'ordre 3, avec une preuve vérifiée en Lean. Si on exclut aussi une involution, **Γ n'a aucune symétrie**. Ce dossier
propose d'exclure l'involution d'un seul coup, par le genre 7-adique des réseaux propres, sans passer par l'énumération
des cas.

## 1. Cadre (tout est standard ; les preuves sont chez nous et disponibles)
- σ est une involution de Γ. Elle fixe exactement un sommet x (Makhnev–Minakova) et retourne chacune des 7 arêtes de
  N(x), appelées **blocs** {(c,+),(c,−)}.
- Les 84 sommets à distance 2 de x sont les **racines** : une par paire d'intérieurs de blocs différents, notée
  (a,s)(b,t). σ change les deux signes.
- Il y a 42 orbites de racines de taille 2 (deux par paire de blocs), et aucune racine fixe.
- **Partie impaire.** L⁻ = {f ∈ Z⁹⁹ : σf = −f}, de base d_o = e_o − e_σo pour les 49 σ-paires (7 blocs et 42 orbites).
  Avec la forme standard divisée par 2, L⁻ = Z⁴⁹ standard.
- A agit sur L⁻ par une matrice **S** (49 × 49, entière, symétrique), avec S_ab = A(rep a, rep b) − A(rep a, σ rep b).
  - **S² + S = 12I exactement** (J est nul sur L⁻). Valeurs propres : 3 (multiplicité 27) et −4 (multiplicité 22).
  - **Entrées.** On prend les représentants de C2 : (a,+)(b,ε) avec a < b, et (c,+) pour les blocs.
    - S_cc = −1 et S_cc′ = 0 ;
    - S_{c,o} = ±1 (le signe de l'intérieur de o en c) si c ∈ supp o, et 0 sinon ;
    - racines : S_oo = 0, S_{oo′} ∈ {0, ±1}, non nul si et seulement si o′ est l'une des 10 orbites à « un voisin »
      de o. Chaque orbite a en plus un unique partenaire D (deux voisins, S = 0).

## 2. Ce qui est démontré (« structures », 26/09 ; contrôlé sur deux vrais graphes)
- **(a) Rang modulo 7.** g⁻ := rang₇(S − 3I). On a (S − 3I)² = −7(S − 3I).
  - Λ₃⁻ = Z⁴⁹ ∩ ker(S − 3) (rang 27) et Λ₋₄⁻ (rang 22) se recollent le long de (Z/7)^{g⁻}, avec
    det Λ₃⁻ = det Λ₋₄⁻ = 7^{g⁻}.
- **(b) Genre.** Formule d'oddité appliquée à Λ₋₄⁻ :
  22 + 6g⁻ + 4[ε⁻ = −1] ≡ |E₋₄1|² = (147 − 1ᵀS1)/7 (mod 8), avec 1ᵀS1 = 497 − 4M. Ici M est le nombre de paires {o,o′}
  à S = −1 (nombre pair, par S² + S = 12I). On en tire :
  - **g⁻ est pair**, et [ε⁻ = −1] ≡ M − g⁻/2 ;
  - **g⁻ ≠ 22** : Λ₋₄⁻ serait 7·(unimodulaire), de signe +1 forcé ;
  - avec g⁻ ≥ 7 (les lignes des blocs sont indépendantes), **g⁻ ∈ {8, 10, …, 20}**.
  - Contrôle : la même formule est vérifiée exactement sur H(2,3) et sur le graphe de Berlekamp–van Lint–Seidel (243
    sommets), munis de x ↦ −x.
- **(c) Partie paire.** g⁺ (la même quantité sur la partie σ-invariante, sans la direction 1) est pair, avec
  2 ≤ g⁺ ≤ 22 et g⁺ + g⁻ = 14 + k, où k = rang₇(A − 3I) − 15 est pair et ≤ 28 (résultats de C1).
- **(d) Vecteurs propres universels.** Pour chaque bloc c, avec les représentants des orbites contenant c pris avec
  (c,+) :
  - z_c = 3d_c + Σ_{o∋c} d_o vérifie Sz_c = 3z_c ;
  - w_c = −4d_c + Σ_{o∋c} d_o vérifie Sw_c = −4w_c.
  Hors de la boule, cela découle de μ = 2 avec (c,±). Les z_c sont orthogonaux entre eux, et les w_c aussi.
- **(e) Minorant local.** Sur les 19 coordonnées d'une « droite-bloc » (d_c, les 6 autres blocs et les 12 orbites
  contenant c), S est déterminé localement. On a S_{o_v,o_w} = [τ(v) = w] − [ι(v) = w], avec τ et ι deux involutions
  des 12 étiquettes. Le rang₇ de S − 3I sur ces coordonnées vaut **toujours au plus 18** (z_c est dans le noyau), et 18
  génériquement.
- **(f) Signe à g⁻ = 18.** Si une droite-bloc a un rang 18, alors la forme B = S − 3I mod 7 est isométrique à sa
  restriction. Donc ε⁻ = χ(disc de B restreinte), et la genèse impose −1. Localement, ce signe vaut ±1 à peu près à
  moitié.
- **(g) Forme « réseau ».** N := le réseau de Gram 3I − S (vecteurs u_a, a ∈ 49 coordonnées) est entier, défini positif,
  de rang 22 et de **déterminant 7^{22 − g⁻}**.
  - Il est engendré par les 42 vecteurs u_o de norme 3, de produits scalaires −S_{oo′} ∈ {0, ±1}. Chacun a un produit
    ±1 avec exactement 10 autres.
  - Les relations 3u_c + Σ_{o∋c} u_o = 0, avec |u_c|² = 4, viennent de z_c.
  - **À g⁻ = 20, det N = 49. À g⁻ = 18, det N = 7⁴.**

## 3. Tâche
**Question principale.** Exclure g⁻ = 20 et g⁻ = 18, ce qui exclut toute involution. Ou prouver g⁻ ≥ 21.

Pistes, à toi de choisir :
1. **Réseau N.** Existe-t-il un réseau entier défini positif de rang 22 et de déterminant 49 (resp. 7⁴), engendré par
   42 vecteurs de norme 3 ayant ces contraintes ? Contraintes : produits 0 ou ±1, dix produits non nuls par vecteur, les
   7 relations 3u_c + Σ_{o∋c} u_o = 0 avec |u_c| = 2, et le motif de supports ci-dessus. Pense aux systèmes de vecteurs
   de norme 3, aux sous-réseaux engendrés, aux bornes de densité et de Hermite, aux congruences sur les nombres de
   vecteurs.
2. **Mineurs.** g⁻ ≤ 20 équivaut à la nullité mod 7 de tous les mineurs 21 × 21 de S − 3I. Les entrées « lointaines »
   (orbites de supports disjoints) sont les seules non déterminées localement. Chaque orbite a exactement 2 voisins
   contenant (t,+) et 2 contenant (t,−) pour tout bloc t hors de son support.
3. **Combinaison avec la partie paire** (c) ou avec les 7 droites-blocs à la fois (e), (f).

Comme d'habitude :
- les maths d'abord, et le calcul exact seulement pour vérifier ;
- distinguer « démontré » et « vérifié » ;
- faire passer les contrôles positifs m = 2 et m = 11. Ils ont une involution x ↦ −x, avec des valeurs propres
  différentes : 1/−2 et 4/−5, et un recollement 3-adique.
