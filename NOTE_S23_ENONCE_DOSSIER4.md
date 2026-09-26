# NOTE S23 : énoncé autonome du dossier n° 4 (carrés disjoints de C2)

Définitions de C2, puis la borne 3h ≤ 2n_dis et la réduction à cinq supports (Astra, dossier 4, §§ 6–7).

Rédigée pour C1, qui ne connaît pas le cadre de C2. Tout se déduit des axiomes de srg(99,14,1,2) et de l'involution, et
chaque fait porte sa preuve courte. Rien ici n'utilise le témoin (2,1) ni le lemme 11.2.

## 0. Cadre de C2
- **Graphe et involution.** Γ = srg(99,14,1,2), et σ est une involution de Γ avec Fix(σ) = {x}.
  - N(x) est formé de 7 arêtes disjointes (λ = 1) : ce sont les **blocs** {(c,+),(c,−)}, c = 1..7.
  - σ échange (c,+) et (c,−) dans chaque bloc (lemme de Makhnev, cadre de C2).
- **Racines** : les 84 sommets à distance 2 de x.
  - Une racine a μ = 2 voisins dans N(x), situés dans deux blocs différents : deux voisins dans le même bloc donneraient
    deux voisins communs, x et r, à une arête.
  - On écrit r = (i,α)(j,β) et **supp r = {i,j}**. L'application r ↦ {(i,α),(j,β)} est une bijection sur les 84 paires
    d'intérieurs de blocs différents : deux tels intérieurs ont exactement 2 voisins communs, x et une racine.
  - σr = (i,−α)(j,−β). Une racine a 2 voisins intérieurs et 12 voisins racines.
- **F0 : u ≁ σu.** Sinon, le triangle unique sur l'arête uσu serait σ-invariant. Son troisième sommet serait alors fixe,
  donc égal à x, or x ≁ u.
- **Orbites.** Ce sont les 42 paires {r, σr}. Le support est constant sur une orbite. Chaque paire de blocs porte
  exactement 2 orbites, dites **sœurs**.
- **Q_UW** := |N(u) ∩ W| pour u ∈ U. La valeur ne dépend pas du représentant et la matrice est symétrique. Q_UU = 0 par F0,
  et chaque ligne a pour somme 12.
- **Partenaire D.** N(u) ∩ N(σu) est σ-invariant, de taille μ = 2, et ne contient pas d'intérieur (u et σu n'en ont aucun
  en commun). C'est donc une orbite {w, σw} =: D(U).
  - Q_{U,D(U)} = 2, et c'est le seul 2 de la ligne : Q_UZ = 2 place z et σz dans N(u) ∩ N(σu).
  - D est un couplage parfait des 42 orbites, soit 21 **carrés** (u, w, σu, σw), qui sont des 4-cycles induits.
  - Type d'un carré : **sib** (supports égaux), **one** (un bloc commun), **dis** (supports disjoints).
  - **n_dis** est le nombre de carrés dis. Il y a donc 2n_dis orbites de type dis.
- **Sh(U,W)** := |supp U ∩ supp W|. En particulier Sh(U,U) = 2.
- **m et π.** Pour u = (i,α)(j,β) :
  - m_i(u) est l'unique voisin de u contenant (i,α) (λ sur l'arête u–(i,α)) ;
  - **π_i(u) est l'unique racine voisine de u contenant (i,−α)** (μ entre u et (i,−α) : les voisins communs sont (i,α) et
    π_i(u)) ;
  - de même en j, et π_i(σu) = σπ_i(u).
- **F1 (sommes par bloc).** Σ_{Z : b ∈ supp Z} Q_UZ vaut :
  - 2 si b ∈ supp U : ce sont m_b et π_b ;
  - 4 sinon : μ = 2 entre u et chacun de (b,±), et ni x ni un intérieur ne sont des voisins communs.
- **Identité d'orbites : Q² + Q = 12I + 4J − Sh.**
  - *Preuve.* Pour u ∈ U, on compte les chemins u–z–w′ avec z racine et w′ ∈ W : il y en a Σ_Z Q_UZ Q_ZW. C'est
    Σ_{w′} |N(u) ∩ N(w′)| moins les chemins par un intérieur. Ces derniers sont au nombre de Sh(U,W) : en un bloc commun,
    l'intérieur de u est dans exactement un des deux sommets w, σw.
  - Pour U ≠ W, on obtient 4 − Q_UW − Sh(U,W).
  - Pour U = W, on obtient 12 + 2 = 14 = 12 + 4 − 0 − 2. ∎
- **Budget de Gram** : c'est la forme hors diagonale de l'identité. **Pour U ≠ W, Σ_Z Q_UZ Q_ZW = 4 − Q_UW − Sh(U,W).**
  - Tous les termes sont ≥ 0 : si des termes connus atteignent déjà le budget, tous les autres sont nuls.

## 1. Apex d'un carré dis (C2, lemme 39.1)
Soit (U, W = D(U)) un carré dis, avec u = (i,α)(j,β), w = (k,γ)(l,δ) et T = {i,j,k,l}.
- **Apex.** c₁ est le troisième sommet du triangle sur l'arête uw, et c₂ celui du triangle sur wσu. Leurs images σc₁ et
  σc₂ sont les apex des arêtes σuσw et σwu.
- **Orbites d'apex** : A = {c₁, σc₁} et B = {c₂, σc₂}.
- **(a)** Un apex est une racine. Il ne contient aucun intérieur des deux sommets du carré qu'il voit. Par exemple, si
  c₁ ∋ (i,α), l'arête uc₁ serait dans deux triangles, {u,w,c₁} et {u,(i,α),c₁}.
- **(b)** Q_UA = Q_WA = Q_UB = Q_WB = 1 : c'est ≥ 1 par les triangles, et ≤ 1 car le seul 2 de la ligne est D.
- **(c)** Q_AB = 0.
  - c₁ ~ c₂ mettrait l'arête wc₁ dans deux triangles.
  - c₁ ~ σc₂ ferait de même avec l'arête uc₁.
- **Apex intérieur** : une orbite d'apex dont le support est contenu dans T.
  - Par (a), ce support est formé d'un bloc de supp U et d'un bloc de supp W. Sinon, l'apex serait σu ou σw, qui ne sont
    pas adjacents à u, ni à w.
  - Quitte à renommer les blocs et à changer de représentants, on prend **A = {i,k} et a = c₁ = (i,−α)(k,−γ)**.
  - On a alors **u = π_i(a) et w = π_k(a)** : u contient (i,α), l'opposé de (i,−α) ∈ a, et de même pour w.

## 2. Lemme : A est de type dis, et Q_UV = Q_WV = 0
Soit A = {i,k} un apex intérieur du carré dis (U,W), et V = D(A).
- **(i) supp V ∩ {i,k} = ∅**, c'est-à-dire que A est de type dis.
  - Supposons i ∈ supp V. Les sommets v et σv sont tous deux voisins de a, et l'un d'eux contient (i,α). Il vaut donc
    π_i(a) = u, et V = U.
  - C'est impossible : Q_AU = 1, alors que Q_AV = 2.
  - Même raisonnement pour k, avec w.
- **(ii) Q_UV = Q_WV = 0.**
  - Le budget de (U,A) vaut 4 − Q_UA − Sh(U,A) = 4 − 1 − 1 = 2. Il est déjà atteint par Z = W, puisque Q_UW·Q_WA = 2·1.
    Donc Q_UV·Q_VA = 2Q_UV = 0.
  - De même, le budget de (W,A) vaut 2 et il est atteint par Z = U.
- **(iii) A est apex intérieur d'un seul carré.**
  - Si A est apex intérieur d'un carré (X,Y), les sommets de ce carré voisins de a′ ∈ A sont π_i(a′) et π_k(a′) (même
    argument qu'au §1).
  - Leurs orbites sont U et W. Donc {X,Y} = {U,W}.

## 3. Théorème (Astra §7) : 3h ≤ 2n_dis
Ici, h est le nombre de carrés dis qui ont un apex intérieur.
- **U n'est apex intérieur d'aucun carré.**
  - Sinon, les deux sommets de ce carré voisins de u seraient π_i(u) et π_j(u), un dans chacune de ses deux orbites.
  - Or π_i(u) = a, car a ~ u, a ∋ (i,−α), et π_i(u) est unique. Ce carré est donc (A, V).
  - On aurait alors π_j(u) ∈ V, ce qui contredit Q_UV = 0.
  - Même raisonnement pour W : π_k(w) = a, et π_l(w) ∈ V contredirait Q_WV = 0.
- **Comptage.** Pour chacun des h carrés, on choisit un apex intérieur.
  - Les 2h orbites de ces carrés sont distinctes, car D est un couplage.
  - Les h apex choisis sont distincts, par (iii).
  - Aucune orbite n'est à la fois orbite d'un carré compté et apex intérieur, par le point précédent (un apex n'est
    jamais une orbite de son propre carré).
  - Les 3h orbites sont toutes de type dis : les carrés sont dis, et les apex le sont par (i).
  - Donc **3h ≤ 2n_dis**. ∎
- **Corollaire.** Au moins ⌈n_dis/3⌉ carrés dis n'ont pas d'apex intérieur.

## 4. L'outil : 𝒫 = 6I − 2Q − Sh + J est semi-défini positif
- **Identités de comptage :**
  - QJ = 12J ;
  - ShJ = 24J, car chaque bloc est dans 12 supports ;
  - Sh² = 10Sh + 8J. En effet, (Sh²)_UW = Σ_{b∈U, b′∈W} #{Z ⊇ {b,b′}} = 12|U∩W| + 2(4 − |U∩W|) ;
  - Q·Sh = Sh·Q = 8J − 2Sh, par F1.
- **En développant**, tous les termes commutent : 𝒫² = 84I − 28Q − 14Sh + 14J = 14𝒫.
  - 𝒫 est symétrique, de valeurs propres 0 et 14.
  - **Toute sous-matrice principale de 𝒫 est donc semi-définie positive.**
- **Entrées** : 𝒫_XX = 5 et 𝒫_XY = 1 − 2Q_XY − Sh(X,Y).

## 5. Réduction à cinq supports (Astra §§ 6.1–6.3)
**Notations.**
- U = {i,j}, W = {k,l}, A = {i,k} est l'apex intérieur, B l'autre orbite d'apex, et V = D(A).
- K = 𝒫[U,W,A] = [[5,−3,−2],[−3,5,−2],[−2,−2,5]]. Ses mineurs principaux dominants valent 5, 16 et 16, et
  K⁻¹ = (1/16)·[[21,19,16],[19,21,16],[16,16,16]].

**(6.1) supp B ∩ T = ∅.**
- Sur (U,W,A), la colonne de B est c = (−1−s, −1−t, 1−e), avec s = Sh(U,B), t = Sh(W,B) et e = Sh(A,B).
- Le résidu 5 − cᵀK⁻¹c vaut :

| support de B | résidu |
|---|---|
| B ∩ T = ∅ | 3 |
| {j,m} ou {l,m}, avec m ∉ T | −21/16 |
| {i,m} ou {k,m}, avec m ∉ T | −101/16 |
| {j,l} | −8 |
| {i,l} ou {j,k} | −15 |
| {i,k} (sœur de A) | −24 |

- Les cas {i,j} et {k,l} sont déjà exclus par le §1(a).
- Il reste **supp B = {p,q} ⊂ complément de T**. On note **r** le septième bloc.

**(6.3) Les supports possibles de V.**
- Par le §2(i), supp V ⊂ {j,l,p,q,r}, soit 10 candidats. De plus, Q_BV ∈ {0,1}, puisque D(V) = A ≠ B.
- Sur (U,W,A), la colonne de V est (a, b, −3), avec a = 1 − [j ∈ V] et b = 1 − [l ∈ V]. On a aussi
  𝒫_BV = 1 − 2Q_BV − |{p,q} ∩ supp V|.
- Le complément de Schur sur (B,V) est S = [[3, x],[x, y]], avec :
  - x = 𝒫_BV − 3 + (3/2)(a + b) ;
  - y = 5 − (21a² + 38ab + 21b² − 96(a + b) + 144)/16.
- Comme K ≻ 0, **𝒫[U,W,A,B,V] ⪰ 0 si et seulement si y ≥ 0 et 3y ≥ x².**

| supp V | y | x (Q_BV = 0 ; 1) | 3y − x² (Q_BV = 0 ; 1) | verdict |
|---|---|---|---|---|
| {j,l} | −4 | — | — | exclu |
| {j,p}, {j,q}, {l,p}, {l,q} | 11/16 | −3/2 ; −7/2 | −3/16 ; −163/16 | exclus |
| {j,r}, {l,r} | 11/16 | −1/2 ; −5/2 | 29/16 ; −67/16 | Q_BV = 0 forcé |
| {p,q} | 3 | −1 ; −3 | 8 ; 0 | reste |
| {p,r}, {q,r} | 3 | 0 ; −2 | 9 ; 5 | reste |

**Conclusion.** supp V ∈ {{p,q}, {p,r}, {q,r}, {j,r}, {l,r}}, avec Q_BV = 0 pour les deux derniers.

**Certificats entiers** (ordre U, W, A, B, V ; on donne x puis xᵀ𝒫x < 0) :
- {j,l} : x = (3,3,3,0,1), valeur −4 ;
- {j,p} et {j,q} avec Q_BV = 0 : x = (41,39,40,8,16), valeur −16 ;
- les mêmes avec Q_BV = 1 : x = (171,165,152,56,48), valeur −7824 ;
- {j,r} avec Q_BV = 1 : x = (147,141,136,40,48), valeur −3216 ;
- pour les cas en l : échanger les deux premières coordonnées.

## 6. Statut
- **Vérifié par « structures »** :
  - §§ 0 à 4 à la main ;
  - § 5 à la main, avec les formules de Schur ci-dessus, et en fractions exactes par
    `scripts/check_dossier4_compressions.py`, dont la sortie coïncide avec les deux tableaux.
- Concorde avec Astra (dossier 4, §§ 6–7) sur tous les points.
- **Portée.** Rien de tout cela ne ferme les carrés dis. Ce sont des filtres de colonnes pour le modèle unifié de C2 :
  supports de l'apex et de D(A), et au moins ⌈n_dis/3⌉ carrés dis sans apex intérieur.
