# NOTE S25 : la famille srg(2m²+1, 2m, 1, 2) vue depuis ses membres connus (9 et 243)

Demande de l'utilisateur (26/09) : étudier à fond les graphes à 9 et 243 sommets, chercher une règle vraie pour eux
qui ne peut pas valoir à 99, et essayer de trouver un nouveau membre. Code dans `famille/`. Recherche
bibliographique en cours : unicité du graphe à 243 sommets et méthode de preuve (§6).

## 1. Arithmétique de la famille
- **Paramètres.** λ = 1 et μ = 2 imposent v = (k² + 2)/2. Avec k = 2m, les valeurs propres sont r, s = (−1 ± t)/2,
  où t² = 8m − 7. Les multiplicités sont entières seulement si **t | 63** : m ∈ {2, 7, 11, 56, 497}, soit
  v ∈ {9, 99, 243, 6 273, 494 019}.
- **Les membres qui existent ont t ∈ {3, 9}** (des puissances de 3, et v est une puissance de 3). **Les trois cas
  ouverts ont t ∈ {7, 21, 63}, c'est-à-dire exactement 7 | t.**
  - Pour p | t (p impair), v ≡ 81/32 (mod p). Donc p | v ⟺ p = 3. Le côté 3-adique (p | v) est celui des graphes
    connus ; le côté 7-adique (p ∤ v) n'existe que dans les cas ouverts.
  - Conséquence 1 : un argument 7-adique qui n'utilise que « 7 | t, 7 ∤ v » et les identités de la famille exclurait
    d'un coup 99, 6 273 et 494 019.
  - Conséquence 2 : aucun membre réel de la famille ne peut servir de contrôle au côté 7-adique. Ces contrôles
    doivent venir d'autres familles (NOTE_S5, registre S1).

## 2. Ce que sont les deux graphes connus (construits et vérifiés, `famille/family_graphs.py`)
- **m = 2 : H(2,3)**, la grille 3×3, c'est-à-dire Cay(F₃², {±e₁, ±e₂}).
- **m = 11 : Berlekamp–van Lint–Seidel.** C'est Cay(F₃⁵, {±h_i}), où les h_i sont les colonnes d'une matrice de
  contrôle du code de Golay ternaire [11,6,5]. Le code (résidus quadratiques mod 11) est reconstruit et vérifié :
  729 mots, distance 5.
- **Mécanisme d'existence.** Le code est parfait, donc chaque classe a un unique représentant de poids ≤ 2, ce qui
  donne λ = 1 et μ = 2 en une ligne. Dans les deux cas, le graphe enraciné en x est la boule de Hamming B₂(0) ⊂ F₃^m,
  et elle se referme en groupe grâce à un code parfait.

## 3. Leurs invariants dans le langage de C1 et C2 (`famille/family_invariants.py`, `famille/control_c1_note125.py`)

| invariant | H(2,3) | BvLS | à 99 |
|---|---|---|---|
| holonomie des droites | triviale partout | triviale partout | tordue en tout sommet (L8) |
| P (prismes) | 6 = P_max | 8 910 = P_max | ≤ 1 151 < 1 386 = P_max (L16, Krein) |
| points de Hamming | 9/9 | 243/243 | aucun (L8) |
| S (événements E ~ Ē) | 0 | 0 | libre, ≤ 4 158 − 3P |
| Σ_L c₄(G(L)) = S + 3P | 18 | 26 730 | identité de C1 (notes/125), contrôlée ici |
| rang₃(A − rI) | 4 | 67 | sans objet (le côté est 7-adique) |
| Smith 3-adique de A − rI | 1⁴ 3¹ 0⁴ | 1⁶⁷ 3³⁴ 9¹⁰ 0¹³² | — |

## 4. Les « règles » vraies pour 9 et 243 qui échouent à 99 (toutes déjà prouvées)
- **Être une puissance de 3 / un graphe de classes d'un code parfait.** À 99, c'est impossible : 99 n'est pas une
  puissance de 3, et 3 n'est pas une valeur propre de H(7,3) (C2 B01 ; C1 L1).
- **Tout-sœur** (cadre de C2, avec involution : toutes les grilles par le point fixe sont complètes). Exclu à 99
  (C2 : DRAT et lemme des étoiles).
- **Un point de Hamming.** Aucun à 99 (C1 L8).
- **P = P_max (grid-full).** Exclu à 99 par la condition de Krein pour les droites (C1 L16). Dans la famille, elle
  n'exclut grid-full qu'à 99, pas à 6 273 ni à 494 019 (C1, notes/70).
- **Vertex-transitivité / Cayley.** Exclu à 99 (|Aut| ≤ 3, littérature).

**Pourquoi cela ne suffit pas.** Toutes ces règles disent la même chose : Γ₉₉ serait fortement tordu, donc très
loin de 9 et 243. Les membres connus vivent à l'extrême opposé (torsion nulle) et ne disent rien du régime où vivrait
Γ₉₉. La contradiction demande une règle qui interdise la torsion forte : un « plafond » (notes/80 et 100 de C1).
Aujourd'hui, seule l'arithmétique 7-adique pousse dans ce sens : un Γ 7-générique aurait k = 30, exclu par L22.

## 5. Nouveau : la famille est rigide pour la commutation de Godsil–McKay
**Lemme.** Soit Γ = srg(v, k, 1, 2) avec k ≥ 8. Γ n'a aucun ensemble de commutation de Godsil–McKay non trivial.
Un tel ensemble serait C ⊂ V avec |C| = 2n ≥ 4, C induisant un graphe d-régulier, et chaque sommet hors de C ayant
0, n ou 2n voisins dans C.

*Preuve.*
- **Comptage autour de a ∈ C.** Ses k − d voisins extérieurs ont chacun au moins n − 1 autres voisins dans C. Pour
  x ∈ C ∖ {a}, a et x ont au plus 1 voisin commun (s'ils sont adjacents) ou 2 (sinon). D'où
  (k − d)(n − 1) ≤ d + 2(2n − 1 − d), soit (k − 4)n − (k − 2) ≤ d(n − 2).
  - Pour n = 2, cela donne k ≤ 6, impossible.
  - Pour n ≥ 3, cela donne d ≥ (k − 4) + (k − 6)/(n − 2), donc d ≥ k − 3 : chaque sommet de C a au plus 3 voisins
    hors de C.
- **Borne spectrale.** e(C, C̄) ≥ (k − r)|C||C̄|/v, où r est la seconde valeur propre. On en tire
  |C̄| ≤ 3v/(k − r), puis n ≥ v(1 − 3/(k − r))/2.
- **Conclusion.** Cette borne vaut 36 pour 99, plus de 101 pour 243, et elle est bien plus grande que k pour 6 273 et
  494 019. Un sommet extérieur voisin de C aurait au moins n > k voisins dans C, ce qui est impossible. C est donc
  une réunion de composantes, et C = V. ∎

**Portée.** On ne peut pas fabriquer un autre srg(243,22,1,2), ni un srg(99,14,1,2), par commutation locale de
Godsil–McKay (à une seule classe) à partir d'un graphe de la famille. Ce n'est pas une preuve d'unicité.

## 6. Suite
- **Unicité de BvLS, et surtout sa méthode.** Si une preuve force « tout sommet est un point de Hamming » à m = 11,
  il faut voir si l'argument se transpose à m = 7, où L8 l'interdit. Ce serait exactement le plafond qui manque. S'il
  utilise 3 | v, il ne se transpose pas, et on saura pourquoi.
- **Un membre tordu serait l'objet le plus utile.** Soit un autre srg(243,22,1,2), soit un « faux monde » fermé aux
  paramètres locaux de 99. Ce serait le premier contrôle positif du régime tordu. Le §5 ferme la voie la plus simple
  (commutation locale).
