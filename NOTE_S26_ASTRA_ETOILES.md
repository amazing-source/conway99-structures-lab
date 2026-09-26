# NOTE S26 : relecture du rapport d'Astra « étoiles d'ombres » (suite du dossier n° 5, 26/09)

Rapport copié tel quel dans `astra/RAPPORT_ASTRA_DOSSIER5_SUITE.md`. Le cadre et les notations sont ceux de NOTE_S24.

## §1. Famille cyclique exclue à k = 22 : VÉRIFIÉ par un code indépendant
Scripts `scripts/star_shadows.py` (test d'étoile réutilisable) et `scripts/verify_astra_d5_cyclic.py`. Log :
`logs/verify_astra_d5_cyclic.log`.
- **Famille.** Les 100 triplets (r, s, q) donnent des boules sans triangle, dont les réunions de deux parties n'ont
  aucun 4-cycle.
- **Filtre de L23.** Les triplets génériques avec τ non carré mod 7 (seuls possibles à k = 22) sont exactement les 6
  triplets {r, s} = {2, 4}, q ∈ {1, 3, 5}. Le τ est calculé exactement par Bareiss pour les 100.
  - τ(2,4,1) = 3 096 828 960 300 ≡ 5 et τ(2,4,3) = 2 905 653 667 500 ≡ 3, comme chez Astra.
- **Recherche sur 7 groupes.** Sur les groupes d'Astra, on retrouve les mêmes candidats (381, 336, 354, 336, 381,
  362, 362 et 375, 348, 405, 348, 375, 392, 351). Aucune solution.
- **Étoile complète** (10 groupes, capacité avec égalité) sur **les 6 triplets**, sans utiliser la réduction à deux
  représentants : aucune solution, entre 20 064 et 22 829 nœuds, environ 10 s par triplet.
- **Réduction par Fourier** : cohérente. Les racines 6ᵉ de l'unité sont dans F₇, les facteurs det H(ω)² sont des
  carrés, et 2 est un carré mod 7. Donc χ(τ) = χ(det H(−1)).

## §2. À k = 24 avec τ carré : |J_bad| ≥ 8 : DÉMONTRÉ, refait à la main (un point à confirmer par C1)
- **Mécanique.** R = uuᵀ + vvᵀ et −1 non carré mod 7. Donc R_yy = u_y² + v_y² = 0 ⟺ u_y = v_y = 0, et la ligne y de
  R est alors nulle. On a aussi R_yy = Z_yy − 6.
- **Noyau.** MR = 0 (automatique, NOTE_S24) donne Mu = Mv = 0, avec u et v portés par J_bad.
- **Indépendance de 5 colonnes de M.** Sur F₇, une dépendance entre au plus 5 arêtes de l'incidence non orientée
  est un C₄ alterné. Les circuits du matroïde des cycles pairs de taille ≤ 5 sont les seuls C₄ (un C₆ ou deux
  triangles liés demandent au moins 6 arêtes).
  - Elle se relève aux entiers en un échange n₁ + n₃ = n₂ + n₄, car les entrées sont dans [−2, 2].
  - Alors, dans chaque partie, n₁ rencontre n₂ (sinon n₁ et n₂, ou n₁ et n₄, partagent deux points dans une partie
    et un point en A). Cela fait 3 voisins communs, contre μ = 2.
- **Sept ombres.** 7 arêtes ont une nullité ≤ 2. Nullité 2 impose un C₄ :
  - biparti de nombre cyclomatique 2 avec ≤ 7 arêtes : K₂,₃ ou le thêta (1,3,3) ;
  - non biparti de nombre cyclomatique 3 : n ≤ 5, avec au moins n + 2 arêtes, or ex(5, C₄) = 6.
  Les noyaux coïncident, d'où l'échange interdit. Donc **|J_bad| ≥ 8**.
- **Point à confirmer par C1 (propriétaire de L23).** « τ carré ⇒ la partie non dégénérée de R est anisotrope ». Cela
  suppose que χ(disc R) = χ(τ)·ε(k), avec ε(22) = −1 (c'est L23) et ε(24) = +1. Le motif ε(k) = χ(−1)^{k/2} est
  cohérent avec les deux affirmations d'Astra (L23 à k = 22, et anisotrope ou scindé à k = 24), mais je n'ai pas
  refait la comptabilité des signes.

## §3. Composantes 4+4+4 des droites dégénérées : VÉRIFIÉ
Script `scripts/verify_astra_d5_444_components.py`.
- **Admissibles.** 66 données normalisées sur 243 (couplages internes distincts sur A–B et A–C, πm_Bπ⁻¹ ≠ m_C, π sans
  point fixe).
- **Arbres couvrants.** τ = 8 640 pour 54 données et 8 100 pour 12 : ≡ 2 et 1 mod 7, deux carrés non nuls. En graphes
  abstraits (sans les couches), on obtient 2 classes d'isomorphisme ; les types II et III d'Astra ont le même τ.
- **Conséquence.** d(L) = 3 pour une boule formée de trois telles composantes, et k ≥ 20 par L21.
- **Lien avec C1 (notes/125).** Par le lemme de Klein, ces composantes ont forcément des 4-cycles de type (ii),
  c'est-à-dire des événements S. Elles existent donc à P = 0 avec S > 0 (Astra), mais sont exclues à P = S = 0 (C1).
  Les deux résultats sont compatibles.

## Bilan
Le rapport est correct. §1 est une vraie exclusion, mais pour une famille de 6 boules seulement. §2 est une
contrainte, pas une contradiction. §3 est une classification. L'exclusion uniforme de k = 22 reste ouverte.
