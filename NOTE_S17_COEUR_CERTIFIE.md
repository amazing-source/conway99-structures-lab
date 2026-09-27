# S17 — Cœur tout-tordu (1 triangle et sans triangle) : calcul certifié

> ⚠ **ERRATUM (27/09, « structures »).** Tout ce qui suit suppose que le graphe D décompose E(K₇) en cycles SIMPLES : les 121 215
> décompositions de `scripts/k7_cycle_decomp.py`. Ce n'est démontré nulle part : le carnet de C2 (§14) ne prouve que la
> 2-régularité du graphe D sur les paires. En tout-tordu, D est un **système de transitions** (en chaque bloc, un couplage des
> 6 paires en 3 transitions), et une chaîne de paires tordues peut repasser par un même bloc, en circuit non simple. Il y a 15⁷
> systèmes étiquetés, dont 121 215 seulement à circuits tous simples.
> **L'exclusion du « tout-tordu pur » ne vaut donc que pour les circuits simples.** Le cas non simple reste ouvert, sauf lemme
> de simplicité. Il est couvert si les 7 types locaux (0,6) (tier06) sont UNSAT.


> **État au 26/09/2026 (fin de matinée) : le tout-tordu pur est ENTIÈREMENT EXCLU.**
> - ≥ 2 triangles : à la main (NOTE_S8, S8.5, NOTE_S11).
> - 1 triangle : 15 classes, 84 certificats entiers au niveau paires avec (P2) (§2.1). Revérifiés par du code séparé ;
>   couverture de 41 790 décompositions.
> - 0 triangle : 14 classes, 14 certificats entiers au niveau ORBITES, à la racine et sans branchement (§2.3).
>   Revérifiés par du code séparé ; couverture de 30 360 décompositions.

26/09/2026, labo « structures ». **Calcul autorisé par l'utilisateur** (réponse « Calcul certifié plafonné », environ
1 h). Seuls les certificats exacts comptent.

## 1. Ce qui est démontré, et comment
**Énoncé visé.** Dans le cas involution de C2, en tout-tordu pur (toutes les orbites « un bloc », D = décomposition de
K₇ en cycles), aucune décomposition à un triangle ou sans triangle n'est réalisable. Avec S8, S8.5 et S11 (au moins
2 triangles, à la main), c'est **tout le tout-tordu**.

**Contraintes utilisées** (toutes démontrées) :
- système local C0–C8 de C2, résolu exactement par énumération (`local_system.py`) ;
- 6-cycle libre en chaque bloc : symétrie des arêtes, coupes de sous-tours de tailles 3, 4, 5 (valides pour un seul
  6-cycle) et exclusion de deux triangles ;
- PPE (S11, tiré de 16.1 de C2) : chaque paire est doublée par exactement 2 orbites ;
- couplage E[P][Q] = E[Q][P] pour les paires disjointes (symétrie de l'adjacence) ;
- **(P2) de C1** (notes/116 §2, vérifiée en S12), transposée au cadre de C2 (S14) : pour deux paires Q1 ≠ Q2,
  Σ_Z val_Z(Q1)val_Z(Q2) + Σ_{Z⊂Q1} val_Z(Q2) = 16 si elles sont disjointes, 12 si elles partagent un bloc.

**Modèle de configurations (S15).** Une colonne est un triplet (orbite, paire de voisins libres, solution locale
exacte). Toute valeur val_Z(Q) est alors une constante de la colonne, et (P2) devient linéaire.

**Certificat.** Des multiplicateurs entiers y, un par ligne et stockés avec l'étiquette de la ligne, tels que
Σ_c min(0, w_c)u_c > Σ_r max(y_r lo_r, y_r hi_r), avec w = Aᵀy et u_c ∈ {0,1} la boîte de la branche. Pour tout x
admissible, Σ_r y_r(Ax)_r = wᵀx serait à la fois ≥ le membre gauche et ≤ le membre droit : c'est absurde.

**Branchement.** Quand le PL d'un nœud est faisable, on sépare sur la paire de voisins d'une orbite, en prenant
**toutes** ses paires. Les feuilles couvrent donc tous les cas.

**Vérification indépendante** (`scripts/verify_core_certs.py`, entiers seuls, aucun PL) :
- les colonnes et la table des val sont refaites avec un code séparé (sa propre fonction de repère) ;
- chaque ligne du certificat est **reconstruite depuis son étiquette** ;
- l'inégalité est vérifiée exactement ;
- l'arbre de branchement doit être complet ;
- tout représentant doit couvrir, par une forme canonique propre, toutes les décompositions étiquetées.

**Contrôle négatif** (`scripts/verify_core_negative.py`) :
- un certificat privé de ses lignes (P2) est rejeté ;
- 4 corruptions sur 5 d'un multiplicateur sont rejetées. La 5ᵉ portait sur un multiplicateur minuscule devant la marge
  de 10⁹, ce qui est normal.

## 2. Résultats

### 2.1 Un triangle : IMPOSSIBLE (certifié, revérifié indépendamment)
- 15 classes non isomorphes, soit 41 790 décompositions étiquetées.
- Calcul : 937 s, RAM max 258 Mo (plafond 4000 Mo), `logs/certify_one.runlog`.

| classe | type | étiquetées | nœuds | certificats |
|---|---|---|---|---|
| one0 | (3,4,4,4,6) | 1260 | 11 | 10 |
| one1 | (3,4,4,4,6) | 5040 | 1 | 1 |
| one2 | (3,4,4,4,6) | 1680 | 1 | 1 |
| one3 | (3,4,4,5,5) | 630 | 8 | 7 |
| one4 | (3,4,4,5,5) | 1260 | 46 | 41 |
| one5 | (3,4,4,5,5) | 2520 | 1 | 1 |
| one6 | (3,4,4,5,5) | 5040 | 1 | 1 |
| one7 | (3,4,4,5,5) | 5040 | 1 | 1 |
| one8 | (3,4,4,5,5) | 1260 | 10 | 9 |
| one9 | (3,4,4,5,5) | 5040 | 1 | 1 |
| one10 | (3,4,7,7) | 1260 | 8 | 7 |
| one11 | (3,5,6,7) | 5040 | 1 | 1 |
| one12 | (3,5,6,7) | 2520 | 1 | 1 |
| one13 | (3,6,6,6) | 2520 | 1 | 1 |
| one14 | (3,6,6,6) | 1680 | 1 | 1 |

- Au total, 84 certificats (`logs/core_certs/one_*.json.gz`).
- Les classes les plus symétriques, étiquetées 630 ou 1260, demandent un branchement. Leur PL moyenne sur la symétrie.
- **Revérification indépendante** (`logs/verify_core_one.log`) :
  - les 84 certificats sont valides ;
  - les 15 arbres sont complets ;
  - **couverture : les 41 790 décompositions étiquetées à 1 triangle sont toutes isomorphes à un représentant
    certifié.**
- **Conclusion : aucune configuration tout-tordue à un seul triangle n'existe.** C'est une preuve assistée par
  ordinateur, vérifiable en entiers.

### 2.2 Sans triangle : non résolu par (P2) ; voir §2.3 pour le niveau orbite
- 14 classes non isomorphes, soit 30 360 décompositions étiquetées :
  - (4,4,4,4,5) : 1 classe ;
  - (4,4,6,7) : 2 ;
  - (4,5,5,7) : 3 ;
  - (4,5,6,6) : 3 ;
  - (5,5,5,6) : 3 ;
  - (7,7,7) : 2.
- **À la racine, avec (P2), les 14 PL sont FAISABLES** (`logs/certify_zero_root.runlog`, points intérieurs, 215 s). La
  version typée (types signés des doubles) reste aussi faisable sur zero0 et zero5.
- **Explication.** Sans triangle, tous les systèmes locaux sont identiques (C6, C8 seulement). Le PL admet des
  solutions fractionnaires symétriques, et seul un branchement les casse. Les classes à 1 triangle qui ont demandé un
  branchement étaient précisément les plus symétriques.
- **Coût.** Chaque nœud coûte environ 8 à 15 s. Avec une profondeur probable de 2 ou 3, il faudrait des centaines à des
  milliers de nœuds par classe, soit des heures au total. C'est au-delà de l'heure autorisée : **décision de
  l'utilisateur requise**, ou bien de nouveaux ingrédients mathématiques.
- **Pistes :**
  - (P2) par ORBITES plutôt que par paires (côté invariant Q_out de C2) : **suivie, voir §2.3** ;
  - le côté signé S_out ;
  - la reformulation en repère ajusté : UᵀU = 16,5 I + 16 J − 4A_T(7), avec pour tout 4-cycle de blocs
    Σ_Z ⟨u_Z, g⟩² = 98.
- Note technique : le simplexe dual de HiGHS n'aboutit pas en 110 à 300 s sur ces PL (12 600 colonnes,
  844 000 non-zéros), alors que les points intérieurs mettent 8 s.

### 2.3 Identité de second ordre au niveau des ORBITES (S19, 26/09) : la racine casse
**Énoncé.** Pour toutes orbites U, W (42 × 42), avec q_U(W) le nombre de voisins d'un sommet de U dans W :
  Σ_Z q_U(Z) q_Z(W) + q_U(W) = 12 [U = W] + 4 − |blocs(U) ∩ blocs(W)|,
soit **Q² + Q = 12 I + 4 J − Sh**, où Sh = Σ_i b_i b_iᵀ.

*Preuve.*
- Soient B l'adjacence lointaine (84 × 84) et N l'incidence lointain–étiquette. A² = 12I + 2J − A, restreint aux
  lointains, donne B² + B = 12I + 2J − NNᵀ.
- σ commute avec B. Avec S les indicatrices d'orbites (84 × 42), on a BS = SQ. Puisque σ échange les deux étiquettes
  de chaque bloc (note 10 de C2 : sept orbites intérieures de signe −1), (NNᵀS)(z, W) = Σ_{w∈W} |lab(z) ∩ lab(w)|, ce
  qui vaut le nombre de blocs communs.
- SᵀS = 2I conclut. ∎

**Faits annexes, généraux dans le cadre de C2.**
- q_U(U) = 0 : les voisins communs de z et σz vont par paires {w, σw}, donc z ~ σz contredirait λ = 1.
- Le seul 2 d'une ligne est le partenaire D : μ(z, σz) = 2 est épuisé par le carré. Un double lointain touche donc les
  deux orbites de sa paire.
- Spectre : 12 sur 𝟏, −2 sur les b_i − b_j (dimension 6, note 10 §5 de C2), puis 3 (×20) et −4 (×15) sur les 35
  dimensions restantes (trace Q = 0).
- (P2) par paires en est la compression : la partie symétrique, de dimension 14, qui est le repère ajusté E₋₂ de C1
  (notes/117 §2 bis). **Le niveau orbite ajoute la partie antisymétrique**, de dimension 21 (différence entre les deux
  orbites d'une paire), ainsi que les termes croisés.
- Sans symétrie, la même identité vaut pour les 42 « cellules antipodales » {(c_s,d_t), (c_{1−s},d_{1−t})}, avec une
  colonne par lointain. Transmis à C1 pour son cadre.

**Modèle** (`scripts/orbit_p2_colgen.py`).
- Une colonne-orbite est une colonne par paires, plus le choix, pour chaque paire lointaine de valeur 1, de l'orbite
  qui reçoit le voisin. Un double touche les deux orbites, et le carré va à D(Z). Cela donne 2^m colonnes-orbites par
  colonne par paires, soit **352 800** pour chaque classe sans triangle.
- Lignes :
  - une unité par orbite ;
  - les 6-cycles libres (symétrie et coupes de sous-tours 3, 4, 5) ;
  - PPE ;
  - couplage q_U(W) = q_W(U) pour les paires disjointes ;
  - les 903 identités ci-dessus.
- Phase 1 par génération de colonnes. Le pricing est exact : on énumère les 2^m choix.
- **Certificat** : multiplicateurs entiers (lignes P2o doublées, échelle 10⁹) tels que
  Σ_{r ≠ O} min(y lo, y hi) > Σ_Z max_{c ∈ Z} (Aᵀy)_c.
  Le maximum est pris sur **toutes** les colonnes-orbites.

**Résultats.** Les 14 classes sont certifiées à la racine, sans branchement. Pour chacune, le PL par paires avec
(P2) était faisable à la racine (§2.2). Ce sont les tests de 5 min, puis le lot autorisé par l'utilisateur (« Lancer
les 13 classes », `scripts/run_orbit_batch.sh`, une classe à la fois, moins de 0,8 Go de RAM).
| classe | type | étiquetées | tour | phase 1 à l'arrêt | écart low − high (unités de y) |
|---|---|---|---|---|---|
| zero0 | (4,4,4,4,5) | 2520 | 8 | 34,30 | 26,32 |
| zero1 | (4,4,6,7) | 2520 | 8 | 17,90 | 9,38 |
| zero2 | (4,4,6,7) | 2520 | 7 | 29,04 | 12,80 |
| zero3 | (4,5,5,7) | 2520 | 6 | 36,73 | 2,43 |
| zero4 | (4,5,5,7) | 2520 | 6 | 34,97 | 0,36 |
| zero5 | (4,5,5,7) | 2520 | 6 | 48,27 | 16,34 |
| zero6 | (4,5,6,6) | 5040 | 7 | 32,94 | 19,11 |
| zero7 | (4,5,6,6) | 2520 | 6 | 35,71 | 2,12 |
| zero8 | (4,5,6,6) | 2520 | 6 | 42,66 | 6,56 |
| zero9 | (5,5,5,6) | 840 | 6 | 50,72 | 20,70 |
| zero10 | (5,5,5,6) | 2520 | 6 | 45,20 | 9,93 |
| zero11 | (5,5,5,6) | 840 | 6 | 56,90 | 31,38 |
| zero12 | (7,7,7) | 840 | 7 | 33,94 | 18,16 |
| zero13 | (7,7,7) | 120 | 5 | 77,77 | 21,43 |

- Total : 30360 décompositions étiquetées.
- zero13 : le test flottant initial (convergé, phase 1 = 69,68) a été refait dans le lot pour extraire le certificat.

**Revérification indépendante et couverture** (`python scripts/verify_orbit_certs.py coverage`,
`logs/verify_orbit_certs.log`) :
- les 14 certificats sont VALIDES : bornes recalculées identiques, et mêmes low et high que le générateur à l'unité
  près ;
- couverture : les **30 360** décompositions étiquetées sans triangle sont toutes isomorphes à un représentant
  certifié (forme canonique propre `canon_any`), et les effectifs s'additionnent exactement ;
- verdict : **ALL VERIFIED**.

**Contre-vérification par C2, avec son propre code** (`c2_experimental/scripts/xcheck_structures_orbit_certs.py`,
`xcheck_structures_coverage.py` ; carnet C2 §32). C2 n'a lu que les fichiers de certificats et la sémantique documentée.
- Colonnes reconstruites avec le système local de C2 (sommes par bloc 4 − [t∈p] − 2[t=k]) : 352 800 colonnes-orbites
  par classe.
- Lignes reconstruites depuis les étiquettes, avec ses propres bornes : 0 écart.
- 14/14 VALIDES, avec low et high identiques.
- Couverture : réénumération propre, qui donne 30 360 décompositions sans triangle. Les 14 représentants sont deux à
  deux non isomorphes, et Σ 5040/|Aut| = 30 360.
- Identité redérivée : Q² + Q = (k − 2)I + 4J − Sh.
**Le tout-tordu pur est donc exclu, avec double vérification indépendante.**

**Vérification de zero0** (`scripts/verify_orbit_certs.py`) :
- Le vérificateur est un code séparé, en entiers seuls. Les colonnes par paires viennent de `my_columns` (propre
  fonction de cadre), et le vérificateur fait lui-même l'expansion en orbites. Chaque ligne est reconstruite depuis
  son étiquette, bornes comprises.
- Il retrouve **exactement** low = 37 079 589 870 et high = 10 760 944 756, soit un écart de 26,3 à l'échelle 10⁹,
  sur 352 800 colonnes et 42 orbites : VALIDE.
- Incident de format, corrigé : le premier fichier stockait les bornes P2o non doublées. Le vérificateur, qui ne les
  utilise pas, les a signalées. Le générateur stocke désormais les bornes des lignes doublées, et le fichier de zero0 a
  été corrigé de la même façon.

**Contrôle sur un vrai graphe.** L'identité d'orbites est exacte sur BvLS, srg(243,22,1,2), avec σ(v) = −v :
Q² + Q = 20I + 4J − Sh. Q y a une diagonale nulle et un seul 2 par ligne (`scripts/orbit_identity_real.py`).

## 3. Portée et réserves
- La preuve est **assistée par ordinateur**. Chaque certificat est un objet fini, vérifiable en entiers ; il ne dépend
  pas du solveur qui l'a trouvé.
- Deux dépendances sont partagées entre le générateur et le vérificateur :
  - le solveur local exact `solve_core`, lui-même recoupé avec doubles_budget3 sur 1340 clés (S11) ;
  - la liste des décompositions de K₇ (`zeros_budget.decompositions`, dont les effectifs par type ont été contrôlés,
    121 215 au total).
