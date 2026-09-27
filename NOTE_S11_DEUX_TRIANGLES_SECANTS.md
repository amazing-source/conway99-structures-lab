# S11 — Tout-tordu à deux triangles sécants : impossible (preuves à la main)

> ⚠ **ERRATUM (27/09, « structures »).** Tout ce qui suit suppose que le graphe D décompose E(K₇) en cycles SIMPLES : les 121 215
> décompositions de `scripts/k7_cycle_decomp.py`. Ce n'est démontré nulle part : le carnet de C2 (§14) ne prouve que la
> 2-régularité du graphe D sur les paires. En tout-tordu, D est un **système de transitions** (en chaque bloc, un couplage des
> 6 paires en 3 transitions), et une chaîne de paires tordues peut repasser par un même bloc, en circuit non simple. Il y a 15⁷
> systèmes étiquetés, dont 121 215 seulement à circuits tous simples.
> **L'exclusion du « tout-tordu pur » ne vaut donc que pour les circuits simples.** Le cas non simple reste ouvert, sauf lemme
> de simplicité. Il est couvert si les 7 types locaux (0,6) (tier06) sont UNSAT.


25/09/2026, labo « structures », pour C2. Cadre et notations de NOTE_S8 : système local C0–C8 par orbite, CL (S8.4),
lemme D, budget 16.1 de C2.

**Résultat.** Aucune décomposition tout-tordue de K₇ à deux triangles sécants n'est réalisable. Il y en a 27 720
étiquetées, en 8 classes. Avec S8.5 (triangles disjoints) et §5 de NOTE_S8 (au moins 3 triangles), **tous les types
tout-tordus à au moins deux triangles sont exclus à la main.** Restent les types à un triangle (41 790 étiquetées) et
sans triangle (30 360).

## 1. Deux outils globaux

**Budget exact paire par paire (PPE).**
- Soit P′ une paire. Parmi les 20 orbites des paires disjointes de P′, **exactement 2** ont φ(P′) = 2, 12 ont
  φ(P′) = 1 et 6 ont φ(P′) = 0.
- *Preuve.*
  - Dans 16.1 de C2, on a n_s(P′) = 0. Un partenaire D w de P′ = {i,k} (appariée par i) a exactement v et σv comme
    voisins touchant le bloc i (lemme 3), donc φ_w(P′) = 2. Son terme C(2,2) − 1 est nul.
  - Les autres étiquettes ont φ ≤ 2. Donc exactement 4 étiquettes, soit 2 orbites, ont φ = 2.
  - Ce sont des orbites de paires disjointes de P′. Une paire proche a φ = n_r ≤ 1 (bloc libre) ou φ = 0 (F1), hormis
    le carré.
  - Chaque étiquette de P′ a exactement 8 voisins dans des paires disjointes : 12 extérieurs, moins 4 proches (m, π de
    ses deux blocs, sans sœur adjacente puisque n_s = 0). La somme vaut donc 32 étiquettes, soit 16 orbites, d'où
    12 uns et 6 zéros.
- Relu par C2 (25/09), qui l'étend aux mélanges.
- Mon calcul exhaustif (NOTE_S8 §3) n'utilisait que la somme Σ ≤ 42. **Le total vaut exactement 42**, réparti à raison
  de 2 par paire.

**Couplage.** Pour deux paires disjointes P, Q, on pose E[P][Q] = φ_X(Q) + φ_{X′}(Q), où X, X′ sont les deux orbites de P.
- Comme e(P,Q) = 2E[P][Q] compte les arêtes entre les 4 étiquettes de P et les 4 de Q, on a **E[P][Q] = E[Q][P]**.
- Ce couplage relie deux systèmes locaux différents, que le calcul exhaustif traitait séparément.

## 2. Lemmes locaux

**Lemme T (triangle dans R).** Si R_x = T ∪ {r} avec T un triangle D, alors **y_r = 1 − n_k − n_r**. En particulier
n_k + n_r ≤ 1, et la paire {k,r} ne porte jamais de double.
*Preuve.*
- On somme C8 sur t ∈ T et on utilise C5 : N_T + Y_T + Φ_r = 6.
- C8(r) donne Φ_r = 4 − n_r − y_r, et C6 donne Y_T = 2 − n_k − y_r.
- Avec N_T = 2 − n_k − n_r, on obtient 2y_r = 2 − 2n_k − 2n_r. ∎

Contrôle : les 430 cas des clés exhaustives.

**Lemme A (deux triangles T1 = {v,a1,a2}, T2 = {v,b1,b2} sécants en v ; e, f les deux autres blocs).**
(a) Les deux orbites de {a1,a2} ont φ(b1b2) = 2, et les deux orbites de {b1,b2} ont φ(a1a2) = 2. Par PPE, ce sont
    les seuls doubles de ces deux paires.
(b) Au bloc v, les orbites e et f ont chacune exactement un voisin dans {a1,a2} et un dans {b1,b2}.

*Preuve.*
- Pour une orbite Y de {a1,a2}, on a k = v et n_v = 0 (C0). C4 sur T2 et C6 donnent φ(b1b2) = 1 + y_e + y_f.
- En sommant sur Y, Y′ : 2 + E[a1a2][ve] + E[a1a2][vf].
- Or E[ve][a1a2] = 1 + (1 − N_{a1a2}) :
  - C1 pour l'orbite appariée par v ;
  - C2 ou C3 pour l'orbite libre en v, avec N ≤ 1.
- La somme vaut donc 6 − N − N′ ≤ 4, ce qui force N = N′ = 1 et φ = 2. ∎

**Corollaire A1.** Soit C_ef = (…, p, e, f, q, …).
- L'orbite W de {e,f} libre en f a k = p. CL donne φ_W = 3 − N sur la paire de base du triangle qui ne contient pas p.
  Par (a), N = 2 : au bloc f, les voisins de e sont exactement cette paire de base.
- Même chose au bloc e, avec q.

**Pas de triangle dans un 6-cycle libre.** Si a ~ b et a ~ c au bloc j, alors b ≁ c. Cela transforme « N ≥ 1 »
en « N = 1 », donc en double par CL.

## 3. Les 8 classes

Étiquettes : v = 0, T1 = {0,1,2}, T2 = {0,3,4}, e = 5, f = 6. Notation (j;i) : orbite de {i,j} libre en j.

| classe | cycles hors triangles | argument final |
|---|---|---|
| (3²,4,4,7) | (1,3,2,5,0,6,4), (1,5,3,6), (2,6,5,4) | E[46][35] ≥ 2+1 > 0+2 ≥ E[35][46] |
| (3²,4,5,6)a | (3,1,4,6,0,5), (3,6,1,5,2), (4,5,6,2) | (3;5) : φ12 = 2 − n4 + y6, donc 5 ~ 4 au bloc 3 ; puis E[36][45] = 4 > 3 |
| (3²,4,5,6)b | (1,5,2,6), (1,3,5,0,6,4), (2,3,6,5,4) | 10 doubles forcés sur les 4 paires croisées, contre 8 |
| (3²,4,5,6)c | (1,5,0,6,3), (1,6,5,3,2,4), (2,5,4,6) | E[36][45] ≥ 2+1 > 0+2 |
| (3²,4,5,6)d | (1,5,0,6,3), (1,6,2,3,5,4), (2,5,6,4) | E[56][34] = 2 donne 4 ~ {5,6} au bloc 3 ; puis E[36][45] = 4 > 2 |
| (3²,4,5,6)e | (1,5,0,6), (1,3,2,5,6,4), (2,6,3,5,4) | E[26][15] ≥ 2+1 > 0+2 |
| (3²,4,5,6)f | (1,5,0,6), (1,3,5,2,6,4), (2,3,6,5,4) | saturations : un troisième double en 36 ou 45 |
| (3²,5³) | (3,6,0,5,1), (4,2,5,6,1), (3,5,4,6,2) | E[36][45] ≥ 2+1 > 1+1 |

**Modèle détaillé, (3²,5³).**
1. Lemme A : {3,4} est saturée par les orbites de {1,2}.
2. W = (6;5), avec k = 2 ∈ T1. CL donne φ34 = 3 − N34 ≤ 1, d'où N34 = 2 : au bloc 6, 5 ~ {3,4}.
3. (6;3), avec k = 1 : CL donne φ45 = 3 − N45. Au bloc 6, 3 ~ 5 et 3 ≁ 4 (pas de triangle), d'où φ45 = 2.
4. (3;6), avec k = 0 : CL donne φ45 ≥ 1. Donc E[36][45] ≥ 3.
5. (4;5), avec k = 3 et R = T1 ∪ {6} : lemme T, φ36 = y6 ≤ 1. (5;4), avec k = 6 et R = T1 ∪ {3} : φ36 = y3 ≤ 1.
   Donc E[45][36] ≤ 2. Contradiction avec le couplage. ∎

**Classe b, la plus longue.**
- Couplage {1,2}–{5,6} : E[56][12] = 1 + 1 = 2. Donc N56 = 2 pour les orbites de {1,2}.
- Couplage {3,4}–{5,6} : E[56][34] = y3 + y4 = 1 + 1 (lemme T pour W, W′). Donc 4 ~ {5,6} au bloc 3 et 3 ~ {5,6} au
  bloc 4.
- Au bloc 1, le cycle est (2,5,a,0,b,6) avec {a,b} = {3,4}, par C0, C3 et le lemme T. De même au bloc 2.
- (1;3) : la somme C8(2) + C8(4), avec C8(6), C8(0) et C6, donne φ24 = 2 + y6, donc φ24 = 2. De même :
  - (1;4) double en 23 ;
  - (2;3) double en 14 ;
  - (2;4) double en 13.
- Six orbites vérifient φ_a + φ_b = 3 sur deux paires croisées, donc un double chacune. C'est C8 avec y = 0 par le
  lemme T :
  - (3;6) : 14 + 24 ;
  - (4;5) : 13 + 23 ;
  - (5;1), (6;1) : 23 + 24 ;
  - (5;2), (6;2) : 13 + 14.
- Cela fait 10 doubles sur {13, 14, 23, 24}, alors que PPE en impose exactement 8. ∎

**Classe f.**
- Saturations par CL :
  - {4,6} par (5;3) et (3;5) ;
  - {3,5} par (6;4) et (4;6) ;
  - {3,6} par (5;4) ;
  - {4,5} par (6;3).
- Au bloc 0, la structure est (1,x,s,2,x′,s′), par le lemme A et C0/C2. C8(3) + C8(4) + C8(6) donne :
  - φ36 + φ46 = 3 pour (0;5), d'où φ36 = 2 ;
  - φ35 + φ45 = 3 pour (0;6), d'où φ45 = 2.
- Pour (0;1), la somme croisée vaut 6 − φ56 ≥ 5. Avec φ35, φ46 ≤ 1, il faut un double en 36 ou 45, donc un troisième. ∎

## 4. Contrôles, et rôle du calcul (tests courts, règle du 25/09)

**Preuves.** Elles sont à la main. `scripts/verify_two_sharing_hand.py` vérifie chaque étape locale avec le solveur
exact du système local (`scripts/local_system.py`), pour **toutes** les paires de voisins permises par les faits
invoqués. Il vérifie aussi le lemme T sur les 430 cas des clés exhaustives (`logs/verify_two_sharing_hand.log`).

**Tests, pour orienter les maths.** Chacun a duré moins de 5 minutes, et aucun n'est utilisé comme preuve.
- `scripts/ppe_test.py` (PLNE, scipy/HiGHS) : les 8 classes sont infaisables avec PPE et couplage, et faisables sans.
  - Contrôles : triangle de triangles, Fano et triangles disjoints infaisables.
  - Variantes (`ppe_variants.py`) : PPE seul ne suffit pas, le couplage seul non plus (sauf pour la classe a).
- Le relâchement linéaire est déjà infaisable (`ppe_lp_relax.py`). Les certificats de Farkas creux
  (`ppe_certificate.py`) et les filtres de suppression (`ppe_minimal_rows*.py`, `ppe_single_extra.py`,
  `ppe_filter_protected.py`, `ppe_cross_bounds.py`) ont indiqué quelles lignes lire. Les arguments ci-dessus en sont la
  lecture humaine, souvent plus courte.
- **Limite mesurée** (`ppe_sample_rest.py`) : pour les types à 1 triangle, des échantillons restent **faisables** avec
  PPE, couplage et système local. Ces outils ne suffiront donc pas pour le cœur.
