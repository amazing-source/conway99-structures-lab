# S8 — Cas « tout-tordu pur » de C2 : le type (3,3,3,3,3,6) est impossible (preuve à la main) ; bornes pour les autres types

*25/09, labo « structures », en appui de C2 (branche involution). Statut : **DÉRIVÉ** pour S8.1–S8.3 (preuves à la
main) ; **ÉTABLI par calcul exhaustif contrôlé** pour le tableau du §3 (système local exact, solveur indépendant au
niveau des étiquettes, 0 désaccord). Toutes les références « § » renvoient au carnet de C2 (`notes/MATH_C2_PROGRAMME.md`, branche
math-c2-programme).*

## 0. Cadre
Γ = srg(99,14,1,2), σ involution à un seul point fixe x ; blocs B_i = {(i,0),(i,1)} ; étiquettes {a,b} ; 42 orbites.
**Tout-tordu pur** : Σ = ∅ (aucune paire frère), et chaque paire ij est tordue : une orbite appariée par i, l'autre par j.
Le graphe D induit alors une **décomposition de E(K₇) en cycles** (§14).
- Le cycle contenant l'arête ij passe par k – i – j ; l'orbite de ij appariée par i l'est vers la paire ik.
- Il y a 121 215 décompositions étiquetées, en 20 types de longueurs (`scripts/k7_cycle_decomp.py`). Recomptage
  indépendant : Fano = 30, (7,7,7) = 960, (3⁵,6) = 420 (`scripts/check_k7_counts.py`).

## 1. Faits utilisés (tous issus de λ = 1, μ = 2 et des lemmes de C2)
Soit X ⊂ ij appariée par i à D(X) ⊂ ik, libre en j ; x = {a,b} ∈ X ; R_x = [7] ∖ {i,j,k}.
- **(F0) n_s = 0.** Une arête sœur serait à la fois dans R_i et dans R_j (§11). Or en i l'orbite appariée par i n'a que
  sa digone, et en j c'est l'autre orbite.
- **(F1)** Les étiquettes voisines de x qui rencontrent le bloc i sont exactement v et σv, dans la paire ik (Lemme 3 :
  un seul voisin contient a, un seul contient a′, et le carré D les fournit).
- **(F2)** En j, les 6 orbites libres (une par paire jr) forment **un seul 6-cycle** de R_j. Les cycles libres sont pairs
  et de longueur ≥ 4, et un 2-cycle serait une paire D passant par j. Les voisins de x qui rencontrent j sont ses deux
  R_j-voisins sur ce cycle ; on note n_r = 1 si l'un est dans la paire jr.
- **(F3)** μ = 2 avec chaque élément non contenu dans x : x a exactement 2 voisins contenant chaque élément d'un bloc
  t ∉ {i,j}, donc exactement 4 étiquettes voisines rencontrent t.
- **(F4)** Donc Σ_{P ⊂ R_x} φ_x(P) = 6 + n_k (c'est ω(x) = n_k de 18.3), et φ ≤ 2 sur ces paires (16.1). D'où
  **#doubles − #zéros = n_k ≥ 0** parmi les 6 paires de R_x.
- **(F5) Comptages exacts d'un triangle D complet τ (C2 §25, preuve par positivité, revérifiée ici) :**
  |N(x) ∩ Lab(τ)| = 2 / 1 / 3 selon que la paire de x est dans τ / rencontre τ en un bloc / est disjointe de τ.
- **(B) Budget (16.1, revérifié ; n_s = 0).** Σ_{84 étiquettes} doubles = 84 au plus, soit **Σ_{42 orbites} doubles ≤ 42**.
- **Zéro exporté (18.3).** Si un R_j-voisin de x est dans la paire jr d'un triangle τ = {j,r,s} avec r, s ∈ R_x, alors
  τ rencontre la paire de x en j seulement. Par (F5) le compte vaut 1, déjà atteint, donc φ_x(rs) = 0.
- **Contraintes d'ordre** (toujours par F5) :
  - (a) si le cycle de x est le triangle {i,j,k}, x n'a pas de R_j-voisin dans jk ;
  - (b) pour un triangle {j,r,s} ∌ i, au plus un des deux R_j-voisins de x tombe dans jr ou js.

## 2. Théorème S8.1 — le type (3,3,3,3,3,6) est impossible
**Structure.** Soit 0 le bloc absent du 6-cycle. Il porte 3 triangles, donc les 6 autres blocs sont appariés en
antipodes a/a′, b/b′, c/c′. Les 12 arêtes restantes forment l'octaèdre sur {a,…,c′}. Deux faces opposées F = {a,b,c}
et F′ = {a′,b′,c′} sont les deux derniers triangles, et l'équateur a–b′–c–a′–b–c′–a est le 6-cycle. Le type est unique
à isomorphisme près (420 = 5040/12).
**Bloc 0.**
- Les 6 orbites libres en 0 sont indexées par a,…,c′. Pour celle de a : k = a′ (triangle propre), R_x = {b, b′, c, c′}.
- Par (a), pas de voisin a′. Par (b), pas deux voisins dans {b,b′} ni deux dans {c,c′}. Donc un voisin dans chacune de
  ces deux paires, d'où 2 zéros ({b,b′} et {c,c′}).
- Idem pour les 6 orbites : **12 zéros**.
**Bloc j = a (idem pour les 5 autres, par symétrie).**
- Orbites libres indexées par 0, a′ (triangle {0,a,a′}), b, c (face F) et b′, c′ (équateur, k_{b′} = c, k_{c′} = b).
- On pose S₁ = {0, a′}, S₂ = {b, c}, H = {b′, c′}. Arêtes interdites : 0–a′ et b–c (a). Aucun sommet n'a ses deux
  voisins dans S₁, ni ses deux voisins dans S₂ (b).
- Zéros :
  - une arête S₁–S₂ donne un zéro à chaque bout ({b,c} côté S₁, {0,a′} côté S₂) ;
  - une arête S₁–H donne un zéro au bout H ({0,a′}) ;
  - les autres arêtes n'en donnent aucun (b, c ∉ R_{b′} ∪ R_{c′} à la fois).
  Donc #zéros = e(S₁,H) + 2e(S₁,S₂) = 4 + e(S₁,S₂), puisque S₁ n'a pas d'arête interne.
- **e(S₁,S₂) ≥ 2.** Si e = 0 ou 1, S₁ occupe au moins 3 des 4 places de H, et les 4 extrémités d'arêtes de S₂ (sans
  arête interne) n'ont plus où aller. Et e = 2 est réalisé (0–b′–b–a′–c′–c–0).
- **Au moins 6 zéros par bloc.**
**Conclusion.** Σ_{orbites} doubles ≥ Σ zéros ≥ 12 + 6·6 = **48 > 42**. Contradiction avec (B). ∎
*Contrôle machine : `scripts/zeros_budget.py` donne 48 exactement, sur les 420 décompositions étiquetées.
`scripts/doubles_budget3.py` (système local exact) donne 54.*

## 3. Système local complet et résultats (calcul exhaustif, contrôlé)
**Système local par orbite** (x ∈ X, rôles i = 0, j = 1, k = 2, R_x = {3,4,5,6}). Inconnues : φ(P) ∈ {0,1,2} pour les 6
paires P ⊂ R_x ; y_r = φ_x(kr) ∈ {0,1,2} ; n_r = 1 pour les deux paires jr des voisins de x sur le 6-cycle libre en j.
| contrainte | origine |
|---|---|
| C0 : n_k = 0 si le cycle propre est le triangle {i,j,k} | F5 (compte 2 dans son propre triangle, atteint par v, σv) |
| C1 : φ(ab) = 1 pour un triangle {i,a,b} | F5 (compte 1) + F1 |
| C2 : n_a + n_b + φ(ab) = 1 pour un triangle {j,a,b} | F5 + F2 |
| C3 : n_k + n_b + y_b = 1 pour un triangle {j,k,b} | F5 + F2 |
| C4 : y_c + y_d + φ(cd) = 3 pour un triangle {k,c,d} | F5 (compte 3) |
| C5 : Σφ = 3 sur un triangle contenu dans R_x | F5 |
| C6 : Σ_r y_r = 2 − n_k | μ = 2 avec (k,0), (k,1) ; v, σv comptent une fois chacun |
| C8 : n_t + y_t + Σ_s φ(ts) = 4 pour t ∈ R_x | μ = 2 avec (t,0), (t,1) ; φ(it) = 0 par F1 |
doubles(x) ≥ #{φ = 2} + #{y = 2}. Budget : Σ_{orbites} doubles ≤ 42.

**Lemme S8.2 (orbite morte).** Supposons que le cycle propre de x soit le triangle {i,j,k} et qu'il existe deux triangles
{j,a,b} et {k,c,d} avec {a,b,c,d} = R_x. Alors le système est impossible.
*Preuve.* On somme C8 sur {a,b} et sur {c,d}, avec X = Σ φ croisés :
  N_ab + (y_a + y_b) + 2φ(ab) + X = 8   et   N_cd + (y_c + y_d) + 2φ(cd) + X = 8.
On y substitue :
- N_cd = 2 − N_ab (C0 : les deux voisins sont dans R_x) ;
- φ(ab) = 1 − N_ab (C2) ;
- y_c + y_d = 3 − φ(cd) (C4) ;
- y_a + y_b = φ(cd) − 1 (C6 avec n_k = 0).
La différence des deux sommes vaut alors **−4 = 0** : contradiction. ∎
*Conséquence.* Toute décomposition contenant une chaîne de trois triangles T₁ – T – T₃ (T₁ et T₃ passent par deux
sommets distincts de T et sont disjoints hors de T) est impossible : l'orbite de la paire {i,j} appariée par i meurt, où
i est le troisième sommet de T.

**Lemme S8.3 (Fano, second mécanisme).** En un bloc porté par 3 triangles de Fano, chaque orbite libre n'admet que ses deux
« paires C1 » comme voisins. Par exemple, pour {3,6} au lieu de {3,4} ou {5,6}, C4 et C8 forcent φ(45) = 3. Suivre une
paire C1 ferme alors un triangle (0–3–4 ou 0–5–6) dans le 6-cycle : aucun ordre possible.

**Résultats** (`scripts/doubles_budget3.py`, toutes les décompositions étiquetées ; `logs/doubles_budget3.log`) :
| type | étiquetées | impossibles d'emblée | borne des survivantes | verdict |
|---|---|---|---|---|
| (3⁷) | 30 | 30 | – | exclu |
| (3⁵,6) | 420 | 420 | – | exclu (et §2) |
| (3⁴,4,5) | 5040 | 5040 | – | exclu |
| (3³,4³) | 5565 | 3780 | 48–51 > 42 | exclu |
| (3³,5,7) | 2520 | 2520 | – | exclu |
| (3³,6,6) | 1890 | 1890 | – | exclu |
| (3²,4,4,7) | 2520 | 0 | 28 | survit |
| (3²,4,5,6) | 25200 | 5040 | 26–34 | 20160 survivent |
| (3²,5³) | 5880 | 840 | 30 | 5040 survivent |
| types à 1 triangle | – | 0 | 9–16 | survivent |
| types sans triangle | – | – | système sans contrainte de triangle | hors de portée |

**Contrôles.**
- Tests unitaires : Fano, clé (3,4) → 2 ; (3,5) et (2,3) → impossibles, refaits à la main.
- Solveur indépendant au niveau des étiquettes (`labellevel_check.py`), comparé sur **toutes les 670 clés locales**
  (`crosscheck_doubles.py`) : 353 faisables avec le même minimum, 317 impossibles des deux côtés, **0 désaccord**
  (`logs/crosscheck_t3.log`, `logs/crosscheck_extra.log`).
- Relabellisation canonique contre relabellisation fixe : 0 écart sur 2 × 300 décompositions au hasard.
- Explications bloc par bloc : `logs/explain_blocks.log`.


## 3 bis. Part prouvée à la main contre part calculée (`scripts/chain_free.py`)
- **Chaînes.** Dès que deux triangles sont disjoints, tout troisième triangle crée une chaîne. Il doit contenir le
  7ᵉ point, et le reste de ses sommets est pris dans les deux triangles.
- **Tuées par le lemme S8.2 (à la main)** : (3⁵,6) 420/420 ; (3³,5,7) 2520/2520 ; (3⁴,4,5) 2520/5040 ;
  (3³,4³) 3150/5565 ; (3³,6,6) 1260/1890.
- **Sans chaîne (triangles deux à deux sécants), tuées par le calcul seulement** :
  - Fano 30 (aussi le Théorème F de C2, et S8.3) ;
  - quadrilatère (3⁴,4,5) 2520, d'emblée ;
  - étoile (3³,6,6) 630, d'emblée ;
  - étoile (3³,4³) 630 d'emblée plus 105 par le budget (48) ;
  - triangle de triangles (3³,4³) 1680, par le budget (51).
- Relecture de C2 (25/09) : C0–C8 et S8.2 re-dérivés à la main, « TIENT », inscrits à son carnet §26. Le calcul
  exhaustif y est mentionné comme calcul du labo « structures », recoupé par un second solveur.

## 4. Erreurs rencontrées en chemin (consignées)
- Premier contrôle Fano à 42, puis 70 au lieu de 84 : il manquait les contraintes (a) et (b). Le contrôle a joué son
  rôle.
- Le premier `doubles_budget2` utilisait « Σφ ≥ 6 » au lieu de « = 6 + n_k », et oubliait les triangles {k,a,b}
  (compte 3). C'était plus faible, mais valide ; corrigé dans doubles_budget3.
- Script de recomptage : formule `ordered*3//6` au lieu de `ordered//6` (2880 au lieu de 960). L'erreur était dans
  le contrôle, pas dans la classification ; corrigée.
