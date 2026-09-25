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

**Lemme S8.4 (du complément, CL).** Si k porte un triangle {k,c,d} avec c, d ∈ R_x, alors, pour {a,b} = R_x ∖ {c,d},
**φ(ab) = 3 − N_ab** (N_ab = n_a + n_b). En particulier N_ab ≥ 1, et N_ab = 1 force le double φ(ab) = 2.
*Preuve.* Demi-sommes de C8 sur {a,b} et sur {c,d}, avec X = Σ φ croisés :
  N_ab + (y_a + y_b) + 2φ(ab) + X = 8   et   N_cd + (y_c + y_d) + 2φ(cd) + X = 8.
On y substitue :
- y_c + y_d = 3 − φ(cd) (C4) ;
- y_a + y_b = φ(cd) − 1 − n_k (C6) ;
- N_ab + N_cd = 2 − n_k.
La différence vaut 2N_ab − 6 + 2φ(ab) = 0. ∎
- S8.2 en est un cas particulier : C2 sur {j,a,b} donne φ(ab) = 1 − N_ab, contre 3 − N_ab. **Le triangle propre n'y sert
  pas.**
- ML1 (§5) est CL appliqué deux fois.
- Relu par C2 (25/09) et inscrit comme lemme général à son §26 (valable en Σ = ∅, carrés « un bloc », triangles D
  complets) ; C2 compte l'employer pour les mélanges.

**Corollaire S8.5 (deux triangles disjoints).** Une décomposition tout-tordue qui contient deux triangles disjoints est
impossible.
*Preuve.* Soient T, T′ les deux triangles et i le 7ᵉ bloc. Les trois cycles par i apparient ses six voisins deux à deux.
Comme |T| = 3 est impair, l'un de ces cycles passe par s – i – t avec s ∈ T et t ∈ T′ (s, t sont des blocs). Au bloc
j = s, l'orbite de la paire {i,s} appariée par i a k = t et R_x = (T ∖ s) ∪ (T′ ∖ t). Pour {a,b} = T ∖ s, C2 sur T
donne φ(ab) = 1 − N_ab et CL sur T′ donne φ(ab) = 3 − N_ab : l'orbite est morte. ∎
- Avec au moins 3 triangles, c'est déjà le cas des chaînes (§3 bis).
- **Nouveau : les types à exactement deux triangles disjoints**, soit (3²,4,5,6) 5040 et (3²,5³) 840. Ce sont
  **exactement** les « impossibles d'emblée » du calcul exhaustif. Aucune décomposition (3²,4,4,7) n'a deux triangles
  disjoints.
- Contrôle (`scripts/two_disjoint_triangles.py`, `logs/two_disjoint_triangles.log`) :
  - la clé morte est infaisable pour les 10 paires de voisins, avec ou sans triangle propre ;
  - l'orbite morte est trouvée et vérifiée dans chacune des 5880 décompositions ;
  - les effectifs coïncident avec ceux du calcul exhaustif.
- Il reste, à deux triangles, les 27 720 décompositions à triangles **sécants**. Leurs bornes locales exactes valent
  26 à 34, en dessous de 42 : le système local seul ne suffira pas.

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
| (3²,4,5,6) | 25200 | 5040 (= triangles disjoints, S8.5) | 26–34 | 20160 survivent |
| (3²,5³) | 5880 | 840 (= triangles disjoints, S8.5) | 30 | 5040 survivent |
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
  - étoile + transversale (3⁴,4,5) 2520, d'emblée. *(Correction du 25/09, relevée par C2 : j'avais écrit
    « quadrilatère, 4 triangles se coupant en 6 points distincts ». C'est faux : les intersections sont 2, 2, 2,
    1, 5, 6. Le vrai quadrilatère complet est impossible en (3⁴,4,5), car le 7ᵉ point n'aurait que le degré 4.
    Le calcul portait bien sur la bonne famille.)*
  - étoile (3³,6,6) 630, d'emblée ;
  - étoile (3³,4³) 630 d'emblée plus 105 par le budget (48) ;
  - triangle de triangles (3³,4³) 1680, par le budget (51).
- Relecture de C2 (25/09) : C0–C8 et S8.2 re-dérivés à la main, « TIENT », inscrits à son carnet §26. Le calcul
  exhaustif y est mentionné comme calcul du labo « structures », recoupé par un second solveur.
- **Mise à jour (25/09, soir).** Les familles sans chaîne ont toutes une preuve à la main, relue par C2 (§5). **Tous les
  types tout-tordus à au moins 3 triangles sont donc exclus entièrement à la main.** Les types à deux triangles
  disjoints aussi (S8.5).

## 4. Erreurs rencontrées en chemin (consignées)
- Premier contrôle Fano à 42, puis 70 au lieu de 84 : il manquait les contraintes (a) et (b). Le contrôle a joué son
  rôle.
- Le premier `doubles_budget2` utilisait « Σφ ≥ 6 » au lieu de « = 6 + n_k », et oubliait les triangles {k,a,b}
  (compte 3). C'était plus faible, mais valide ; corrigé dans doubles_budget3.
- Script de recomptage : formule `ordered*3//6` au lieu de `ordered//6` (2880 au lieu de 960). L'erreur était dans
  le contrôle, pas dans la classification ; corrigée.

## 5. Versions à la main des familles sans chaîne (pour la consigne « sans solveur » de C2) — 25/09 soir
Relues et acceptées par C2 au fur et à mesure (§26 de son carnet) : Fano, étoile (6,6), étoile (4,4,4) classe B,
étoile + transversale, ainsi que les classifications à la main (octaèdre : types I et II des 6-cycles, équateurs et
cerfs-volants, impossibilité de trois cerfs-volants ; étoile + transversale : 4-cycle forcé).
- **ML1** : si k porte deux triangles qui couvrent R_x, alors n_k = 0, il y a exactement un voisin dans chaque paire, et
  φ(ab) = φ(cd) = 2.
- **Lemme D** : #doubles − #zéros = n_k parmi les paires de R_x.
- **Fano** : les paires admissibles sont exactement les paires C1 ; elles ferment un triangle dans le 6-cycle.
- **Étoile (6,6), bloc 0** : faits F_a à F_e, puis six cas sur les voisins de l'orbite 2.
- **Étoile B, bloc 6** : ML1 pour l'orbite 4, exclusions pour les orbites 0, 1, 3, 5 ; fermeture en deux cas.
- **Étoile + transversale, bloc 1** : ML1 et C2 pour l'orbite 0 ; exclusions pour 2, 3, 4, 5, 6 ; trois cas.
- **Étoile A (trois équateurs)** : budget.
  - Centre : au moins 2 doubles par orbite (C2 et lemme D), soit 12.
  - Chaque bloc non central j : l'orbite de l'antipode de j (ML1) a 2 doubles, et ses deux voisins, qui ont n_k = 1,
    en ont au moins 2 chacun, soit au moins 6.
  - Total au moins 48 > 42.
- **Triangle de triangles (famille 4) : au moins 43 > 42.** Relu en entier par C2, énumérations comprises (« TIENT »,
  25/09), inscrit à son §26.
  - Cadre : T1 = {0,1,2}, T2 = {1,5,6}, T3 = {2,4,6}, et w = 3 sur aucun triangle. Z₃ agit par ρ = (1 6 2)(0 5 4).
    Le budget exact vaut 51 ; il suffit d'atteindre 43.
  - Outil : CL (S8.4), avec les exclusions immédiates C0/C2/C3 et la borne immédiate max(#φ forcés à 2,
    n_k + #zéros exportés) du lemme D.
  - **Bloc w = 3 : au moins 4.**
    - Règles CL : 0 doit avoir un voisin dans {1,5}, 1 dans {0,5}, 2 dans {0,4}, 4 dans {0,2}, 5 dans {4,6},
      6 dans {4,5}. Une orbite qui n'a qu'un voisin dans sa paire gagne un double.
    - On note z le nombre d'orbites qui ont leurs deux voisins dans leur paire. On a z ≤ 3, sinon on obtient un
      triangle 0-1-5, 0-2-4 ou 4-5-6.
    - z = 3 force le cycle (0,1,5,6,4,2). Le fait X1 (orbite 5, paire {1,6}, au moins 2) donne alors 1 + 2 + 1.
  - **Bloc d'intersection 1 : au moins 7.**
    - C0, C2 et C3 placent 5 et 6 à distance 3, ainsi que 0 et 2.
    - Avec les règles CL 0{3,5} et 5{0,3}, il reste 3 cycles, chacun avec les bornes 1+2+1+1+1+1.
  - **Bloc privé 0 : au moins 6.**
    - 1 et 2 sont à distance 3. Avec les quatre règles CL, il reste 5 cycles, chacun à 6.
    - Deux faits servent : X2 (orbite 3, paire {1,4}) et X3 (orbite 5, paire {1,6}), chacun au moins 2.
  - Total : 4 + 3·7 + 3·6 = **43**. ∎
  - *X3 sans énumération.* On a y1 = 0, φ14 = 2, φ24 = 1 + y6 et φ26 = 3 − y2 − y6 (C8(2)). Supposons qu'il n'y ait
    pas de second double. Alors φ24 ≤ 1, donc y6 = 0 ; puis φ26 ≤ 1 force y2 = 2, ce qui est un double.
    Contradiction.
  - Contrôle (`scripts/verify_pasch_hand.py`, avec son propre solveur `scripts/local_system.py` en noms d'origine ;
    `logs/verify_pasch_hand.log`) :
    - CL tient sur les 420 clés du triangle de triangles et sur les 1340 clés canoniques des calculs exhaustifs. Le
      solveur propre concorde avec doubles_budget3 sur toutes ces clés.
    - X1, X2 et X3 tiennent pas à pas, et leur minimum exact vaut 2.
    - Les règles laissent exactement 8, 3 et 5 cycles. Les bornes données sont celles des règles et ne dépassent pas
      les minima exacts. Tout cycle exclu est réellement infaisable.
    - ρ est un automorphisme. Les minima exacts par bloc sont 6 pour w, 8 pour les intersections et 7 pour les privés.
Contrôles machine de chaque ensemble de faits : `scripts/check_star66_facts.py`, `scripts/check_facts_generic.py`,
`scripts/budget_block_detail.py`. Pour le triangle de triangles, la recherche du plan le moins coûteux se trouve dans
`scripts/pasch_greedy2.py` et `scripts/pasch_rules.py`. `scripts/pasch_cheapest.py` n'a servi que de recherche exacte, abandonnée
au-delà de 5 minutes. `scripts/show_system.py` affiche le système local d'une orbite.
