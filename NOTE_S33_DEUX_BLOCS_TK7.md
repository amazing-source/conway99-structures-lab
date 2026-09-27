# NOTE S33 : deux blocs couplés pour T = K₇ (le tout-tordu non simple)

27/09, ~8 h. « structures ». Code : `involution/sat_twoblock.py`, liste `involution/twoblock_instances.json`.
DÉMONTRÉ = preuve écrite ci-dessous ; VÉRIFIÉ = contrôle machine sur un vrai graphe ou un monde synthétique.

## 0. Pourquoi
- Le tout-tordu pur (T = K₇ : les 21 paires tordues, les 7 blocs de profil (0,6)) n'est exclu que pour les circuits D SIMPLES
  (erratum S31).
- Le cas non simple est couvert si les systèmes locaux (0,6) sont infaisables. Mais un bloc seul résiste :
  - le PL de l'étape 2 est faisable sans SEC ;
  - le SAT mono-bloc avec `--sec --twist` est inconnu à 240 s (glucose, portable).
- Or en T = K₇ **deux blocs quelconques c, d sont tous deux (0,6)**, et leurs systèmes voient des sommets communs. On les
  impose donc ensemble.

## 1. Repère monde et descriptions (DÉMONTRÉ)
- Sommets de N(x) : 2k (+) et 2k+1 (−) du bloc monde k, avec c = bloc 0, d = bloc 1 et σ(v) = v ⊕ 1.
- Chaque sommet extérieur est une paire de sommets de N(x) de blocs différents. On distingue :
  - 20 c-racines non partagées, 20 d-racines non partagées ;
  - 4 partagées (c±, d±) ;
  - 40 « lointains-lointains », qui ne touchent ni c ni d.
- **Le système en c décrit adj(X, Y) dès que X ou Y est une c-racine** :
  - c-racine × c-racine : par le type (μ_c, h_c), avec (c+,v) ~ (c+,μv), (c+,v) ~ (c−,hv) et les images par σ ;
  - c-racine × sommet lointain en c : par U_c, avec U[z,b] = [z ~ (c+,b)] et (c−,u) ~ Y ⟺ U[σY, πu].
- Idem en d, à travers le recollement G : étiquettes du repère d → N(x) ∖ d, avec G(γ) = c+.

## 2. Couplage : toute paire décrite deux fois (DÉMONTRÉ)
Chaque paire décrite par les deux systèmes donne une équivalence entre les deux descriptions. Modulo σ :

| paires | nombre | contenu |
|---|---|---|
| c-racine × d-racine | 200 | U_c[(d+,w), v] = U_d[(c+,v), w] et U_d[(c−,u), w] = U_c[(d−,πw), πu] (C1/C2) |
| c-racine × partagée | 40 | unités sur U_d données par le type de c : [μ_c(d+) = v], [h_c(v) = d+], [h_c(d+) = u], [μ_c(d−) = πu] |
| partagée × d-racine | 40 | unités sur U_c données par le type de d |
| **partagée × lointain-lointain** | **80** | **U_c[Y, d±] = U_d[Y, c±]** : (c±, d±) est racine des deux côtés. Absent de la spec initiale |
| partagée × partagée | 4 | compatibilité des types ; automatique pour les 7 types (0,6), qui n'ont ni h fixe ni μ = π |

- Toutes les autres paires ne sont décrites que par un système (c × c, c × lointain, d × d, d × lointain), ou par aucun
  (lointain × lointain).
- **Cohérence des pivots de la paire cd** : [d+ pivot en c] = [c− pivot en d].
  - L'orbite o_{d+} en c est l'orbite de l'étiquette c+ en d.
  - Une orbite qui pivoterait en c ET en d aurait son partenaire D dans la paire cd, donc serait frère. C'est exclu en (0,6).

## 3. CROSS : le (iii) à deux blocs (DÉMONTRÉ)
Soient r une c-racine et s une d-racine, non partagées. Un voisin commun Z de r et s est :
- un sommet de N(x) (celui que r et s partagent, s'il existe) ;
- ou une c-racine : adj(r,Z) par le type c, adj(s,Z) par U_d, ou par le type d si Z est partagée ;
- ou une d-racine : symétrique ;
- ou un lointain-lointain : U_c ∧ U_d.

L'adjacence lointain × lointain, inconnue, n'intervient jamais, puisque r et s sont des racines. Le compte est donc EXACT :
  **[r, s partagent un sommet de N(x)] + Σ_Z adj(r,Z)·adj(s,Z) + adj(r,s) = 2**   (λ = 1, μ = 2).

Cela fait 200 équations modulo σ. Pour les paires (racine, lointain-lointain), on n'a qu'une inégalité : non codée.

## 4. Réduction des instances (DÉMONTRÉ)
- On choisit c = un bloc de type MINIMAL (indice de classe parmi les 7) et d = son bloc 0 dans l'étiquetage canonique. Le
  label 0 est pivot pour les 7 représentants, donc γ est non-pivot en d.
- Il reste t_c ≤ t_d, soit 28 paires. γ est pris modulo Stab(t_d) :
  - les stabilisateurs sont d'ordre 6 (types 0, 3, 6), avec 1 orbite de γ ;
  - ou d'ordre 2 (types 1, 2, 4, 5), avec 3 orbites.
- Le recollement des 5 autres blocs est une VARIABLE SAT : 50 variables B[β,k,o], bijection signée.
- Le stabilisateur de γ est trivial partout : aucune symétrie résiduelle. Au total, **60 instances**.

## 5. Contrôles (VÉRIFIÉ, `python sat_twoblock.py --selftest`)
- **BvLS, graphe réel m = 11**, 3 repères aléatoires (c, d, p, γ et recollement tirés au hasard) :
  - le couplage (4 + 72 + 72 + 288 + 648 paires) et les 648 équations CROSS tiennent EXACTEMENT sur la vraie structure ;
  - contrôles négatifs : un mauvais recollement donne 441 violations, un seul U_c retourné en donne 4.
  - ⚠ BvLS est dégénéré (μ = τ = π) : les unités y sont triviales.
- **Monde synthétique m = 7, avec les vrais types (0,6)** : types exacts par définition, tout le reste aléatoire et
  σ-invariant.
  - Le couplage tient avec le vrai recollement fixé.
  - Avec le recollement en variables, le solveur trouve le vrai recollement, et lui seul.
  - Un U retourné est détecté.
- **Positif SAT** : (i) + ponctuelles aux deux blocs + couplage, sans (ii) ni (iii), est SAT en 0,5 s.

## 6. Calibrage (portable, glucose, un seul processus, ≤ 240 s) sur l'instance 0 (types (0,0), γ = 1, `--sec --twist`)

| variante | statut |
|---|---|
| recollement « identité » (les deux types alignés sur les 5 autres blocs), avec CROSS | **UNSAT en 20,5 s** |
| même recollement, sans CROSS | UNSAT en 226 s |
| même recollement, sans SEC ni twist | inconnu à 240 s |
| mono-bloc t06_0 seul, SEC + twist | inconnu à 240 s |
| 3 recollements aléatoires fixés | inconnus à 80 s |
| 8 recollements aléatoires sous hypothèses | inconnus à 30 s |
| instance complète (recollement en variables) | inconnue à 240 s |

**Lecture.**
- Le couplage mord : un recollement meurt alors que le mono-bloc résiste. CROSS accélère d'un facteur 10 ; SEC et twist
  sont nécessaires à cette vitesse.
- Ce n'est pas une mort instantanée : le recollement identité est un cas facile, pas typique.
- Suite : kissat sur le PC du labo (sonde identité pour calibrer, une sonde aléatoire, puis l'instance 0 complète). Si c'est
  encore trop faible : trois blocs (triangle c, d, e), ou le modèle complet avec les variables lointain × lointain, où tous
  les comptes λ/μ deviennent exacts.
