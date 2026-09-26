# S17 — Cœur tout-tordu (1 triangle et sans triangle) : calcul certifié par (P2)

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

### 2.2 Sans triangle
(en cours)

## 3. Portée et réserves
- La preuve est **assistée par ordinateur**. Chaque certificat est un objet fini, vérifiable en entiers ; il ne dépend
  pas du solveur qui l'a trouvé.
- Deux dépendances sont partagées entre le générateur et le vérificateur :
  - le solveur local exact `solve_core`, lui-même recoupé avec doubles_budget3 sur 1340 clés (S11) ;
  - la liste des décompositions de K₇ (`zeros_budget.decompositions`, dont les effectifs par type ont été contrôlés,
    121 215 au total).
