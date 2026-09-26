# NOTE S28 : les configurations AD meurent au niveau « paires »

26/09, soir. Demande de C2 (note 41 §4) : montrer que toute configuration AD meurt. AD signifie Σ = ∅ et que chaque paire
est soit ancrée, soit doublement disjointe (dd). Code dans `ad/`.

## 1. Outil : l'équilibre par composante, forme générale
**Rappels.**
- En un bloc k, les composantes de G_k sont les orbites de ⟨σ, μ_k, τ_k⟩ sur les étiquettes. Ce sont des réunions de blocs.
- En AD, aucun bloc n'a de paire à moitié. Donc (39.4) M_k (les paires ancrées en k) et F_k sont des réunions de composantes.
  Les composantes contenues dans M_k sont exactement les **D-cycles au pivot k** : digone si d_k = 2, triangle si d_k = 3,
  et l'un des 4 types si d_k = 6.

**Équilibre (note 23 §8, démontré).** Soit S une réunion de composantes de G_k, vue comme ensemble de blocs, et w une
racine qui évite le bloc k. Alors

  #(intérieurs de w dans S) + Σ_{m ∈ S} φ_w(km) = |S|,

où φ_w(P) désigne le nombre de voisins de w parmi les 4 racines de la paire P.
- Pour un digone, cela redonne 16.2 (2/0/1/2).
- Pour un triangle, cela redonne le lemme du triangle de C2 (0/2/1/3), que j'ai relu : f = 3e_x + 3·1_{(c,±)} + 1_I + 1_Lab,
  ‖E₋₄f‖² = (24 − 2E_L)/7.

**Fermeture des classes ancrées** (le point qui rend l'équilibre applicable ; relecture de C1). Soit m ∈ In(k), avec l'orbite
{k±, m±}.
- Ses deux voisins étiquetés k sont les deux racines de son partenaire D, qui est dans une paire km′ avec m′ ∈ In(k).
- Donc τ_k et μ_k (m_k et π_k) envoient les étiquettes de In(k) dans In(k). σ les préserve aussi, car la paire km est ancrée :
  ses deux orbites sont appariées via k.
- In(k) est donc stable. Les D-cycles en k en sont les composantes, et F_k, son complément, est une réunion de composantes.

**Zéros par bloc (démontrés).** Soit W une orbite qui passe par le bloc t.
- **Si W est appariée via t**, ses deux voisins par t sont les deux racines de D(W) ⊂ tm′. Donc φ_W(tm′) = 2, et 0 sur les
  autres paires par t.
- **Sinon**, ses voisins par t sont dans des orbites libres en t. Donc φ_W(tm) = 0 pour m ∈ In(t).

## 2. Lemme T⁺ (« queue chargée », démontré)
**Énoncé.**
- **Hypothèses.**
  - Soit X ⊂ ij une orbite appariée par le pivot i, avec D(X) ⊂ ik.
  - Soit S la composante de i dans G_k.
  - Soit Π l'ensemble des blocs p ∈ S ∖ {i, j} en lesquels {i, k} est une composante de G_p. En AD, cela veut dire que
    p porte le digone {i,k}.
- **Conclusion.** |S| ≥ 3 + [j ∈ S] + |Π|.

**Preuve.** Soit x ∈ X, et p ∈ Π.
- x évite p. L'équilibre en p sur {i,k} donne φ_x(pi) + φ_x(pk) = 2 − 1 = 1, car l'intérieur de x au bloc i est dans {i,k}.
- Or φ_x(pi) = 0 : les voisins de x par i sont les deux racines de D(X) ⊂ ik, et p ≠ k. Donc φ_x(pk) = 1.
- L'équilibre en k sur S compte alors, parmi ses |S| voisins, au moins :
  - l'intérieur de x en i ;
  - les 2 racines de D(X) ;
  - l'intérieur de x en j si j ∈ S ;
  - une racine de pk pour chaque p ∈ Π.
  Ces sommets sont deux à deux distincts. ∎

La queue 25.2 de C2 est le cas Π = ∅.

## 3. AD2 : les 14 survivants de C2 meurent, et 12 d'entre eux par T⁺ en une ligne
Rappel du cadre : les paires dd sont 01 et 23, appariées en croix, et le bloc A = 6 a d₆ = 6. **Aucune de ces preuves n'utilise
la structure D du bloc A.** Elles valent donc pour ses 130 options. Notation : X ⊂ ij pivote en i, avec D(X) ⊂ ik.

| rep | i, j, k | S = F_k | j ∈ S | Π (digones {i,k}) | 3 + [j∈S] + \|Π\| > \|S\| = 4 |
|---|---|---|---|---|---|
| 69 | 1, 4, 5 | {1,3,4,6} | oui | {3} | 5 |
| 130 | 0, 2, 3 | {0,2,5,6} | oui | {5} | 5 |
| 132 | 0, 2, 4 | {0,2,5,6} | oui | {5} | 5 |
| 183 | 0, 3, 2 | {0,3,5,6} | oui | {5} | 5 |
| 191 | 5, 0, 4 | {0,1,5,6} | oui | {1} | 5 |
| 214 | 4, 0, 5 | {2,3,4,6} | non | {2,3} | 5 |
| 223 | 0, 2, 3 | {0,2,5,6} | oui | {5} | 5 |
| 235 | 0, 2, 3 | {0,2,4,6} | oui | {4} | 5 |
| 307 | 3, 0, 5 | {0,3,4,6} | oui | {4} | 5 |
| 313 | 3, 0, 4 | {0,1,3,6} | oui | {1} | 5 |
| 314 | 0, 2, 5 | {0,2,3,6} | oui | {3} | 5 |
| 320 | 5, 0, 1 | {0,2,5,6} | oui | {2} | 5 |

Dans chaque ligne, S = F_k est une seule composante. En AD2, une composante de F_k à 2 blocs exigerait deux paires dd par k
(cas S de 11.6). Il n'y a donc que des composantes à 3 blocs ou plus.

**Reps 251 et 263** : une racine, avec 8 équations d'équilibre, et une identité linéaire contradictoire.
- **Rep 251** : w ∈ 05, pivot 0, D(w) ⊂ 02.
  - Zéros : φ(0·) = 0 hors 02 ; φ(35) = φ(45) = 0.
  - Équilibres :
    - en 1 : M = {3,4,5}, φ(13)+φ(14)+φ(15) = 2, et F = {0,2,6}, φ(12)+φ(16) = 2 ;
    - en 2 : M = {1,4,5}, φ(12)+φ(24)+φ(25) = 2, et F = {0,3,6}, φ(23) = φ(26) = 0 ;
    - en 3 : F = {1,2,5,6}, φ(13)+φ(36) = 3 ;
    - en 4 (d = 0) : φ(14)+φ(24)+φ(46) = 3.
  - F1 en 6 : φ(16)+φ(36)+φ(46)+φ(56) = 4. F1 en 5 : φ(15)+φ(25)+φ(56) = 2.
  - En éliminant, on obtient 4 = 6.
- **Rep 263** : même schéma, et l'on obtient φ(13)+φ(14)+φ(15) = 4 contre 2.

**Contrôle.** Les 14 sont retrouvés par `ad/proofs_ad2.py`, qui extrait un sous-système minimal infaisable.

**Relecture indépendante (C1, 26/09 soir) : les 14 preuves sont correctes.** C1 a tout refait à la main, arithmétique et
applicabilité aux données In et D de chaque rep. Pour 251, c'est un certificat linéaire pur :
- posons A = φ13, B = φ14, C = φ15, D = φ12, E = φ24, F = φ25, G = φ16, H = φ36, I = φ46, J = φ56 ;
- on a (1) C+F+J = 2, (2) A+B+C = 2, (3) D+G = 2, (4) D+E+F = 2, (5) A+H = 3, (6) B+E+I = 3, (7) G+H+I+J = 4 ;
- la combinaison (1) − (2) + (3) − (4) + (5) + (6) − (7) donne 0 = 2.

263 se traite de même. Les 14 contradictions sont toutes infaisables au niveau PL.

## 4. Toutes les classes AD, n_dd ≥ 2 : balayage au niveau paires
**Modèle** (`ad/pair_level_relaxed.py`). Il ne contient que des faits démontrés :
- F1 et les zéros par bloc ;
- l'équilibre sur les D-cycles connus et sur F_k ;
- les partenaires dis, φ ∈ {2,3} ;
- la symétrie e(P,P′) ;
- la **dichotomie par paire**. Elle est identique au budget μ 16.1 et vient de (Q²)_{P⁺P⁻} = 2 − n_s :
  #{Z hors partenaires D : φ_Z(P) = 2} + 2·#{Z : φ_Z(P) = 3} = 2 − n_s(P).

Les structures libres (bloc A, blocs de degré 4, graphe D des paires dd) sont en variables de partenaire. L'équilibre n'y est
imposé que sur les réunions certaines. C'est donc une relaxation valide, et une classe infaisable est morte pour tous les
choix.

**Énumération.** C'est celle de C2 (`ad_enum_all.py`), avec les lemmes de degré et de queue. Je l'ai relue.
- Le test « admits » est une permutation sans point fixe qui envoie chaque paire dd sur une paire dd disjointe. C'est exactement
  l'existence du graphe dis.

**Résultats.**

| n_dd | classes | mortes | remarque |
|---|---|---|---|
| 2 | 14 | 14 | toutes à une racine (§3) |
| 3 | 8 | 8 | une seule demande la dichotomie (rep 20 : 3-cycle dd, A en 6-cycle) |
| 4 | 510 | 510 | 2 demandent la dichotomie |
| 5 | 853 | 853 | |
| 6 | 4 594 | 4 594 | |
| 7 à 21 | en cours | | 7 : quelques vivantes dans la version relâchée, à retester en énumérant les structures |

## 5. Hors AD : feuilles de mélanges de C2 (Σ = ∅, tout « un bloc »)
**Modèle général** (`ad/general_model.py`) : squelette complet quelconque.
- Pour les paires à moitié, l'équilibre ne porte que sur les D-classes **fermées**, c'est-à-dire sans paire à moitié au
  bloc, et sur leur complément.
- On ajoute les plafonds des orbites libres : φ_W(tm) ≤ #orbites de tm qui ne pivotent pas en t.

**Certificats exacts** (`ad/certify.py`).
- Au niveau PL : un Farkas rationnel, vérifié en exact.
- Sinon : un arbre de branchement dont chaque feuille porte un Farkas.
- Le contrôleur `certify_leaves.py --check` reconstruit tout depuis la feuille.

**Résultats.**
- **Classe 9, job 42** (10 feuilles au hasard sur 146 331) : **10/10 mortes**, 9 au niveau PL et 1 par un arbre à 7 nœuds.
  Tout est revérifié. Sans équilibre, les 10 vivent ; sans dichotomie, 3 seulement meurent.
- **M1 n = 6, job 12** : 1 feuille morte sur les 6 premières. C'est le « cœur », où G_c est connexe partout et l'équilibre
  vide. Le PL d'orbite de C2 y reste nécessaire.

## 6. Statut
- **Démontré à la main** : le lemme T⁺, et la mort des 14 survivants AD2.
- **Démontré, par calcul exact certifiable** :
  - AD pour n_dd ≤ 6 ;
  - l'échantillon de la classe 9.
- **Suite** : AD pour n_dd = 7 à 21, puis la famille dis générale (paires mixtes « un bloc + dis », tordues avec dis),
  qu'il faut énumérer.

## 7. Mise à jour (26/09, nuit) : résidus AD, (★) de C2, le cœur commun
**Balayage complet.**
- 26 524 classes AD de n_dd = 2 à 21.
- Mortes dans la version relâchée (structures libres en variables de partenaire) : toutes jusqu'à n_dd = 6.
- Survivantes :
  - n_dd = 7 : 6 ;
  - n_dd = 8 : 20 ;
  - n_dd = 9 : 14 ;
  - n_dd = 10 : balayage en cours ;
  - n_dd = 18 et 19 : vivantes ;
  - n_dd = 21 (« tout dd ») : inconnue.

**Retest en énumérant les 130 structures de A.**
- Les survivantes de 7 et 8 ont TOUTES un bloc A. La plupart des structures de A meurent. Il reste A en 6-cycle ou en 2+4.
- 10 classes sont exportées pour le modèle U de C2 (ad_residues_n7_8.json) ; le 2e passage en trouve d'autres.
- Tout bloc A force g⁻ = 20 (droite-bloc de rang 18, signe +1).

**Lemme (★) de C2 (équilibre général, pour tout ensemble S de blocs).** Re-dérivé ici, JUSTE :
Σ_{W∌k} δ_W² + Σ_{W∋k, W∉Lab} ε_W² = Σ_{o∈Lab} (2 − g_o)(g_o + 5).
- Preuve par racine : BBᵀ = 12I + 2J − A − Cls − A² sur N(ℓ_k), d'où
  Σ_z (n_C(z) − s)² = Σ_{v∈C} L_v(7 − L_v) − Σ_{v∉C} d_C(v)².
- Linéarisation exacte aux points entiers : tangentes de δ² et ε², et r_o ≤ min(10 − 4g, 12 − 6g). Code : ad/star_lemma.py
  et ad/star_generic.py.
- Résultats :
  - « tout dd » : faisable en PL, inconnu en entier à 300 s ;
  - résidus AD : en cours.

**Le cœur commun : l'étoile totalement libre.** Un bloc où aucune orbite ne pivote (μ ∩ τ = ∅). On la retrouve dans :
- les résidus AD (blocs à d = 0) ;
- le K₃ ancré de C1 (bloc 7) ;
- « tout dis ».

Rien de ce qui suit ne la tue :
- (iii) seul, faisable sur 6 tirages sur 10 ;
- (i)+(ii) seuls, faisable au moins une fois ;
- l'entrelacement à une ombre ;
- les coupes PSD d'orbite (𝒫² = 14𝒫) ;
- (★) en PL.

Feuille de problèmes pour les modèles extérieurs : astra/PROBLEMES_OUVERTS_INVOLUTION.md.
