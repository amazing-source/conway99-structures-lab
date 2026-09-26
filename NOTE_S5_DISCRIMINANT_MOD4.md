# S5 — Un critère 2-adique de non-existence pour les graphes fortement réguliers (discriminant mod 4 des réseaux propres)

*25/09, labo « structures », méthode C3 (D2 : généraliser, balayer les tables de Brouwer, ne jamais exclure un graphe
connu). Statut : **DÉRIVÉ + OBSERVÉ** — preuve courte, contrôlée sur de vrais graphes et sur toute la table de Brouwer
(v ≤ 1300). **Nouveauté à confirmer par un expert** : le critère n'apparaît ni dans la liste des conditions de
Brouwer–Van Maldeghem (*Strongly Regular Graphs*, 2022, §1.1 et table 1.2), ni dans les commentaires de la table en
ligne, où 35 des ensembles exclus sont marqués « ? ».*

## 1. Énoncé
Soit Γ un srg(v, k, λ, μ) à valeurs propres entières k > r > s, et d = r − s. Pour θ ∈ {r, s}, notons θ′ l'autre valeur
propre, m_θ la multiplicité de θ et Λ_θ = ℤ^v ∩ ker(A − θI).

**Théorème.** Supposons d impair et v₂(k − θ′) ≥ v₂(v) (automatique si v est impair). Alors :
1. m_θ est pair ;
2. (−1)^{m_θ/2} · det Λ_θ ≡ 1 (mod 4) ;
3. det Λ_θ = ∏_{q | d} q^{c_q} · P_θ, avec P_θ = ∏ q^{v_q(v)} sur les premiers q | v, q ∤ d, v_q(k − θ′) < v_q(v).

**Corollaire A.** Si de plus tout premier divisant d est ≡ 1 (mod 4), alors (−1)^{m_θ/2} P_θ ≡ 1 (mod 4).

**Corollaire B (famille infinie).** Si tout premier divisant n est ≡ 1 (mod 4) et si e est impair, il n'existe pas de
graphe de type Latin négatif NL_e(n) = srg(n², e(n+1), −n + e² + 3e, e² + e).
En effet v = n² n'a que des premiers divisant d = n, donc P = 1, alors que m_s = e(n + 1) ≡ 2 (mod 4).

**Corollaire C.** Si v est impair, que tout premier de v divise r − s et que tout premier de r − s est ≡ 1 (mod 4),
alors f ≡ g ≡ 0 (mod 4).

## 2. Preuve
1. **Pair.** Λ_θ ⊥ 1, donc pour x ∈ Λ_θ on a Σ x_i = 0, et x·x ≡ Σ x_i ≡ 0 (mod 2).
2. **Déterminant impair.**
   - E_θ = (A − θ′I)/(θ − θ′) − (k − θ′)J/((θ − θ′)v).
   - Si d est impair et v₂(k − θ′) ≥ v₂(v), E_θ est 2-entier, donc E_θ(ℤ₂^v) ⊆ Λ_θ ⊗ ℤ₂.
   - Or E_θ(ℤ^v) = Λ_θ^#, puisque ℤ^v est unimodulaire et Λ_θ primitif. Donc Λ_θ ⊗ ℤ₂ est unimodulaire.
3. **Réseau pair de déterminant impair.** Sur ℤ₂, c'est une somme de plans [[2a, b], [b, 2c]] avec b impair, de
   déterminant 4ac − b² ≡ −1 (mod 4). Le rang est donc pair et det ≡ (−1)^{rang/2} (mod 4), à un carré d'unité près
   (≡ 1 mod 8).
4. **Support premier du déterminant.** Le groupe discriminant est tué par le dénominateur de E_θ. Pour q ∤ d, avec
   W = ℤ_q^v ∩ 1^⊥ :
   - (A − θ′)/d est entier sur W, donc W = Λ_r ⊥ Λ_s sur ℤ_q ;
   - disc W ≅ ℤ/v est cyclique, donc un seul des deux porte toute la partie ℤ/q^{v_q(v)} ;
   - c'est Λ_θ exactement quand E_θ n'est pas q-entier, c'est-à-dire v_q(k − θ′) < v_q(v).
   D'où l'expression de det Λ_θ. Les q-parties, pour q | d, sont des puissances de q ; elles sont ≡ 1 (mod 4) si q ≡ 1 (mod 4). ∎

*Remarque.* Pour srg(99, 14, 1, 2), d = 7 ≡ 3 (mod 4) : le critère ne fait que forcer la parité du rang de recollement 7-adique
α = 14 + k, c'est-à-dire k pair (L2 de C1). C'est cohérent, et rien de neuf pour 99.

## 3. Contrôles
- **Vrais graphes, réseaux calculés exactement** (`scripts/genus_validate.py`, `scripts/disc4_validate.py`) :
  - Petersen, T(7), T(9), T(13), T(15), T(17), T(19), T(21), Paley(25), Paley(169), LS₃(5), LS₃(13), LS₄(13).
  - Partout : Λ pair, déterminant impair, parties cycliques et signes prédits, (−1)^{m/2} det ≡ 1 (mod 4).
  - Cas non triviaux :
    - T(7) : Λ₋₂ de rang 14 ≡ 2 (mod 4), sauvé par la partie cyclique 3 ;
    - T(15) : Λ₋₂ de rang 90, sauvé par 7.
- **Table de Brouwer complète** (HTML brut téléchargé, v ≤ 1300 ; `scripts/brouwer_crosscheck.py`) :
  - 4538 lignes : 1183 existantes, 2889 inconnues, 466 inexistantes.
  - **Aucun graphe existant exclu.** Les lignes complémentaires concordent (36 + 36).
- **Test décisif.** À v = 625 et 841, les graphes NL_e(n) connus (quadriques elliptiques VO⁻, van Lint–Schrijver,
  codes projectifs, Paley) ont **tous e pair**. Tous les e impairs sont « ? » dans la table, et le critère les exclut.
- **Résultat connu retrouvé** : srg(289, 54, 1, 12) = NL₃(17) (Bondarenko–Radchenko, preuve longue) tombe en cinq lignes.

## 4. Ensembles exclus (v ≤ 1300), statut actuel dans la table de Brouwer
Déjà connu inexistant : (289, 54, 1, 12).
**Marqués « ? », exclus par le théorème (35) :**
| famille NL_e(n), e impair | sporadiques |
|---|---|
| (169,42,5,12), (169,70,27,30) | (209,52,15,12), (235,52,9,12), (375,154,53,70), (375,170,85,70) |
| (289,90,23,30), (289,126,53,56) | (391,140,39,56), (430,39,8,3), (441,190,89,76), (494,85,12,15) |
| (625,130,15,30), (625,182,45,56) | (507,44,1,4), (616,205,90,57), (630,68,1,8), (690,318,121,168) |
| (625,234,83,90), (625,286,129,132) | (721,220,69,66), (750,210,55,60) [pg(14,14,4)?], (750,214,63,60) |
| (841,150,11,30), (841,210,41,56) | (780,369,158,189), (871,348,137,140), (875,184,33,40), (931,180,29,36) |
| (841,270,79,90), (841,330,125,132), (841,390,179,182) | (1026,451,212,187), (1148,481,210,195), (1248,203,22,35) |

Au-delà de 1300 (hors table), jusqu'à v = 5000 : 185 ensembles exclus en tout. Parmi eux, 63 sont de type NL_e(n) avec
e impair (n = 13, 17, 25, 29, 37, 41, 53, 61, 65) et 122 sont sporadiques (`logs/disc4_scan_5000.log`).

## 5. Limites et suite
- Pour les premiers q | d avec q ≡ 3 (mod 4), l'exposant c_q est inconnu ; le critère ne donne alors que sa parité.
  C'est le cas de 99.
- Le niveau mod 8 (formule d'oddité complète, `scripts/genus_scan.py`) n'ajoute aucune exclusion entière dans la classe
  « propre ». Il ne donne que des plafonds de type L22 (44 ensembles, dont 99 et srg(57, 14, 1, 4)).
- Extension possible : v pair avec une partie 2-adique cyclique (E_θ non 2-entier). Le type 2-adique y est encore
  calculable si e ≥ 2.
- **Avant toute annonce** : audit indépendant (Astra), puis avis d'expert (A. E. Brouwer tient la table).

## 6. Relecture de C2 (25/09) et compléments
- **Mécanisme confirmé** à partir des principes, y compris la règle du facteur de v. (k − r)(k − s) = μv et
  (k − s) − (k − r) = d, donc p ∤ d divise exactement l'un des deux facteurs, et disc(ℤ1) est cyclique.
- **(289, 54, 1, 12) et la famille NL** : recalculés à la main par C2, d'accord.
- **La condition de rang (m pair) ne mord jamais.** Sur les 6849 ensembles à multiplicités entières avec v et d impairs
  (v ≤ 3000, sans filtre de Krein), f et g sont toujours pairs : c'est une conséquence de l'intégralité. Les 36
  exclusions viennent toutes du discriminant mod 4.
- **r et s ensemble n'apportent rien de plus** : det_r · det_s · v est un carré, et (−1)^{(v−1)/2} v ≡ 1 (mod 4).
  Le mod 8 (Arf) n'ajoute rien non plus.
- **Filtre pour v pair** : on exige 2^{v₂(v)} | k − θ′. Vérifié par exemple sur (1248, 203, 22, 35), où k − θ′ = 224 = 2⁵·7
  et v₂(1248) = 5.
- **Antériorité.** C2 ne connaît pas l'énoncé et y voit une condition « Bruck–Ryser 2-adique ». Recherche bibliographique
  (web, livre de Brouwer–Van Maldeghem, articles récents de non-existence) : rien trouvé. À confirmer par un expert.
- **Extension (C2, point 5)** : pour p ≡ 3 (mod 4) divisant d, il faudrait connaître la parité de α_p par un autre
  argument (borne de p-rang). Rare, mais cela rouvrirait des cas.

## 5. Audit adverse par Astra (26/09) : le critère TIENT (`astra/RAPPORT_ASTRA_SUIVI4_ET_AUDIT_S5.md`)
- **Preuve : tient étape par étape.** Précisions à intégrer à la rédaction :
  - m_θ pair (pour que m_θ/2 ait un sens) ;
  - E_θℤ^v = Λ_θ^# par surjectivité (une forme entière sur un sous-réseau primitif s'étend à ℤ^v) ;
  - le « tout ou rien » des q | v, q ∤ d : une somme de deux q-groupes non triviaux n'est pas cyclique. A agit par k
    sur le discriminant de W et par θ sur la partie de Λ_θ, donc q^t | k − θ. De (k − r)(k − s) = μv et
    (k − r) − (k − s) = −d, exactement un facteur est divisible par q, avec valuation ≥ t.
- **Vrais graphes** (réseaux saturés, Bareiss) :
  - Petersen, Λ₁ de déterminant 162, hypothèse 2-adique fausse (elle ne doit donc pas être supprimée quand v est pair) ;
  - Petersen, Λ₋₂ de déterminant 405 ;
  - T(7), Λ₃ de déterminant 109 375 ;
  - T(7), Λ₋₂ de déterminant 46 875.
  Tout est conforme. Il y a aussi une formule analytique pour T(n), n impair.
- **Balayage** : deux implémentations (Python et JavaScript), avec les 26 pages de la table relues.
  - 4 538 lignes, 1 183 existantes, 0 existante exclue.
  - **35 « ? » exclus**, 0 manquant, 0 en trop.
  - **Deux inexistants connus retrouvés** : (289,54,1,12), qui est le cas de Bondarenko–Radchenko, et
    **(154,51,8,21)**, un contrôle supplémentaire que nous n'avions pas cité.
- **Antériorité.**
  - Le lemme de réseau (réseau pair de déterminant impair ⟹ (−1)^{m/2} det ≡ 1 mod 4) est **classique** : J. Martinet,
    *On parity classes*, notes de Lausanne (2003), corollaire 6.4.
  - Proches mais différents : Koyama 1991 (Bruck–Ryser–Chowla pour les srg d'ordre impair) et Bagchi (arXiv:2102.05576,
    designs quasi symétriques).
  - **La nouveauté du critère paramétrique et du corollaire NL_e(n) reste INDÉCISE** (non retrouvés dans une recherche
    bornée). À ne pas présenter comme nouveau sans l'avis d'un expert.
