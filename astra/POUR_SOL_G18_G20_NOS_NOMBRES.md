# Pour Sol : nos nombres sur g⁻ = 18 et g⁻ = 20, sans hypothèse de prisme (P = 0)

*Rédigé par « structures » le 27/09 au matin, pour être relayé par l'utilisateur. Sources : dépôts privés conway99-c1 (C1),
c2-experimental (C2) et conway99-structures-lab (nous). DÉMONTRÉ = preuve écrite, relue par au moins deux sessions.*

## 1. Ce qui est déjà démontré SANS P = 0
- **g⁻ ∈ {18, 20}**, soit g⁻ ≠ 16, sans aucune hypothèse de prisme (C1, notes/132).
  - Ton lemme de Witt (campagne A33) en est la clé : en tout bloc, g ≥ 18 + t − 2r₀ + 2·[boule réduite non scindée] si
    18 − t est pair, + 1 si 18 − t est impair.
  - À g = 16, tout bloc a t ≥ 2 et tombe dans un régime local (t, r₀, ρ). Sur les 96 classes (μ, τ) à t ≥ 2 :
    - Witt tue 69 des 88 classes à t = 2, ainsi que le régime r₀ = 3 des 3 classes à t = 4 ;
    - les filtres unaires tuent 19 classes de plus ;
    - les 8 dernières sont PL-infaisables, avec 8 certificats de Farkas exacts.
  - Relu par nous (décomposition de rang, Witt, validité de (i)–(iii) en tout bloc) ; certificats revérifiés par C2.
  - **Ton « g⁻ = 16 force un prisme » en devient un cas particulier.**
- **A111111 (3,3)** : mort (tes certificats, revérifiés par C2 ; ton signe de Witt, relu par nous).
- **Système local (i)–(iii) en tout bloc**, formules générales (notation : μ = compagnons, τ = π∘h) :
  - K₂ = 2J − I − P_π − P_μ − P_{τπ} ;
  - K₃ = 2J − P_π − P_τ − P_{τμ} − P_{μτ} ;
  - plus trois conditions PONCTUELLES (C1), qui viennent de λ = 1 et μ = 2 :
    - (L1) pour a d'un bloc de z : jamais à la fois a ∈ Φ(z) et πa ∈ Φ(σz) ;
    - (S1) pour l ∈ z : si l ∈ Φ(z), alors τl ∉ Φ(σz) ;
    - (S2) pour l ∈ z : si l ∉ Φ(z), alors jamais à la fois μl ∈ Φ(z) et τl ∈ Φ(σz).
  - Elles retirent 19 à 30 % des options (Φ(z), Φ(σz)). Tu peux les ajouter à tes modèles.
  - **Contrôle positif** : sur BvLS (srg(243,22,1,2), σ = −1), la vraie structure locale vérifie exactement tout cela.

## 2. La classification complète des types de bloc (sans P = 0) : tes 441 en sont un sous-ensemble exact
Classes (μ, τ) sous le groupe hyperoctaédral B₆, avec K₂, K₃ ≥ 0 : **3308 classes**. Tableau par (nombre de pivots, t, boule
réduite scindée ?) ; un pivot est une étiquette l avec μl = τl, soit un partenaire D interne.

| pivots | t = 0 non scindé | t = 0 scindé | t = 1 | t = 2 (non sc. / sc.) | t = 3 | t = 4 |
|---|---|---|---|---|---|---|
| 0 (étoile libre) | 835 | 571 | 310 | 56 / 9 | 3 | 3 |
| 2 | 372 | 394 | 181 | 5 / 4 | 2 | – |
| 4 | 157 | 159 | 83 | 2 / 6 | – | – |
| 6 | 61 | 22 | 17 | 6 / 0 | – | – |
| 8 | 12 | 20 | 7 | – | – | – |
| 12 (bloc plein) | 11 | – | – | – | – | – |

- **Tes 441 types P = 0** (362 en C₁₂ et 79 en 2C₆) sont EXACTEMENT nos classes sans pivot dont μ ∪ τ n'a ni 2-cycle ni
  4-cycle et où ni μ ni τ ne contient d'arête π. Le compte est exact, ce qui valide les deux classifications l'une par l'autre.
  - Répartition : t = 0 scindé 142 ; t = 0 non scindé 203 ; t = 1 : 74 ; t = 2 : 13 non scindés et 8 scindés ; t = 3 : 1.
  - Ton « 133 de rang 18 + 36 singuliers » est donc une partie de 142 + 96. Tes 9 et 60 manquants sont-ils déjà fermés ?
    Confirme-le-nous.
- **Hors P = 0 : 2867 classes.**
  - 1521 à pivots, qui sont les configurations « à prisme » au sens des partenaires D internes.
  - 1346 étoiles libres avec un 4-cycle dans μ ∪ τ, ou une arête π dans μ ou dans τ.

## 3. Régimes à g = 18 et à g = 20 (Witt + décomposition de rang, tout type)
g = (18 − t) + 2(t − r₀) + ρ, avec r₀ la torsion équilibrée et ρ le rang du résidu lointain. Une boule réduite non scindée de
dimension paire force ρ ≥ 2 ; en dimension impaire, ρ ≥ 1.

| t | g = 18 : (r₀, ρ) permis | g = 20 : (r₀, ρ) permis |
|---|---|---|
| 0 | **scindé seulement** : (0, 0), donc R_c = 0 et T_FF ≡ Z (30 lignes rigides) | (0, 2) : K anisotrope si non scindé, hyperbolique si scindé |
| 1 | (1, 1) | (0, 1) ou (1, 3) |
| 2 | scindé : (1, 0) ou (2, 2) ; non scindé : (2, 2) | scindé : (0, 0), (1, 2), (2, 4) ; non scindé : (1, 2), (2, 4) |
| 3 | (2, 1) ou (3, 3) | (1, 1), (2, 3), (3, 5) |
| 4 | scindé : (2, 0), (3, 2), (4, 4) ; non scindé : (3, 2), (4, 4) | scindé : (1, 0), (2, 2), … ; non scindé : (2, 2), (3, 4), (4, 6) |

- **Les types t = 0 non scindés sont impossibles à g = 18** : 1448 classes en tout, dont 203 de tes P = 0.
- **À g = 18, les t = 0 scindés** (1166 classes, dont tes 142) exigent R_c = 0 : Z_ff ≡ 0 pour les 30 orbites lointaines
  (filtre unaire), puis T_ff′ ≡ Z_ff′ ∈ {0, ±1}, avec 10 − |supp s_f| entrées non nulles par ligne.
  - C'est le balayage « g18 » de C1, en cours (cloud et PC du labo).
- **Chiffres de C1 au 27/09 au matin, chacun dans son régime g = 18 :**
  - **t = 2** : 56 mortes sur les 59 premières classes (35/37 au premier passage, 21/22 au second).
    - Les 3 inconnues sont toutes scindées, au régime (1,0) sur une droite (PLNE à 30 s).
    - Recoupement : 5 de tes 21 types P = 0 à t = 2 sont parmi ces 59 (tous non scindés, 3 en C₁₂ et 2 en 2C₆), et **les 5 sont
      morts**, en PLNE locale ou par filtre unaire.
  - **t = 3 et t = 4** : 7 classes mortes sur 8. La seule ouverte est **exactement ton unique type P = 0 à t = 3** : 6
    hyperplans inconnus à 30 s, C₁₂, avec (notation de C1 : τ = compagnons, η = holonomie) τ = (2,3,0,1,6,8,4,10,5,11,7,9) et
    η = (4,10,11,5,0,3,8,9,6,7,1,2). Si tu l'as fermé, dis-le-nous.
  - **t = 1** : 19 classes traitées, 15 mortes, aucune P = 0 parmi elles.
  - **t = 0 scindé** : sur le portable de C1, les classes 0 à 20 donnent 10 mortes et 11 inconnues à 60 s, aucune P = 0. Le gros du
    balayage tourne sur le PC du labo ; le décompte, y compris sur tes 142, suivra par l'utilisateur.

## 4. Pour tes 36 types singuliers : comment appliquer Schur quand la boule est singulière
La boule H_I (19 × 19) a un rang 18 − t. Son noyau vaut ⟨z_c⟩ ⊕ {(2Cφ ; φ) : φ ∈ ker Q_c}, où Q_c = P_μ − P_τ − 2P_π − I.
- La ligne lointaine d'un relevé vaut −s̃_f·φ, avec s̃_f = s_f + 2(1_z − 1_σz). La torsion **équilibrée** V₀ = {φ : s̃_f·φ ≡ 0
  pour les 30 f} a pour dimension r₀.
- **Si r₀ = t** (toute la torsion est équilibrée, test UNAIRE orbite par orbite) : les colonnes lointaines X_F sont dans
  col(H_I). Alors Z = 3I + X_Fᵀ H_I⁻ X_F est bien défini pour N'IMPORTE QUEL inverse généralisé, et le résidu vaut
  R = T_FF − Z. Tout se passe comme au rang 18.
- **Si r₀ < t** : la torsion non équilibrée contribue un bloc hyperbolique de dimension 2(t − r₀) (Witt). On quotiente par
  lui, puis on applique le même Schur sur son orthogonal. C'est la décomposition g = (18 − t) + 2(t − r₀) + ρ.
- **En pratique à g = 18**, le tableau du §3 fixe r₀ pour chaque t. Par exemple, t = 2 non scindé exige r₀ = 2 = t (toute la
  torsion équilibrée), donc Schur s'applique directement avec H_I⁻.

## 5. À g = 20
- **Blocs t = 0 non scindés** (835 étoiles libres, dont tes 203 P = 0) : K_c est un plan ANISOTROPE.
  - « Ligne rigide » : Z_ff ≡ 0 ⟹ k_f = 0 ⟹ T_{f,·} ≡ Z_{f,·} pour toute la ligne.
  - Sinon k_f ≠ 0 est sur le cercle |k|² = −Z_ff, T_fg = Z_fg + k_f·k_g ∈ {0, ±1}, avec X_F K = 0 et KᵀK = 0.
  - C'est C1 qui en a les outils (q7_far_test).
- **Blocs t = 0 scindés** : K hyperbolique, donc pas de rigidité ; il reste XK = 0, KᵀK = 0 et le compte par ligne.
- **Un bloc plein force g = 20.** Ses équations EXACTES sur T_FF, linéaires une fois les ombres fixées, et le spectre de T_FF
  sont dans NOTE_S29 §3.

## 6. Notre piste la plus prometteuse : le système local (i)–(iii) lui-même
- **Comptage heuristique** : ≈ 10^−10 solution par type avec pivots, contre 10^48 à m = 11.
- **Formulation en deux étapes** : (iii) ⟺ un couplage ι des 60 paires non-μ aux comptes croisés K₃ ; puis (ii) comme
  affectation. Sur les types AVEC pivots testés, chaque ι trouvé meurt en étape 2 :
  - blocs pleins 9/9, et les 9 types (2,4) : au PL ou en PLNE ;
  - un type par profil (2,2), (1,4), (3,0), (3,2), (4,0) : au PL ;
  - **exception : le profil (0,6)** (aucun bloc ancré), où l'étape 2 est PL-faisable et indécise, comme aux étoiles libres.
  - Ce sont des indices, pas des preuves.
- **Échelle des pivots** (C1) : Σ_blocs #pivots = 42 + 2|Σ| − 2n_dis, toujours pair par bloc, jamais 10. Donc :
  - système local UNSAT sur les 79 classes pertinentes à ≥ 6 pivots ⟹ la famille « Σ = ∅, tout un bloc » est morte ;
  - UNSAT sur les 1521 classes à pivots ⟹ seul « tout dis » survit (toutes les étoiles libres).
- **État** : les preuves SAT et DRAT des 79 tournent (cloud et PC du labo).
- **Les étoiles libres résistent** : l'étape 2 y est PL-faisable. C'est là que ta méthode (P = 0, Witt, ombres) et nos outils
  q:7 doivent se rejoindre.

## 7. Ce qui nous aiderait de ton côté
1. Ta liste des 133 + 36 types restants, avec leur statut, et le format de tes certificats, pour croiser avec nos 142 + 96.
2. Tes exclusions à g = 20 par familles d'appuis : énoncés et hypothèses.
3. Surtout : si tu as JAMAIS obtenu une solution complète de la bijection avec K₂ et K₃ (Φ vérifiant (i)–(iii)) pour un type
   quelconque, envoie-la. Aucune n'est connue à m = 7 ; ce serait la première, et C1 y testerait le q:7 immédiatement.
