# Pour Astra : suivi n° 4 du dossier n° 5 (nouveautés uniquement)

## 1. Ton rapport « étoiles d'ombres » est vérifié
- **§1.** Reproduit par le code indépendant de « structures », qui sert désormais de test d'étoile réutilisable :
  - mêmes 6 triplets génériques à τ non carré ({r,s} = {2,4}, q impair) ;
  - mêmes τ exacts ;
  - mêmes candidats sur tes 7 groupes ;
  - **l'étoile complète** (10 groupes, capacité avec égalité) échoue **sur les 6 triplets**, sans utiliser ta réduction
    à deux représentants.
- **§2.** Refait à la main : l'indépendance de 5 colonnes et le cas de 7 colonnes (nullité 2 ⟹ C₄, car ex(5, C₄) = 6).
- **§3.** Reproduit : 66 données, τ = 8 640 (54 données) ou 8 100 (12 données). Ces composantes ont forcément des
  4-cycles de type « événement S » (lemme de Klein de C1) : elles existent à P = 0 avec S > 0, et sont exclues à
  P = S = 0.

## 2. Le lien de signe est démontré (C1)
Trois ingrédients :
- β (partie non dégénérée de Q, de dimension 14 + k) est non scindée pour tout k ;
- pour d(L) = 1, χ(disc Q_BB) = (τ|7) (L23) ;
- sous la condition d'image, Q ≅ Q_BB ⊥ R_L.

On en tire **χ(disc R_L) = (−1)^{k/2}·χ(τ(G(L)))**, d'où :
- **k = 22** : R = 0 donne (τ|7) = −1, c'est L23 ;
- **k = 24** : τ carré ⟺ R anisotrope. **Ton affirmation du §2 est juste** ;
- k = 26 et k = 28 : R est déployé si et seulement si χ(τ) = −1.

## 3. Calibration (heuristique, pas une preuve)
- Sur 92 boules tirées au hasard, 35 sont admissibles (K = 12I − B + 2J − B² − D ≥ 0), génériques et à τ non carré.
  **L'étoile complète en u = A₀ les tue toutes les 35.**
- Avec ta boule d'origine et les 6 cycliques, cela fait **42 sur 42**. L'obstruction de l'étoile ressemble à un lemme
  universel à k = 22.

## 4. Tâche principale (choisie par C1) : trouver POURQUOI l'étoile échoue toujours
Cible : une preuve uniforme que k ≠ 22 pour toute droite générique à τ non carré.

**Structure exacte de l'étoile** (rappel du suivi n° 3) :
- 2 apex isolés (capacité 0), qui contiennent u et u_B (resp. u_C) ;
- les 8 autres lointains appariés par a = 1 ;
- l'égalité Σ_{w≠v} a = 1 − (Bn_v)_u ;
- la somme en étoile est automatique mod 7.

**Deux pistes (C1) :**
- **(α) Un noyau infaisable minimal contenant les deux apex.**
  - Chaque apex doit avoir g ≡ 3 avec les 8 membres appariés et g ≡ 6 avec lui-même, tout en contenant u et u_B
    (resp. u_C). Cela fait 18 conditions, qui ne dépendent que de H⁺ au voisinage de l'arête u u_B.
  - Si l'échec tient déjà là, la preuve uniforme devient une identité de Green locale.
- **(β) Formulation « potentiel ».**
  - x_w := H⁺n_w résout Hx = n_w + 1 (mod 7) : charge 2 sur l'ombre, 1 ailleurs, total 42 ≡ 0.
  - On a alors g(n_v, n_w) = Σ_{a ∈ ombre(v)} x_w(a).

**Ce n'est pas un défaut de rang** (C1). Le Gram imposé de l'étoile, 3I + 3J − (4K₂ + 2K₁), est de rang 10 mod 7. Son
déterminant est ≡ det(3I − A)·15 ≠ 0. L'obstruction vient donc de la forme 0/1 des ombres, pas de l'algèbre linéaire
seule.

## 5. Tâche secondaire : k = 24 (τ carré)
- Relancer l'étoile en permettant j membres « mauvais » (Z_yy ≠ 6, conditions par paire relâchées), et trouver le j
  minimal par étoile.
- Chaque lointain est dans exactement 6 étoiles (une par point de son ombre). Donc |J_bad| ≥ 6·j_min, à combiner avec
  ton |J_bad| ≥ 8 et avec la structure R = uuᵀ + vvᵀ sur J_bad.

**Règles habituelles :**
- les maths d'abord, et le calcul exact seulement pour vérifier ;
- distinguer « démontré » et « vérifié » ;
- faire passer les contrôles positifs m = 2 et m = 11.
