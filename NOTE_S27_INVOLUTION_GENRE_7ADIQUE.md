# NOTE S27 : q:7 dans le cadre de C2 : genre des réseaux propres sous l'involution

26/09, soir. L'utilisateur demande des paris. La piste q:7 est proposée par C1 : C1 cherche les minorants,
« structures » traite le genre.

Cadre : σ est une involution de Γ = srg(99,14,1,2), avec Fix(σ) = {x}. Elle retourne les 7 blocs. Il y a 42 orbites de
racines. Code dans `involution/`.

## 1. Partie impaire : g⁻ = rang₇(S − 3I) est PAIR, et g⁻ ≤ 20
**Réseau.**
- L⁻ = {f ∈ Z⁹⁹ : σf = −f}, de base d_o = e_o − e_σo (7 blocs et 42 orbites).
- Avec la forme standard divisée par 2, c'est Z⁴⁹ standard : unimodulaire, impair, de vecteur caractéristique 1.
- A y agit par S, symétrique, avec S² + S = 12I, de valeurs propres 3²⁷ et (−4)²².
**Recollement.**
- Λ₃⁻ et Λ₋₄⁻ se recollent le long de (Z/7)^{g⁻}, et g⁻ = rang₇(S − 3I) puisque (S − 3I)Z⁴⁹ = 7Λ₋₄^#.
- det Λ₃⁻ = det Λ₋₄⁻ = 7^{g⁻}.
**Formule d'oddité pour Λ₋₄⁻** (défini positif, de rang 22) : 22 + 6g⁻ + 4[ε⁻ = −1] ≡ |E₋₄1|² (mod 8).
- L'oddité₂ est donnée par la projection du vecteur caractéristique, et |E₋₄1|² = (147 − 1ᵀS1)/7.
- **Calcul de 1ᵀS1** dans la jauge de C2 (représentants (a,+)(b,ε)) : 1ᵀS1 = −7 + 84 + 420 − 4M = 497 − 4M, où M est le
  nombre de paires {o,o′} avec Q = 1 et S = −1.
  - Cette valeur est invariante de jauge mod 8.
  - M est pair : cela découle aussi directement de 1ᵀ(S² + S)1 = 588.
**Conséquences.**
- **g⁻ est pair**, et [ε⁻ = −1] ≡ M − g⁻/2 (mod 2).
- **g⁻ = 22 est impossible.** Λ₋₄⁻ serait 7·(unimodulaire) de déterminant exactement 7²², ce qui force ε⁻ = +1, alors que
  la formule donne ε⁻ = −1.
- Avec g⁻ ≥ 7 (C2), on obtient **g⁻ ∈ {8, 10, …, 20}**.
- **Positivité dans la jauge de C2** : |E₋₄1|² = (4M − 350)/7 ≥ 0 et |E₃1|² = 99 − 4M/7 ≥ 0, donc **88 ≤ M ≤ 173**.

## 2. Partie paire : g⁺ est pair, et g⁺ ≠ 0
- L⁺ a pour base e_x, b_c, r_o, de Gram diag(1, 2I₄₉). E₁₄ = 1·1ᵀG⁺/99 et E₃ = (T + 4)(T − 14)/(−77) sont 2-entiers, donc
  L⁺ ⊗ Z₂ = Z₂1 ⊥ Λ₃⁺ ⊥ Λ₋₄⁺, avec pour Jordan 1¹ ⊥ 2⁴⁹.
- Donc Λ₋₄⁺ = 2·(unimodulaire) sur Z₂. On note Λ̃₋₄⁺ sa version de forme divisée par 2.
- Sur Z₃, Λ₋₄⁺ = 1²¹ ⊥ 9¹ (recollement Z/9 avec 1, de norme 99), d'excès 0. Il n'y a pas de partie 11-adique.
- **Oddité₂ de Λ̃₋₄⁺.** 1^⊥/2 ≅ Z₂⁴⁹ de Gram I + 2J, de caractéristique 1₄₉, ce qui correspond à v_c = 1 − 99e_x. Sa
  projection est −99E₋₄e_x, avec |E₋₄e_x|² = 4/9, d'où une demi-norme 2178 ≡ 2 (mod 8).
- **Formule** : 22 + 6g⁺ + 4[ε⁺ = −1] ≡ 2, donc **g⁺ est pair** et [ε⁺ = −1] ≡ 1 − g⁺/2. Un constituant vide aurait ε = +1,
  donc **g⁺ ≠ 0**.

## 3. Recollement global (L2 et L3 de C1 retrouvés)
- Λ₋₄ ⊗ Z₃ = 1⁴³ ⊥ 9¹, et la caractéristique projetée est nulle (E₋₄1 = 0). On en tire que k est pair et que
  [ε_g = −1] ≡ k/2.
- 7-adiquement, Λ₋₄ = Λ₋₄⁺ ⊥ Λ₋₄⁻ (χ₇(2) = 1), donc ε_g = ε⁺ε⁻ et **g⁺ + g⁻ = 14 + k**. Cela redonne M pair : c'est
  cohérent.
- On retrouve aussi k ≤ 22 + 20 − 14 = 28, c'est-à-dire L22.

## 4. Contrôles positifs de la machinerie (`involution/oddity_control.py`)
Sur de vrais graphes avec l'involution x ↦ −x, la formule « rang + excès_p ≡ |E_s1|² (mod 8) » est vérifiée exactement :
- **H(2,3)** : Λ₋₂⁻ de rang 2, det 9, Jordan 3-adique 3^{+2}, 6 ≡ 6 ✓ ;
- **BvLS** : Λ₋₅⁻ de rang 55, Jordan 3-adique 1^{−10} 3^{−10} 9^{+35}, 7 ≡ 7 ✓.

## 5. Minorant par les droites-blocs (`involution/blockline_rank.py`, `blockline_kernel.py`)
- La partie impaire de la boule de L_c est de dimension 19 : d_c, les 6 autres blocs et les 12 orbites de racines
  contenant c. S y est déterminé localement :
  - S_{o_v,o_w} = [τ_u v = w] − [ι v = w], où ι = mate ∘ σ_U est une involution, par σ-symétrie de la boule.
- **Rang₇ de S − 3I sur ces 19 coordonnées : toujours ≤ 18.** Sur 40 000 configurations aléatoires, il vaut 18 dans
  environ 81 % des cas, et descend jusqu'à 14.
- **La déficience est structurelle.** z_c = 3d_c + Σ_{o∋c} d_o est un vecteur propre EXACT de S de valeur 3. En effet,
  en dehors de la boule, μ = 2 avec (c,±) donne 2 − 2 = 0. De même, w_c = −4d_c + Σ_{o∋c} d_o est un vecteur propre de
  valeur −4.
  - Les 7 vecteurs z_c sont orthogonaux entre eux, de norme 21, et sont dans Λ₃⁻. Les 7 vecteurs w_c sont orthogonaux,
    de norme 28, dans Λ₋₄⁻.
  - On a z_c ≡ w_c mod 7 : ils fournissent les 7 éléments de recollement du minorant g⁻ ≥ 7 de C2.
- **Conséquence** : g⁻ ≥ max_c rang(N_c). Dès qu'une droite-bloc est « générique » (rang 18), **g⁻ ∈ {18, 20}**.

## 6. Ce qui tuerait l'involution
- Un minorant **g⁻ ≥ 21**. Les sous-matrices sans entrée « lointaine » plafonnent à 18 : il faut la structure lointaine
  (les squelettes de C2, où Q est connu et seuls les signes lointains sont inconnus).
- Ou une contradiction propre aux cas g⁻ = 18 et g⁻ = 20. Leurs signes sont fixés : [ε⁻ = −1] ≡ M − 9 et M − 10, avec
  M pair.

## 7. Le signe local à g⁻ = 18 (`involution/blockline_sign.py`)
- Si g⁻ = 18 et qu'une droite-bloc E a un rang 18, alors E/rad_E ≅ V/rad(B) comme espaces quadratiques. Le rang
  permet de le déduire : l'image de E dans V/rad est de dimension ≥ 18 = g, et E ∩ rad = rad_E.
- Donc ε⁻ = χ((−1)^g disc B) = χ(disc B|_E), le déterminant du bloc 18×18 obtenu en retirant d_c, qui est la
  coordonnée du radical z_c.
- La genèse impose ε⁻ = −1 (M pair, g = 18).
- **Local** : sur ~33 000 configurations de rang 18, le signe vaut −1 dans environ 45 % des cas et +1 dans environ 55 %.
- **Contrainte** : à g⁻ = 18, **toute droite-bloc générique doit avoir le signe local −1**, et les 7 droites-blocs de x
  doivent concorder.
- **Cohérence avec L7 de C1** : la configuration non tordue (τ = ι = mate) donne +1. Elle est donc interdite à
  g⁻ = 18, ce qui rejoint « un point de Hamming n'est pas le point fixe d'une involution ».
- **g⁻ = 20** : aucune droite-bloc ne détermine B (rang 18 < 20). Il faut la structure lointaine.
- **Cas restants pour l'involution** :
  - (a) g⁻ = 20 ;
  - (b) g⁻ = 18, avec toutes les droites-blocs génériques de signe −1 ;
  - (c) g⁻ ≤ 16, avec toutes les droites-blocs de rang ≤ g⁻, ce qui est rare localement (~1,5 % des configurations
    aléatoires ont un rang ≤ 16).

## 8. Mise à jour du 26/09, soir : g⁻ ∈ {16, 18, 20}
**Astra, dossier 6** (relu : `astra/RAPPORT_ASTRA_DOSSIER6.md`).
- min Λ₋₄⁻ ≥ 8 : preuve combinatoire, avec D² = 2I + 2J et 4 blocs contre Cb = −3a.
- det L_T = 2¹²·7^{g⁻−7}.
- Blichfeldt exclut g⁻ = 8 (rang 22) et g⁻ = 10 (L_T, rang 15).
- Réduction à un N_max de déterminant 49, à discriminant anisotrope.

**C1, recoupé par nous** (`involution/blockline_exhaustive.py`, `blockline_t4_digons.py`).
- **t_c = dim ker Q_c ≤ 4**. On énumère les 11 classes de τ et les 10 395 ι sans point fixe (ι = mate∘ρ est sans point fixe
  par F0). Histogramme des rangs de la droite-bloc : {14 : 800, 15 : 124, 16 : 11 363, 17 : 15 691, 18 : 86 367}.
- Les rangs 14 sont exactement les cas où τ a 2, 4 ou 6 digones avec mate et ι aucun.
- Cohérence : ι_c a un digone en d ⟺ τ_d a un digone en c. C'est la même arête : une racine voisine de sa sœur qui partage (d,s).
- Donc g⁻ = 14 est impossible : tous les blocs seraient de rang 14, un digone de τ_c en d forcerait ι_d à en avoir un en c,
  et le bloc d aurait un rang ≥ 15.

**Bilan : g⁻ ∈ {16, 18, 20}.**
- g⁻ = 16 exige de la torsion équilibrée dans les 7 blocs, par g⁻ = 18 + t_c − 2r₀_c + ρ_c.
- 18 et 20 demandent de l'information globale.
