# NOTE S24 : relecture du rapport d'Astra sur le dossier n° 5 (26/09)

Rapport copié tel quel dans `astra/RAPPORT_ASTRA_DOSSIER5.md`. Cadre : droite générique L (d(L) = 1), boule G = M_T
cubique sur 36 sommets (3 parties de 12), H = 3I − B, 60 lointains. M est l'incidence boule × lointains (10 par ligne,
2 par partie et par colonne). R = A_F − 3I − 3J + Z, avec Z = MᵀH⁺M (mod 7). Enfin, k est pair (L2, L4), L21 donne
k ≥ 22 et L23 donne k ≥ 24 si τ(G) est un carré non nul mod 7.

## §1. Une boule exclue à k = 22 : VÉRIFIÉ par un code indépendant
Script `scripts/verify_astra_d5_seven_shadows.py`, sortie dans `logs/verify_astra_d5_seven_shadows.log`. Il a son
propre Green mod 7 (résolution de HX = I − J, puis normalisation 1ᵀX = 0), τ exact (Bareiss) et sa propre recherche.
- **Boule.** G est cubique, sans triangle, et les trois réunions de deux parties sont des C₂₄ ; rang₇ H = 35.
  τ(G) = 1 928 719 172 295 ≡ 6 (mod 7), un non-carré : L23 ne s'applique pas.
- **Candidats.** Les nombres de candidats pour v = 1..6, 8 sont 350, 356, 338, 356, 307, 294 et 345, identiques à ceux
  d'Astra.
- **Recherche.** Aucun choix simultané n'existe (27 626 nœuds).
- **Contrôle positif** (v = 1..7). Un choix existe, et c'est exactement le témoin d'Astra. Revérifié séparément : ombres
  candidates, conditions par paire, étoile.
- **Modélisation relue.**
  - Les 7 lointains existent pour toute complétion : les A-ombres parcourent une fois chacune les 60 paires non
    partenaires.
  - Filtres : λ/μ contre la boule, K = MMᵀ (= nombre de voisins lointains communs) > 0 sur le support, et g(n,n) ≡ 6
    (diagonale de R).
  - Conditions par paire : g ∈ {2, 3}, a = 3 − g, et nᵀn′ ≤ 2 − a.
  - Étoile de u : Σ_w a ≤ 1 − (Bn_v)_u, parce que l'arête u y_v a exactement un voisin commun. Tout cela est correct.
- **Portée.** Toute SRG contenant cette boule vérifie k ≥ 24. Le résultat porte sur une boule, pas sur toutes.

## §2. R = 0 ⇒ SRG exacte : DÉMONTRÉ, preuve refaite à la main
- **Étape clé.** Les colonnes de H⁺M ont une somme nulle dans chaque partie.
  - B e_t = 1, car chaque sommet a exactement un voisin dans la partie t. Donc H e_t = 3e_t − 1, puis
    H⁺e_t = 3⁻¹(e_t − 12·1) ≡ 5e_t + 3·1 (mod 7).
  - Pour une ombre n (2 points par partie, 6 en tout), on obtient e_tᵀH⁺n ≡ 5·2 + 3·6 = 28 ≡ 0.
- **Conséquences mod 7**, en utilisant HH⁺ = I − J (avec 36 ≡ 1), H⁺1 = 0 et MMᵀ ≡ −H² − D + 2J :
  - MZ = −HM ;
  - Z² = J − MᵀM, puisque MᵀJM = 36J ≡ J ;
  - ZJ = 0.
- **Si R = 0**, alors A_F = 3I + 3J − Z. On en tire (B + I)M + MA_F ≡ 7M + 3·MJ ≡ 2J, et A_F² + A_F ≡ 5I + 2J − MᵀM,
  soit 12I + 2J − MᵀM.
- **Relèvement aux entiers.** F = (A² + A − 2J) sur boule × lointains et E = (A² + A − 12I − 2J) sur lointains × lointains
  ont leurs entrées (hors diagonale) ≥ −2 et ≡ 0 mod 7, donc ≥ 0. Leurs sommes de lignes sont nulles
  (40 + 80 − 120 et 64 + 8 − 12 − 120 + 60), donc F = E = 0.

## §3. MR = 0 et rang₇ R ≤ 2 ⇒ SRG exacte : DÉMONTRÉ, refait
- **E ≡ R².** On a R1 ≡ 0 (8 − 3 − 180 = −175) et RZ = (MR)ᵀH⁺M = 0.
- **tr R ≡ 0.** tr Z = −tr H − Σ_t e_tᵀH⁺e_t = −108 − 288, d'où tr R = −360 − 396 = −756 ≡ 0.
- **Cayley–Hamilton en dimension 2** donne T² = 0, donc rang R² ≤ 1. Une matrice symétrique de rang 1 à diagonale nulle
  est nulle, donc R² = 0.
- **Le contre-exemple de rang 4** est vérifié : rang 4, R1 = 0, tr R = 0, diag(R²) = 0 et R² ≠ 0.
- **Complément de « structures » : MR = 0 est automatique dans une vraie SRG.**
  - Q = A − 3I + 3J vérifie Q² ≡ 959J ≡ 0 (mod 7).
  - Avec la condition d'image qui rend la formule de Schur valide (celle de l'additivité des rangs de notes/102 §8), les
    blocs de Q² donnent Q_{B′F} R ≡ 0, puis MR ≡ −3JR = 0.
  - **Conséquence : à k = 24 aussi, le problème d'ombres est un problème de construction complet.**

## §4. Un échange de ≤ 6 ombres force deux prismes : DÉMONTRÉ, refait
- (A − 3I)z ≡ 0 mod 7. Hors du support, la valeur est dans [−6, 6], donc 0. Sur le support, les valeurs sont dans
  {3z, −4z}.
- L'égalité exacte Az = 3z₃ − 4z₋₄ force z₃ ∈ V₃ et z₋₄ ∈ V₋₄. Pour la partie −4, qC(4,2) ≤ 2C(q,2) donne q ≥ 7 > 6.
- Donc Az = 3z. Un graphe sur ≤ 6 sommets de degré minimal ≥ 3 avec λ ≤ 1 et μ ≤ 2 est le prisme :
  - K₄ et K_{1,2,2} ont une arête dans 2 triangles ;
  - K_{3,3} a μ = 3 ;
  - un sommet de degré ≥ 4 est exclu.
- Le degré interne vaut 3, donc il n'y a pas d'arête entre S et T.

## Bilan
Le rapport est entièrement correct. Il n'y a ni obstruction uniforme ni progrès sur les cas dégénérés.
