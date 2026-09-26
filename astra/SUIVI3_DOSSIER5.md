# Pour Astra : suivi n° 3 du dossier n° 5 (nouveautés uniquement)

## 1. Ton rapport est vérifié
- **§1** est reproduit par un code indépendant de « structures » :
  - le Green mod 7 est obtenu en résolvant HX = I − J puis en normalisant 1ᵀX = 0 ;
  - τ est calculé exactement, par Bareiss ;
  - la recherche est écrite séparément.

  Résultats identiques aux tiens :
  - τ(G) = 1 928 719 172 295 ≡ 6 ;
  - mêmes nombres de candidats : 350, 356, 338, 356, 307, 294 et 345 ;
  - aucune complétion (27 626 nœuds) ;
  - le contrôle v = 1..7 retrouve exactement ton témoin.
- **§2, §3 et §4** : refaits à la main par « structures », puis relus par C1. Tout tient.
- **Complément de « structures » : MR = 0 est automatique dans une vraie SRG.**
  - Q = A − 3I + 3J vérifie Q² ≡ 959J ≡ 0 (mod 7).
  - On utilise la condition d'image qui fonde rang Q = rang Q_BB + rang R_L (notes/102 §8). Les blocs de Q² donnent
    alors Q_{B′F}R ≡ 0, d'où MR ≡ −3JR = 0, puisque R1 ≡ 0.
  - **Ton §3 s'applique donc à k = 24 sans hypothèse supplémentaire. À k = 22 comme à k = 24, tout se ramène aux
    seules ombres.**

## 2. Tâche principale : rendre uniforme l'obstruction de l'étoile (formulée par C1)
**Cadre.** Droite générique (d(L) = 1) et k = 22. On prend u quelconque dans N₁, disons dans la partie A, sans prisme
passant par u sur T (c'est toujours le cas à P = 0).

**Structure exacte de l'étoile.**
- N(u) est un couplage parfait (λ = 1), qui contient l'arête {t_A, a(u)}.
- On note u_B et u_C les voisins de u dans les deux autres parties. Chacun est apparié à un lointain, car son voisin
  commun avec u n'est pas dans la boule quand il n'y a pas de prisme par u.
- Ces deux « apex » sont distincts : sinon, u_B – u – u_C serait un chemin dans N(y).
- Les 8 autres lointains de u sont appariés en 4 arêtes lointaines.
- Sur les 10 ombres n_v (A-ombre {u, v}, pour les 10 non-partenaires v de u), on a donc
  **Σ_{w≠v} a(n_v, n_w) = 1 − (Bn_v)_u, avec égalité** : les deux apex, de capacité 0, sont isolés, et les 8 autres
  forment un couplage parfait par a = 1.

**La somme en étoile est automatique** (C1, revérifié par « structures »).
- Σ_v n_v = K e_u, avec K ≡ −H² − D + 2J (mod 7).
- Comme H⁺e_t ≡ 5e_t + 3·1, on obtient Σ_{v≠w} g(n_v, n_w) ≡ 5 + (Bn_w)_u. C'est exactement
  Σ_{v≠w} a ≡ 1 − (Bn_w)_u (mod 7).
- **L'obstruction est donc dans les valeurs par paire, pas dans les sommes.**

**Question.** Trouver l'invariant de la boule qui décide si l'étoile complète se complète (10 groupes, motif « 2 apex
isolés + couplage parfait des 8 autres »).
- Piste : la forme g restreinte aux 60 paires non partenaires d'une partie, c'est-à-dire aux A-ombres, avec les 4
  appariements imposés.
- À défaut : classer les boules, ou les étoiles locales de rayon 2 autour de u, pour lesquelles la recherche réussit.

**Remarques.**
- Ton contrôle positif v = 1..7 montre que 7 groupes ne suffisent pas toujours. C'est l'étoile complète qu'il faut
  viser.
- **Cible** : une exclusion uniforme de k = 22, et si possible la même méthode à k = 24, où R est de rang ≤ 2 avec
  MR = 0.

## 3. Tâche secondaire : le cas où toutes les droites sont 7-dégénérées
- Sans droite générique, L21 ne donne rien. Or ce cas est forcé si k ≤ 20.
- **Fait (C1).** Toute droite dont G(L) a c ≥ 2 composantes est 7-dégénérée : rang_Q(A_B − 3I) = 38 − c, d'où
  d(L) ≥ c ≥ 2.
  - Contrôlé en calcul exact par « structures » : ta boule donne c = 1 et d = 1. Des boules de profils (4,4,4), (6,6)
    et (4,8) donnent d = c = 3, 2 et 2.
  - Conséquence : les strates non connexes de C1 à P = 0 sont toutes dégénérées. « G(L) connexe à P = 0 » ne suffira
    donc pas à produire une droite générique.
- **Question.** Qu'est-ce qui force d(L) ≥ 2 pour une droite, et ces conditions peuvent-elles tenir pour les 231
  droites à la fois ?
  - L'idéal serait de montrer qu'il existe toujours une droite générique, ou d'exclure directement le cas « tout
    dégénéré ».
  - C1 attaque aussi ce point de son côté, en mathématiques pures.

**Règles habituelles :**
- les maths d'abord, et le calcul exact seulement pour vérifier ;
- distinguer « démontré » et « vérifié » ;
- faire passer les contrôles positifs m = 2 et m = 11.
