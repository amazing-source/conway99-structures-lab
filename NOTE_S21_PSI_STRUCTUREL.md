# S21 — Pourquoi Ψ_R annule les comptes permis : preuve à la main (Fano tendu, C1)

26/09/2026, labo « structures ». C'est l'item 1 de C1 (notes/118) : démontrer de façon structurelle que
⟨c(z;a,b), Ψ_R⟩ = 0 dès que c_R(z) > 0, sans contrôle exhaustif. **Mathématiques pures.** Le seul calcul est un
contrôle arithmétique de la forme close du lemme 1 (`scripts/fano_gtable_closed_form.py` : 168 cas sur 168).

## 0. Cadre
- Plan de Fano sur les 7 classes. Pour une paire de blocs R = {P,Q}, on note M = Φ(P,Q) et C4 = [7] ∖ {P,Q,M}.
- C4 est un quadrangle : il ne contient aucun triplet aligné. Ses 6 paires se répartissent en trois couplages M_P,
  M_Q et M_M, selon que le troisième point de la paire est P, Q ou M. Deux couplages distincts n'ont aucune paire en
  commun ; donc toute paire de l'un contient exactement un point de chaque paire de l'autre.
- Un lointain z du bloc B a pour atterrissages a (via une classe) et b (via l'autre). Ils sont permis si
  b ∈ {Φ(P_B,a), Φ(Q_B,a)}, c'est-à-dire si {a,b} ∈ M_{P_B}(B) ∪ M_{Q_B}(B), jamais dans M_M(B). **Donc a et b sont
  dans deux paires M_M différentes.**
- Signes : s = +1 sur une paire de M_M(R) et −1 sur l'autre. Les paires de M_P(R) et de M_Q(R) sont donc « mixtes »
  (s de somme nulle).
- Ψ_R = Σ_{c∈C4} s(c)[9e_{Mc} + e_{Pc} + e_{Qc}] + 4(e_{M_M⁺} − e_{M_M⁻}). Ψ_R est nul sur R, PM, QM et sur les
  paires M_P(R), M_Q(R).

## 1. Lemme (forme close des comptes au Fano tendu)
Pour z dans le bloc {P,Q}, avec atterrissages a (via P) et b (via Q), et a′, b′ les partenaires M_M de a, b :
- spéciaux : 1 sur PM, QM, Pa et Qb ;
- g = 2 sur les deux paires M_M ;
- g = 1 − [a ∈ paire] sur les paires M_P ;
- g = 1 − [b ∈ paire] sur les paires M_Q ;
- g(Ma) = [b = Φ(Q,a)], g(Ma′) = [b = Φ(P,a)], g(Mb) = [b = Φ(P,a)], g(Mb′) = [b = Φ(Q,a)].

*Preuve.* Comptes des familles tendues (2 / 1 / 3) et fentes (4 voisins par classe hors de P et Q).
- **Droite par P**, {P,x,y} avec {x,y} ∈ M_P : [a ∈ {x,y}] + g_{xy} = 1.
- **Droite par Q** : même raisonnement, avec b.
- **Fente de c ∈ C4** : [c=a] + [c=b] + g_{Mc} + (1 − α(c)) + (1 − β(c)) + g_{c,π_M(c)} = 4. Ici α(c) indique si a
  est dans la paire M_P de c, et β(c) si b est dans sa paire M_Q.
- **Paire {x,y} ∈ M_M.** La droite par M donne g_{Mx} + g_{My} + g_{xy} = 3. En sommant les deux fentes et en
  utilisant α(x) + α(y) = β(x) + β(y) = 1 (paires mixtes), on obtient g_{xy} = 3 − [a ∈ {x,y}] − [b ∈ {x,y}].
  Comme a et b sont dans deux paires M_M différentes, cela vaut 2.
- **Valeurs g(M·).** Elles se lisent ensuite sur les fentes de a et de b. ∎

## 2. Proposition
Pour tout z libre (blk(z) ∩ R = ∅) et tout atterrissage permis :
**⟨c(z), Ψ_R⟩ = 0 si c_R(z) > 0, et ±6 si c_R(z) = 0.**
Le bloc B de z est disjoint de R. Il y a 3 positions relatives (10 blocs = 2 + 4 + 4).

**(1) B est une paire de M_M(R)** (troisième point de B = M).
- Dans le repère de B, M_M(B) = {R, l'autre paire M_M(R)}. Donc c_R(z) = g = 2, toujours.
- ⟨c, Ψ_R⟩ = 9 + 9 (spéciaux B·M) + ([a ∈ R] + [b ∈ R]) − 9(g_{Mc3} + g_{Mc4}) − 4g_{c3c4} − Σ(4 paires mixtes).
- Les termes valent : [a ∈ R] + [b ∈ R] = 1, g_{Mc3} + g_{Mc4} = 1, g_{c3c4} = 2, et 2 pour la somme des 4 paires
  mixtes.
- Donc ⟨c, Ψ_R⟩ = 18 + 1 − 9 − 8 − 2 = **0**.

**(2) M ∈ B**, avec B = {M,c}, c ∈ C4(R), c′ son partenaire M_M(R) et {d,d′} l'autre paire M_M(R).
- Troisième point de B : c′. C4(B) = {P,Q,d,d′}. Dans le repère de B, R ∈ M_M… plus précisément R est la paire de
  troisième point M. Donc **c_R(z) = 1 − [a_M ∈ R]**, où a_M est l'atterrissage via M.
- On pose σ = s(c). Alors ⟨c,Ψ_R⟩/σ = 13 − 9[a_M ∉ R] + [b_c ∈ R] + (g_{c′P} + g_{c′Q}) − 4g_{dd′} − 5.
  - Le −5 vient des quatre paires {P,Q}×{d,d′}. Deux sont des paires M_M(B) (2 + 2), les deux autres sont les paires
    M_c(B) (0 + 1).
  - g_{dd′} = 1 − [a_M ∈ {d,d′}].
- Sur les 8 atterrissages : **0 quand a_M ∈ {d,d′}** (c_R = 1), et **6σ quand a_M ∈ {P,Q}** (c_R = 0).

**(3) Le troisième point de B est dans R** (par exemple B = {x,y} ∈ M_P(R), le cas M_Q étant symétrique).
- C4(B) = {Q, M, x′, y′}, avec x′ le partenaire M_M(R) de x, s(x′) = s(x) = σ, et s(y) = s(y′) = −σ.
- Dans le repère de B :
  - M_M(B) = {QM, x′y′} ;
  - M_x(B) = {Mx′, Qy′} ;
  - M_y(B) = {My′, Qx′} ;
  - R = {P,Q} est le bloc « M_B·Q », donc c_R(z) = g_{PQ} ∈ {0,1}.
- Les 8 atterrissages donnent ⟨c,Ψ_R⟩/σ :

| a | b | c_R | valeur/σ |
|---|---|---|---|
| Q | y′ | 0 | 6 |
| Q | x′ | 1 | 0 |
| M | x′ | 1 | 0 |
| M | y′ | 0 | 6 |
| x′ | M | 0 | −6 |
| x′ | Q | 0 | −6 |
| y′ | Q | 1 | 0 |
| y′ | M | 1 | 0 |

**Conclusion.** Dans les trois positions, la valeur est 0 exactement quand c_R(z) > 0, et ±6 sinon. ∎

## 3. Portée
- Cela démontre à la main que **I_R = 0 est nécessaire en tout sommet de Fano tendu**, pour chaque bloc R :
  - l'identité cellcross appariée à Ψ_R fait disparaître toutes les inconnues δ_R(z), car les z avec c_R = 0 ont
    δ = 0 ;
  - il reste ⟨rhs_R, Ψ_R⟩ = 0.
- Ce n'est **pas** une obstruction en soi. C1 a trouvé des retournements qui satisfont les 21 conditions (480
  affectations avec le σ du témoin). La contradiction éventuelle est ailleurs (items 2 et 3 de C1).
