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
- Troisième point de B : c′. C4(B) = {P,Q,d,d′}. Dans le repère de B, R = {P,Q} a pour troisième point M, qui est
  une classe de B. R est donc une paire du couplage de la classe M de B, et le lemme 1 donne
  **c_R(z) = 1 − [a_M ∈ R]**, où a_M est l'atterrissage via M.
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

## 4. Relecture du lemme P-Ψ de C1 (26/09) : TIENT, étapes refaites à la main
**Énoncé (C1).** Au Fano tendu, pour toute paire y d'un bloc R = {X,Y}, avec W = Φ(X,Y) :
θ(z_X(y)) + θ(z_Y(y)) = 0, où θ(z) = s_R(atterrissage de z via W).

**Étape 0 : simplification de θ.** Pour z dans le bloc {X,W}, θ = s(a)·η se simplifie.
- Si b = Φ(W,a), alors {a,b} ∈ M_M(R), donc s(b) = s(a) et θ = s(b).
- Si b = Φ(X,a), alors {a,b} est une paire mixte, donc s(b) = −s(a) et θ = s(b).
- Dans les deux cas, θ = s_R(atterrissage via W).

**Étapes 1 et 2.**
- pairblk moins les 4 spéciaux donne Σ_{z gén. ~ y} c(z) = t(y).
- Chaque générique est libre vis-à-vis de R, avec c_R ≥ 1. La proposition de §2 donne donc ⟨t(y),Ψ_R⟩ = 0.

**Étape 3.** On a r(y)_{R′} = 12[R′ = R] + 8 − 2|R∩R′| − c_{R′}(y).
- Σ_{R′} Ψ_R(R′)(8 − 2|R∩R′|) = 8·11·Σs − 4Σs = 0, car Σ_{c∈C4} s(c) = 0.
- ⟨c(y),Ψ_R⟩ = s(a_y) + s(b_y) + 9(s(a) + s(b)) + 4(2 − 2) = 0, par le lemme 1 et parce que la paire {a,b} est mixte.

**Étape 4 : les spéciaux de y**, avec le lemme 1 dans le repère du spécial.
- **τ-spécial dans le bloc XW** (troisième point Y ; M_M(z) = M_Y(R), M_X(z) = M_X(R), M_W(z) = M_M(R)) :
  ⟨c,Ψ_R⟩ = s(a) + 9s(b) − 4s(b) + Σ_c s(c)g_{Yc}.
  - Si b = Φ(W,a), le dernier terme vaut 0, d'où s(a) + 5s(b) = 6θ.
  - Si b = Φ(X,a), il vaut 2s(b), d'où s(a) + 7s(b) = 6θ.
- **σ-spécial dans le bloc {X,A}**, A = a_y (atterrissage via X égal à Y, par réciprocité de σ_X ; troisième point
  A* = π_X(A)). Les deux atterrissages via A, W ou π_Y(A), donnent **−18 s(A)** (sommes détaillées au JOURNAL).
- **Somme des deux σ-spéciaux** : −18(s(a_y) + s(b_y)) = 0, car la paire est mixte.

**Étape 5.** ⟨t(y),Ψ⟩ = 0 − 6(θ(z_X) + θ(z_Y)) + 0, d'où le lemme. ∎

**Portée.** P-Ψ est bien plus fort que I_R = 0 : il donne une condition par paire. C1 observe que pour le σ du témoin
et pour 8 σ tendus frais, 2 à 6 droites sur 7 n'admettent aucun réglage de retournements. Il reste à le démontrer
pour tout σ tendu (invariant par droite).

## 5. Critère exact de P-Ψ par droite (26/09, maths pures ; conventions de C1, notes/119)
**Codage.** w_Z(u,v) = c_v ⊕ m_v·u, avec C0 = (0,0), C1 = (1,0), M0 = (0,1), M1 = (1,1). Les types sont pris par
rapport au troisième point du bloc.

**Bloc B = XY, de troisième point Z.** Étiquettes x_a et y_b ; flips via τ_x(Y_b) = Z_{b⊕f_x}. On pose
p_a = c_{x_a} ⊕ m_{x_a}f_{x_a} et q_b = c_{y_b} ⊕ m_{y_b}f_{y_b}.
- P-Ψ s'écrit p_a ⊕ q_b ⊕ m_{x_a}b ⊕ m_{y_b}a = 1 pour a, b ∈ {0,1}.
- La somme des 4 équations donne Σm = 0 (mod 2). C'est automatique : les comptes circulants donnent #M = 4 − 2m_{Z,L}.
- Les solutions forment une famille à un paramètre π_B : p_{x0} = π_B, p_{x1} = π_B ⊕ m_{y0}, q_{y0} = 1 ⊕ π_B et
  q_{y1} = 1 ⊕ π_B ⊕ m_{x0}.

**Droite L = {X,Y,Z}.** Chaque étiquette de L est dans 2 blocs de L et porte un seul bit.
- Type C dans un bloc : une équation affine en π_B.
- Type M dans les deux blocs : une équation qui lie π_B et π_{B′}.
- Type M dans un seul bloc : le flip est fixé par ce bloc, et seule reste l'équation du bloc où l'étiquette est de
  type C.

**Critère.** La droite admet un réglage de retournements si et seulement si le système affine en 3 inconnues
(π_XY, π_XZ, π_YZ) est cohérent. Il a une équation par occurrence C et une par étiquette MM.
- Bloc tout-C : il est vivant si et seulement si les deux étiquettes d'une classe ont le même type C, et celles de
  l'autre classe le type C opposé.
- Tout-M : on retrouve la condition de cycle de C1.

