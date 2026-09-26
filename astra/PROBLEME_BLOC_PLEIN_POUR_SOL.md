# Problème fini autonome : le « bloc plein » d'une involution de srg(99,14,1,2)

*Rédigé par « structures » le 27/09 (nuit), pour un modèle extérieur. Tout ce qui suit est démontré chez nous (NOTE_S29,
notes 129 à 132 de C1). La question est purement combinatoire et finie. Répondre par une preuve d'impossibilité (avec certificat
vérifiable) ou par une solution explicite (à vérifier par programme).*

## 1. Données
Étiquettes : 0, …, 11. L'étiquette v = 2i + s désigne le bloc i (i = 0..5) et le signe s (0 pour +, 1 pour −).
π(v) = v xor 1 (opposé dans le bloc).

Une **structure D** est une réunion de cycles sur les 6 blocs qui couvre chaque bloc une fois. Les « digones » (cycles de
longueur 2) sont permis. Il y a 4 types : 6-cycle, 4+2, 3+3, 2+2+2. Pour un cycle (p₁, …, p_ℓ), on pose
m(2p_j) = 2p_{j+1} + 1 (indices mod ℓ). m est l'involution sans point fixe correspondante ; aucune arête de m n'est une arête π.

- **Racines lointaines** : les 60 paires z = {a, a′} d'étiquettes de blocs différents (les arêtes de K₁₂ − π).
  σz = {π a, π a′}.
- **Paires non-m** : les 60 arêtes de K₁₂ − m.

## 2. Problème A (le cœur : Φ seul)
Trouver une bijection Φ : {racines lointaines} → {paires non-m} telle que :
- (ii) pour toutes étiquettes a, b : #{z ∋ a : b ∈ Φ(z)} = K₂[a,b] = 2 − [a=b] − [b=πa] − [b=m a] − [b=m πa] ;
- (iii) pour toutes b, d : #{z : b ∈ Φ(z), d ∈ Φ(σz)} = K₃[b,d] = 2 − [d=m b] − [d=π b] − 2[d=b] (sur les 60 z) ;
- (iv) pour z = {a, a′} et chaque a ∈ z : on n'a jamais à la fois a ∈ Φ(z) et πa ∈ Φ(σz), ni à la fois πa ∈ Φ(z) et
  a ∈ Φ(σz) ;
- (v) pour chaque cycle C de la structure D et chaque z : |{blocs de z} ∩ C| + |Φ(z) ∩ lab(C)| + |Φ(σz) ∩ lab(C)| = |C|, où
  lab(C) est l'ensemble des 2|C| étiquettes des blocs de C. Pour un 6-cycle, cette condition est vide.

Conséquences déjà connues, utiles comme contrôles : Φ(z) ∩ Φ(σz) = ∅ (K₃ diagonal nul), et les matrices de Gram des moments
d'ordre 2 sont PSD (aucune contradiction de ce côté).

**Statut chez nous.** Aucune solution n'est connue, pour aucun des 4 types. HiGHS et CP-SAT sont indécis à 100-280 s. Φ
équivariant sous la rotation d'ordre 6 du 6-cycle est infaisable. Personne n'a jamais exhibé de solution de (ii)+(iii) à ce
paramètre, pour aucun type de bloc.

## 3. Problème B (le secteur impair complet ; plus contraint)
Pour chaque σ-paire Z = {z, σz}, fixer un représentant z = {a₁, a₂}, et poser :
- l'ombre s_Z = 1_{Φ(z)} − 1_{Φ(σz)} ∈ {0, ±1}¹² ;
- d_Z = 1_z − 1_{σz} ;
- B_Z ∈ {0, ±1}⁶, signe de z sur chacun de ses deux blocs (+1 si l'étiquette est paire).

Trouver Φ comme en A, plus une matrice T symétrique 30 × 30 (entrées −1, 0, 1, diagonale nulle) telle que :
- (E1) Σ_W T_{ZW} s_W = −(s_Z + d_Z) pour tout Z ;
- (E2) Σ_W T_{ZW} B_W = −(s_Z(e+) − s_Z(e−))_{e=0..5} ;
- (E3) B_Z·B_W + s_Z·s_W + (T²)_{ZW} + T_{ZW} = 12 [Z = W].
  - En particulier, chaque ligne de T a exactement 6 entrées non nulles.

B est exactement « S² + S = 12I » pour la matrice impaire 49 × 49 de l'involution.

**Fait démontré (utile).** Si B a une solution, T a le spectre exact 3⁹ (−4)⁴ (−1)⁵ (−2)⁶ 1⁶, et la matrice
X = (B ; s) de taille 18 × 30 est de rang 17.

## 4. Ce qui compte comme réponse
- **Impossibilité de A ou de B pour un type** : un certificat vérifiable (Farkas exact sur un PL ou un arbre de branchement,
  ou un argument à la main). Pour B, cela tuerait tout bloc plein de ce type dans une involution.
- **Solution** : l'écrire en JSON (Φ comme liste de 60 couples, et T si B) ; nous la vérifierons.
- Pistes que nous n'avons pas essayées : symétries d'ordre 2 et 4 (le groupe d'automorphismes de (π, m) est d'ordre 12, 32,
  72 ou 384 selon le type) ; recherche par sommet (les 10 images des racines contenant a forment un graphe de degrés K₂[a,·],
  deux tels graphes H_a, H_{a′} partagent exactement une arête si a′ ≠ πa, aucune si a′ = πa, et toute paire non-m est dans
  exactement deux H_a).
