# S38 — État au 27/09 ~14 h 30 et pistes pour conclure vite (structures)

## Ce que je fais en ce moment (instantané)
- **Route (C)** : 156 enveloppes explicites (S35). Exclues : **19/156** (piège, packing, (E), décomposable, Sol).
- **O₂₄** : 141 classes (y₁,y₂).
  - **86 exclues** par PL de classe entière (certificats entiers, vérificateur autonome).
  - 55 sont dans la DFS de C2 : 0 et 12 sont EMPTY. La classe 0 est donc exclue par deux voies indépendantes.
- **Étage PL général** `frame_lp.py` + `run_frame_lp.py` (tranches parallèles) + `verify_frame_cert.py`.
  - **132** : paires complètes (C1) = 2416. Premières 154 : 144 exclues, 10 PL-faisables. Le lot complet tourne au PC (C2).
- **Réseaux pauvres en racines** : des millions de paires, il faut Γ = Aut(U)/W.
  - **id 1 (a₁⁸)** est construit dans le repère de Golay (`golay_family.py`) : a ∈ ℤ²⁴, norme a·a/8. G₁ = Stab_M24(octade K) × signes.
  - Contrôles : |G₁|/|W| = 10 321 920, soit ×2 avec l'échange des voisins = 20,6 M, comme la table.
  - Forme des y (impairs sur Kᶜ, pairs sur K, W(a) ∈ C, Σ_{Kᶜ} a ≡ 4 mod 8) : à Γ près, un y = un motif de valeurs absolues. Classification des paires à faire.

## Pistes pour conclure vite (à discuter)
1. **Transférer la preuve de C3 (ordre 3, résolu « par les maths »).** Quel en est le mécanisme ? Est-ce qu'il s'applique à l'involution (un point fixe, 7 blocs, spectre 27/22) ?
2. **Certificat en forme close.** Dans O₂₄, le dual optimal Y vaut EXACTEMENT 1 sur TOUS les u-candidats, avec la même borne 38 < 42 dans toutes les classes testées.
   - Une forme Y explicite (Golay + positions des ±3 des y) donnerait une preuve à la main des 86 classes.
   - Son analogue « racines d'ombre + y » pourrait se généraliser aux autres réseaux.
3. **Réduction à une matrice W = 2S + I** : symétrique, entière, W² = 49I, diagonale (1⁴², −1⁷), hors diagonale dans {0,±2}, 12 non nuls par ligne.
   - Autrement dit, une « matrice orthogonale rationnelle de dénominateur 7 », avec la structure de blocs du graphe.
   - Chercher une obstruction de type Hasse-Minkowski / code 7-aire autodual (lignes de W mod 7 : code auto-orthogonal de dim a ≤ 24, mots de poids 13).
4. **Ne plus énumérer les paires : ne fixer que y₁ à Γ près, et traiter y₂ dans le PL.**
   - Y₂ = y₂y₂ᵀ, PSD de rang 1, trace 7, Y₂y₁ = 0, avec les valeurs exactes en racines, racines d'ombre et normes 3.
   - Sur O₂₄ cela ne peut pas tout tuer (55 classes PL-faisables), mais sur les réseaux riches en racines, les contraintes en racines pourraient suffire.
5. **Symétrie partout.** Pour chaque famille de voisins (A₁²⁴, A₂¹², A₃⁸, A₄⁶, D₄⁶…), Aut(N0) = W ⋊ Aut(code) est explicite, donc Γ(U) = Stab([v]).
   - Classer les paires à Γ près, comme les 141 de O₂₄.
