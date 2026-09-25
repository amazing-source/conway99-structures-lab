# Vérification indépendante du certificat d'Astra (mission géométrique, reçu le 25/09)

*Labo « structures ». Code écrit de zéro, en entiers exacts : `scripts/verify_astra_certificate.py`.*

## Ce qui est vérifié
1. **Algèbre.**
   - N^T N = 11I + J − Π, vérifié sur les 84 sommets.
   - C_yy = 20.
   - C N = 0, refait à la main à partir de HN = 2J − N(I + Π). Attention : 1_{84}ᵀN = 12·1ᵀ (sommes de colonnes), pas
     2·1ᵀ (sommes de lignes). Avec la bonne valeur, le terme en J s'annule : −16 − 6 − 2 + 24 = 0.
   - Spectre de H sur Im N : 12, 0⁷, (−2)⁶. Sur l'orthogonal : 3⁴⁰, (−4)³⁰. Donc C = 56·P₋₄, de rang 30.
2. **Certificat sur 23 sommets.**
   - Branche : réflexion Q du témoin Z₇ (b ↔ 14 − b + 2(b mod 2)) dans les étoiles de 0 et de 1, g = (2 3 … 13), et
     N(y₂,₄) ∩ (S₀ ∪ S₁) = {y₀,₂, y₀,₅, y₁,₆, y₁,₁₀}.
   - Résultat : **wᵀC₂₃w = −2**, avec Σ_{i<j} w_i w_j H_ij = 1568 et 25086 − 16·1568 = −2.
3. **Certificat sur le graphe (38 sommets : x, 14 étiquettes, 23 extérieurs).**
   - 1ᵀq = 0, ‖q‖² = 251 994, qᵀAq = 755 998, donc **qᵀAq − 3‖q‖² = 16 > 0**.
   - Or A ≤ 3 sur 1^⊥ (spectre 14, 3⁵⁴, (−4)⁴⁴). La branche est donc impossible, quelles que soient les arêtes
     inconnues : elles ne figurent pas dans cette forme.
4. **Cohérence de la branche.**
   - y₂,₄ a 2 voisins contenant 0 et 2 contenant 1 (règle HN).
   - L'arête y₂,₄ ~ y₀,₂ est celle du témoin en l'étiquette 2 (f₂(4) = 0).
   - Aucune violation de λ = 1 parmi les données fixées.

## Lecture
- C'est le **principe de Rayleigh** (entrelacement) appliqué à une branche partielle : un vecteur q ⊥ 1, porté par
  38 sommets entièrement connus, a un quotient > 3. Les propagations de comptage (degrés, λ, μ, quotas HN) ne le voient
  pas.
- **Les « 8 dimensions restantes » d'Astra sont exactement le projecteur caché Π₋₄ de rang 8 du X12 de C1**, pour la
  droite ℓ = {x, 0, 1}.
  - Les 60 sommets lointains de ℓ sont les étiquettes évitant 0 et 1.
  - Le complément orthogonal, dans l'espace −4 de H (dim 30), de l'espace engendré par les deux étoiles est
    {f ∈ E₋₄(A) : supp f ⊆ F}. En effet, un vecteur propre −4 porté par F a des sommes nulles aux étiquettes, donc il
    est dans K₋₄(x).
  - Donc S/56 = Π₋₄ de C1, exactement.
- **Portée (Astra le dit lui-même).** Une branche d'un témoin prescrit est exclue. Ce n'est ni une exclusion du témoin
  affine, ni une non-existence.
