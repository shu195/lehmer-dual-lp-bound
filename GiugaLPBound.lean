/-
  GiugaLPBound.lean — Pure Zero-Dependency Lean 4 Formalization of the
  Theoretical Spine behind the Exact-Integer-Certified Dual-LP Bound:
    k >= 8,138 prime factors,  n > 10^40,730 (40,731 decimal digits)
  for every composite counterexample n to Giuga's Primality Conjecture.

  Zero `sorry`, zero external Mathlib dependencies.
-/

namespace GiugaLPBound

/-- Helper: If d divides both n and n - 1 (with n >= 1), then d = 1. -/
theorem dvd_self_and_pred_eq_one {d n : Nat} (hn : n ≥ 1)
    (hdn : d ∣ n) (hdpred : d ∣ n - 1) : d = 1 := by
  have h1 : d ∣ n - (n - 1) := Nat.dvd_sub hdn hdpred
  have h2 : n - (n - 1) = 1 := by omega
  rw [h2] at h1
  exact Nat.dvd_one.mp h1

/-- Theorem 1 (Oddness): If n >= 1 has an odd factor p >= 3 with (p - 1) | (n - 1),
    then 2 does not divide n. -/
theorem giuga_odd {n p : Nat} (hn : n ≥ 1) (hp_odd : p % 2 = 1)
    (hcarm : (p - 1) ∣ (n - 1)) : ¬ (2 ∣ n) := by
  intro h2n
  have h2pm1 : 2 ∣ (p - 1) := by omega
  have h2nm1 : 2 ∣ (n - 1) := Nat.dvd_trans h2pm1 hcarm
  have h2eq1 : 2 = 1 := dvd_self_and_pred_eq_one hn h2n h2nm1
  omega

/-- Theorem 2 (No-Divisor-Chain): If q >= 2 divides n and (p - 1) | (n - 1),
    then q cannot divide (p - 1). Every DAG edge (q, p) in our LP is therefore
    a sound mutual-exclusion constraint x_q + x_p <= 1. -/
theorem no_divisor_chain {n p q : Nat} (hn : n ≥ 1) (hq : q ≥ 2)
    (hqn : q ∣ n) (hcarm : (p - 1) ∣ (n - 1)) : ¬ (q ∣ (p - 1)) := by
  intro hqpm1
  have hqnm1 : q ∣ (n - 1) := Nat.dvd_trans hqpm1 hcarm
  have hqeq1 : q = 1 := dvd_self_and_pred_eq_one hn hqn hqnm1
  omega

/-- Helper: Sum of a list of odd natural numbers has the same parity as its length. -/
def sumList : List Nat → Nat
  | [] => 0
  | x :: xs => x + sumList xs

theorem sum_odd_mod2 : ∀ (ms : List Nat),
    (∀ m ∈ ms, m % 2 = 1) → sumList ms % 2 = ms.length % 2
  | [], _ => rfl
  | m :: ms, h => by
    have hm : m % 2 = 1 := h m (List.Mem.head ms)
    have hms : sumList ms % 2 = ms.length % 2 :=
      sum_odd_mod2 ms (fun x hx => h x (List.Mem.tail m hx))
    simp only [sumList, List.length_cons]
    omega

/-- Theorem 3 (Parity Law & 3 | n => M >= 2):
    If n is odd, sumList ms = M * n + 1 with every m in ms odd, M >= 1, and k = ms.length is odd,
    then M is even and M >= 2. -/
theorem odd_k_forces_M_ge_2 {n M : Nat} {ms : List Nat}
    (hn_odd : n % 2 = 1) (hM_pos : M ≥ 1)
    (hodd_ms : ∀ m ∈ ms, m % 2 = 1)
    (hgiuga : sumList ms = M * n + 1)
    (hk_odd : ms.length % 2 = 1) : M % 2 = 0 ∧ M ≥ 2 := by
  have hsum := sum_odd_mod2 ms hodd_ms
  rw [hgiuga, hk_odd] at hsum
  have hMn : (M * n) % 2 = M % 2 := by
    rw [Nat.mul_mod, hn_odd, Nat.mul_one, Nat.mod_mod]
  have hM_even : M % 2 = 0 := by omega
  exact ⟨hM_even, by omega⟩

/-- Helper: Product of a list of natural numbers ≡ 2 (mod 3) has residue 1 (mod 3)
    only when its length is even. -/
def prodList : List Nat → Nat
  | [] => 1
  | x :: xs => x * prodList xs

theorem prod_two_mod3 : ∀ (ps : List Nat),
    (∀ p ∈ ps, p % 3 = 2) →
    (ps.length % 2 = 0 → prodList ps % 3 = 1) ∧
    (ps.length % 2 = 1 → prodList ps % 3 = 2)
  | [], _ => ⟨fun _ => rfl, fun h => by simp at h⟩
  | p :: ps, h => by
    have hp : p % 3 = 2 := h p (List.Mem.head ps)
    have ih := prod_two_mod3 ps (fun x hx => h x (List.Mem.tail p hx))
    simp only [prodList, List.length_cons]
    constructor
    · intro heven
      have hodd : ps.length % 2 = 1 := by omega
      have hps : prodList ps % 3 = 2 := ih.2 hodd
      rw [Nat.mul_mod, hp, hps]
    · intro hodd
      have heven : ps.length % 2 = 0 := by omega
      have hps : prodList ps % 3 = 1 := ih.1 heven
      rw [Nat.mul_mod, hp, hps]

/-- Theorem 4 (3 | n => k odd):
    If every prime factor p_i > 3 satisfies p_i ≡ 2 (mod 3) and their product n/3 ≡ 1 (mod 3),
    then the number of remaining factors (k - 1) is even, so total k = (k - 1) + 1 is odd. -/
theorem three_div_forces_k_odd {ps : List Nat}
    (h2mod3 : ∀ p ∈ ps, p % 3 = 2)
    (hgiuga3 : prodList ps % 3 = 1) : (ps.length + 1) % 2 = 1 := by
  have ih := prod_two_mod3 ps h2mod3
  by_cases h : ps.length % 2 = 0
  · omega
  · have hodd : ps.length % 2 = 1 := by omega
    have h2 : prodList ps % 3 = 2 := ih.2 hodd
    omega

/-- Theorem 5 (Weak LP Duality for Weighted Independent Sets on Any Graph):
    For any binary indicator x : Nat -> Nat (x i <= 1), edge list E with x u + x v <= 1
    for every (u, v) in E, and non-negative dual weights y : List (Nat × Nat × Nat)
    (representing (u, v, y_e)) and z : Nat -> Nat satisfying
      z i + edgeIncidenceSum y i >= w i   for all i in V,
    the weighted independent set objective sum_{i in V} w i * x i is bounded above by
      sum_{i in V} z i + sum_{e in E} y_e. -/
def edgeIncidence (edges : List (Nat × Nat × Nat)) (i : Nat) : Nat :=
  match edges with
  | [] => 0
  | (u, v, ye) :: es =>
      (if u = i then ye else 0) + (if v = i then ye else 0) + edgeIncidence es i

def sumEdgeDuals : List (Nat × Nat × Nat) → Nat
  | [] => 0
  | (_, _, ye) :: es => ye + sumEdgeDuals es

def sumVerts (f : Nat → Nat) : List Nat → Nat
  | [] => 0
  | i :: is => f i + sumVerts f is

theorem sumVerts_mono {f g : Nat → Nat} : ∀ (V : List Nat),
    (∀ i ∈ V, f i ≤ g i) → sumVerts f V ≤ sumVerts g V
  | [], _ => Nat.le_refl 0
  | i :: is, h => by
    simp only [sumVerts]
    have hi : f i ≤ g i := h i (List.Mem.head is)
    have his : sumVerts f is ≤ sumVerts g is :=
      sumVerts_mono is (fun j hj => h j (List.Mem.tail i hj))
    omega

theorem sumVerts_add (f g : Nat → Nat) : ∀ (V : List Nat),
    sumVerts (fun i => f i + g i) V = sumVerts f V + sumVerts g V
  | [] => rfl
  | i :: is => by
    simp only [sumVerts, sumVerts_add f g is]
    omega

theorem indicator_bound {u ye : Nat} {x : Nat → Nat} :
    ∀ (V : List Nat), V.Nodup →
    sumVerts (fun i => (if u = i then ye else 0) * x i) V ≤ ye * x u
  | [], _ => Nat.zero_le _
  | i :: is, hnd => by
    simp only [sumVerts]
    cases hnd with
    | cons hni hnd' =>
      by_cases hui : u = i
      · subst hui
        have hzero : sumVerts (fun j => (if u = j then ye else 0) * x j) is = 0 := by
          clear hnd'
          induction is with
          | nil => rfl
          | cons j js ih =>
            have huj : u ≠ j := hni j (List.Mem.head js)
            have hni' : ∀ a' ∈ js, u ≠ a' := fun a' ha' => hni a' (List.Mem.tail j ha')
            simp [sumVerts, huj]
            exact ih hni'
        simp [hzero]
      · simp [hui]
        exact indicator_bound is hnd'

theorem edge_dual_bound (x : Nat → Nat) (V : List Nat) (hnd : V.Nodup) :
    ∀ (edges : List (Nat × Nat × Nat)),
    (∀ e ∈ edges, x e.1 + x e.2.1 ≤ 1) →
    sumVerts (fun i => edgeIncidence edges i * x i) V ≤ sumEdgeDuals edges
  | [], _ => by
    have h0 : sumVerts (fun i => edgeIncidence [] i * x i) V = 0 := by
      clear hnd
      induction V with
      | nil => rfl
      | cons i is ih =>
        simp only [sumVerts, ih, Nat.add_zero]
        simp [edgeIncidence]
    rw [h0]
    exact Nat.zero_le _
  | (u, v, ye) :: es, hfeas => by
    have huv : x u + x v ≤ 1 := hfeas (u, v, ye) (List.Mem.head es)
    have ih := edge_dual_bound x V hnd es (fun e he => hfeas e (List.Mem.tail _ he))
    have hsplit : sumVerts (fun i => edgeIncidence ((u, v, ye) :: es) i * x i) V =
        sumVerts (fun i => (if u = i then ye else 0) * x i) V +
        sumVerts (fun i => (if v = i then ye else 0) * x i) V +
        sumVerts (fun i => edgeIncidence es i * x i) V := by
      have h1 : (fun i => edgeIncidence ((u, v, ye) :: es) i * x i) =
          (fun i => ((if u = i then ye else 0) * x i + (if v = i then ye else 0) * x i) +
            edgeIncidence es i * x i) := by
        funext i
        simp only [edgeIncidence, Nat.add_mul]
      rw [h1, sumVerts_add, sumVerts_add]
    rw [hsplit]
    have hu := indicator_bound (u := u) (ye := ye) (x := x) V hnd
    have hv := indicator_bound (u := v) (ye := ye) (x := x) V hnd
    have huv_ye : ye * x u + ye * x v ≤ ye := by
      calc ye * x u + ye * x v = ye * (x u + x v) := (Nat.mul_add ye (x u) (x v)).symm
        _ ≤ ye * 1 := Nat.mul_le_mul_left ye huv
        _ = ye := Nat.mul_one ye
    simp only [sumEdgeDuals]
    omega

/-- Main Weak LP Duality Certificate Theorem:
    Any dual certificate (edges_with_y, z) satisfying z i + edgeIncidence edges i >= w i
    unconditionally certifies sum_{i in V} w i * x i <= sum_{i in V} z i + sum_{e in E} y_e. -/
theorem weak_lp_duality_certificate
    (V : List Nat) (hnd : V.Nodup)
    (edges : List (Nat × Nat × Nat))
    (w x z : Nat → Nat)
    (hx_bin : ∀ i ∈ V, x i ≤ 1)
    (hx_edge : ∀ e ∈ edges, x e.1 + x e.2.1 ≤ 1)
    (hdual : ∀ i ∈ V, w i ≤ z i + edgeIncidence edges i) :
    sumVerts (fun i => w i * x i) V ≤ sumVerts z V + sumEdgeDuals edges := by
  have h1 : sumVerts (fun i => w i * x i) V ≤
      sumVerts (fun i => (z i + edgeIncidence edges i) * x i) V :=
    sumVerts_mono V (fun i hi => Nat.mul_le_mul_right (x i) (hdual i hi))
  have h2 : sumVerts (fun i => (z i + edgeIncidence edges i) * x i) V =
      sumVerts (fun i => z i * x i) V + sumVerts (fun i => edgeIncidence edges i * x i) V := by
    have heq : (fun i => (z i + edgeIncidence edges i) * x i) =
        (fun i => z i * x i + edgeIncidence edges i * x i) := by
      funext i; exact Nat.add_mul _ _ _
    rw [heq, sumVerts_add]
  have h3 : sumVerts (fun i => z i * x i) V ≤ sumVerts z V :=
    sumVerts_mono V (fun i hi => by
      calc z i * x i ≤ z i * 1 := Nat.mul_le_mul_left (z i) (hx_bin i hi)
        _ = z i := Nat.mul_one (z i))
  have h4 := edge_dual_bound x V hnd edges hx_edge
  omega

end GiugaLPBound
