===================================
SAGBI and Gröbner Bases Detection
===================================

| This page features the source code and explanations for the paper's computational results: 
| Viktoriia Borovik, Timothy Duff and Elima Shehu
| SAGBI and Gröbner Bases Detection
| ARXIV: https://arxiv.org/abs/
| CODE: https://mathrepo.mis.mpg.de/SAGBI_GB_detection/index.html

**Abstract of the paper.** We introduce a *SAGBI detection* algorithm, analogous to the *Gröbner basis detection* algorithm of Gritzmann and Sturmfels. We also present the accompanying package `SagbiGbDetection` implemented in both  ``Macaulay2`` and  ``Julia``. The package provides functions to find a term order :math:`\omega` such that the given generators of an ideal form a Gröbner basis or the given generators of a finitely generated subalgebra of a polynomial ring are the SAGBI basis with respect to :math:`\omega`. In addition, we investigate the complexity and provide several benchmarks that showcase the practical application of this package.

In the following we present the codes for the examples from the applications section. Our code is written either in the programming language `julia <https://julialang.org/>`_ or `macaulay2 <https://macaulay2.com/>`_. 

Cyclic Gaussian graphical model
---------------------------------------
We provide the code for the example on the Sullivant-Talaska ideal :math:`I_5`. First, we have to install all needed packages, including SagbiGbDetection then load our *pkg* and define :math:`I_5`.

.. code-block:: Julia

 using Singular, SagbiGbDetection
 using Combinatorics


 R, v = Singular.polynomial_ring(Singular.QQ, ["s_11", "s_14", "s_15", "s_21", "s_23", "s_24", "s_25", "s_31", "s_33", "s_34", "s_35"]);
  (s_11, s_14, s_15, s_21, s_23, s_24, s_25, s_31, s_33, s_34, s_35) = v;

  A = [s_11 s_31 s_14 s_15;
       s_21 s_23 s_24 s_25;
       s_31 s_33 s_34 s_35];

function define_ring(n)
    # Generate variables for the polynomial ring
     variables_R = vec(["s_$(i)$(j)" for i in 1:n, j in 1:n])

    # Initialize the polynomial ring
    return R, _ = Singular.polynomial_ring(Singular.QQ, variables_R)

end


function symmetric_generic_matrix(v, n)
    # Reshape the variables into an n x n matrix
    generic_matrix = reshape(v, n, n)

    # Make the matrix symmetric
    for i in 1:n
        for j in 1:i-1
            generic_matrix[i, j] = generic_matrix[j, i]
        end
    end

    return generic_matrix
end

  # Define a function to compute the determinant using Singular's det function
   function determinant_spoly(matrix::Matrix{spoly{n_Q}}) where {n_Q}
      # Convert the matrix elements to polynomials
      det_poly = Singular.det(Singular.matrix(matrix))
    
      return det_poly
   end

  function compute_minors(matrix::Matrix{T}, s::Int) where T
      m, n = size(matrix, 1), size(matrix, 2)
      if m < s || n < s || s < 1
          throw(ArgumentError("Invalid values for m, n, or s"))
      end

      determinants = []

      for rows_combination in combinations(1:m, s)
          for cols_combination in combinations(1:n, s)
              minor = matrix[rows_combination, cols_combination]
              push!(determinants, determinant_spoly(copy(minor)))
          end
      end

      return determinants
  end

function k_submatrices_minors(generic_matrix, rows, cols, k)
    minors_vec = []
    for i in 1:k
        sub_matrix = generic_matrix[rows[i], cols[i]]
        minors = compute_minors(sub_matrix, 3)
        append!(minors_vec, minors)
    end
    return minors_vec
end


n = 4 #cycle
# Define the ring
R, v = define_ring(n)

generic_matrix = symmetric_generic_matrix(v, n)

rows = Vector{Vector{Int}}(undef, n)
  for i in 1:n
     rows[i] = [(j + i - 2) % n + 1 for j in 1:3]
  end


cols = Vector{Vector{Int}}(undef, n)
  for i in 1:n
     cols[i] = [(j + i + 1) % n + 1 for j in 0:2]
  end 
  
G = k_submatrices_minors(generic_matrix, rows, cols, n)

weight_vectors = weightVectorsRealizingGB(G,R)


(Vector{ZZRingElem}[[8, 1, 1, 1, 20, 15, 1, 1, 14, 5, 21, 1, 20, 24, 5, 15], [15, 1, 1, 1, 5, 21, 1, 1, 24, 5, 15, 1, 20, 14, 20, 8], [20, 1, 1, 1, 8, 27, 1, 1, 12, 8, 20, 1, 16, 18, 16, 9], [21, 1, 1, 1, 5, 15, 1, 1, 14, 20, 8, 1, 5, 24, 20, 15], [9, 1, 1, 1, 16, 20, 1, 1, 18, 8, 27, 1, 16, 12, 8, 20], [15, 1, 1, 1, 20, 8, 1, 1, 24, 20, 15, 1, 5, 14, 5, 21], [26, 1, 1, 1, 12, 26, 1, 1, 22, 12, 26, 1, 12, 22, 12, 26], [20, 1, 1, 1, 16, 9, 1, 1, 12, 16, 20, 1, 8, 18, 8, 27], [27, 1, 1, 1, 8, 20, 1, 1, 18, 16, 9, 1, 8, 12, 16, 20]], false)

length(weight_vectors[1])

 9

Grassmannian :math:`Gr(2, 4)`
--------------------------------
We provide the code for the example on :math:`(2\times 2)`-minors of a general~ :math:`(3\times 3)`-matrix. Our code is written in the programming language ``Julia``. 

.. code-block:: Julia
 
 T, v = Singular.polynomial_ring(Singular.QQ, ["t11", "t12", "t13", "t21", "t22", "t23", "t31", "t32", "t33"]);
 (t11, t12, t13, t21, t22, t23, t31, t32, t33) = v;

  t_ij = [t11 t12 t13;
          t21 t22 t23;
          t31 t32 t33];
 
 minors = compute_minors(t_ij, 2);
 weightVectorsRealizingSAGBI(minors,T)

 (Vector{ZZRingElem}[[4, 4, 9, 4, 9, 4, 9, 4, 4], [4, 9, 4, 4, 4, 9, 9, 4, 4], [4, 4, 9, 9, 4, 4, 4, 9, 4], [9, 4, 4, 4, 4, 9, 4, 9, 4], [4, 9, 4, 9, 4, 4, 4, 4, 9], [9, 4, 4, 4, 9, 4, 4, 4, 9]], false)

    

Truncation variety :math:`V_{\{ 1, 3\}}`
-------------------------------------------
We present the code for the example on truncation variety in the programming language ``Julia``. 

.. code-block:: Julia

 S, v = Singular.polynomial_ring(Singular.QQ, ["s", "z1", "z2", "z3", "z4", "z5", "z6", "z7", "z8", "z9", "z10"]);
 (s, z1, z2, z3, z4, z5, z6, z7, z8, z9, z10) = v;

 Q = [s, s*z1, s*z2, s*z3, s*z4, s*z5, s*z6, s*z7, s*z8, s*z9, 
      s*(z1*z5 - z2*z4), s*(z1*z6 - z3*z4), s*(z2*z6 - z3*z5),
      s*(z1*z8 - z2*z7), s*(z1*z9 - z3*z7), s*(z2*z9 - z3*z8),
      s*(z4*z8 - z5*z7), s*(z4*z9 - z6*z7), s*(z5*z9 - z6*z8),
      s*(z10 + z1*(z5*z9 - z6*z8) - z2*(z4*z9 - z6*z7) + z3*(z4*z8 - z5*z7))];

 weightVectorsRealizingSAGBI(Q, S)

 Vector{fmpz}[]

Coupled cluster equations (CC equations) on a Grassmannian :math:`V_{1} Gr(2, 5)` in `\mathbb P^9`
---------------------------------------------------------------------------------------------------
We present the code for the example on Coupled cluster equations (CC equations) on a Grassmannian :math:`V_{1} Gr(2, 5)` in `\mathbb P^9`.

.. code-block:: Julia

 S, v = Singular.polynomial_ring(Singular.QQ, ["t", "l", "x13", "x14", "x15", "x23", "x24", "x25"]);
 (t, l, x13, x14, x15, x23, x24, x25) = v;

 Q = [t, t*x13, t*x14, t*x15, t*x23, t*x24, t*x25, 
      t*(x14*x23 - x13*x24), t*(x15*x23 - x13*x25), t*(x15*x24 - x14*x25),
      t*l, t*l*x13, t*l*x14, t*l*x15, t*l*x23, t*l*x24, t*l*x25];
     
 weightVectorsRealizingSAGBI(Q, S)
 (Vector{ZZRingElem}[[1, 1, 2, 3, 4, 4, 3, 2], [1, 1, 3, 2, 4, 3, 4, 2], [1, 1, 2, 4, 3, 4, 2, 3], [1, 1, 4, 2, 3, 2, 4, 3], [1,  1, 3, 4, 2, 3, 2, 4], [1, 1, 4, 3, 2, 2, 3, 4]], true)
 
Algebras Generated by Principal minors.
---------------------------------------------------------------------------------------------------
We present the code for the example on Algebras Generated by Principal minors.

.. code-block:: macaulay2

needsPackage "SagbiGbDetection"
R = QQ[t,a,b,c,d,e,f,MonomialOrder => Eliminate 1]
M = matrix{{a,b,c},{b,d,e},{c,e,f}}
lst = apply(subsets 3, S -> t*det(M_S^S))
S = subring lst

--principal minors are not sagbi w.r.t. any term order
elapsedTime goodWs = weightVectorsRealizingSAGBI lst

-- extract all equivalence classes w.r.t. principal minors
elapsedTime ws = extractWeightVectors lst

-- compute dim, degree and ehrhart polynomial for every weight
ehrharts = apply(ws, w -> (
	Rw := QQ[gens R, MonomialOrder => {Weights => w}];
        lst = apply(lst, f -> substitute(f, Rw));
        param = apply(lst, leadTerm);
        Q = QQ[z_1..z_8];
        phi = map(Rw, Q, param); I = kernel phi;
    	hilbertPolynomial(I, Projective=> false), dim (I) - 1, degree (I)
	)
    )
netList ehrharts

-- compute final table with degrees and number of generators for every weight
netList apply(ws, ehrharts, (w,e) -> {e#1, e#2, isSAGBI gens sagbi(sub(gens S, QQ[gens R, MonomialOrder => {Weights => w}]), Limit => 5),
        isSAGBI gens sagbi(sub(gens S, QQ[gens R, MonomialOrder => {Weights => w}]), Limit => 6),
        elapsedTime #(flatten entries gens sagbi(sub(gens S, QQ[gens R, MonomialOrder => {Weights => w}])))})

-- the union of some equivalence classes w.r.t. principal minors form equivalence class w.r.t. S

flatten entries gens sagbi(sub(gens S, QQ[gens R, MonomialOrder => {Weights => ws#6}])) 
flatten entries gens sagbi(sub(gens S, QQ[gens R, MonomialOrder => {Weights => ws#9}]))
flatten entries gens sagbi(sub(gens S, QQ[gens R, MonomialOrder => {Weights => ws#12}]))

flatten entries gens sagbi(sub(gens S, QQ[gens R, MonomialOrder => {Weights => ws#8}]))
flatten entries gens sagbi(sub(gens S, QQ[gens R, MonomialOrder => {Weights => ws#10}]))
flatten entries gens sagbi(sub(gens S, QQ[gens R, MonomialOrder => {Weights => ws#11}]))