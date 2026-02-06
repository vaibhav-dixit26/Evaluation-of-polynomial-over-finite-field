load("Poly_helper.sage")

#Gen evaluation
def EVA_Gen(arrays, x_vec, H_n, H_eps, d, k_list, q):
    """
    arrays : dict from prepare_all_memory
    x_vec  : (x_n, ..., x_1) in F_q
    H_n      : Tw_n(x)
    H_eps  : Tw_eps(x) with H_eps[0] = 0
    d      : degree
    k_list : (k_h, ..., k_1,0)
    q      : field size
    """
    count_ops=0
    h = len(k_list)
    # --------------------------------------------------
    # Case H_n = 1
    # --------------------------------------------------
    if H_n == 1:
        # print("H=1, only single step calculation")
        A0 = arrays[(0,)]
        
        k1 = k_list[-1]
        
        A0[k1 + 1] = (A0[1] + arrays[(k1,)][0]) % q
        # print(f"A{arrays[(0,)]}[k1+1]=A0[1] + A{arrays[{}]} ")
        count_ops+=1
        # print(A0, k1,arrays[(k1,)][0])
        return A0[k1 + 1], arrays, count_ops
    
    # --------------------------------------------------
    # Main loops
    # --------------------------------------------------
    for eps_idx in range(h, 0, -1):
        k_eps = k_list[h-eps_idx]
        η = H_eps[eps_idx - 1]
        x_keps = x_vec[len(x_vec) - k_eps]
        for δ in range(min(x_keps-1,d-1-η), -1, -1):
            step = η + δ
            # --------------------------------------------------
            # Case step = H - 1
            # --------------------------------------------------
            if step == H_n-1 :

                if x_keps > 1:
                    label = build_label(k_list, x_vec, eps_idx, delta=δ)
                    label_next = build_label(k_list, x_vec, eps_idx, delta=δ + 1)
                    # print(f"label1={label}, label2={label_next}")
                    arrays[label][1] = (arrays[label][0] + arrays[label_next][0]) % q
                    count_ops+=1
                    
                else:
                    label1 = build_label(k_list, x_vec, eps_idx, delta=0 )
                    label2 = build_label(k_list, x_vec, eps_idx, delta=1)
                    num1 = 0 if eps_idx==1 else 1
                    # print(f"label1={label1}, label2={label2}")
                    idx = k_eps - k_list[(h-eps_idx+1)*num1]*num1 + 1

                    arrays[label1][idx] = (arrays[label1][0] + arrays[label2][0]) % q
                    count_ops+=1
            # --------------------------------------------------
            # Case step <= H_n - 2
            # --------------------------------------------------
            elif step <= H_n - 2:
                num = 0 if step == d - 1 else 1
                # print("going to multiply by num=", num)
                # ------------------------------
                # δ > 0
                # ------------------------------
                if δ > 0:
                    a = x_keps - δ
                    label = build_label(k_list, x_vec,  eps_idx, delta=δ)
                    label_next = build_label(k_list, x_vec,  eps_idx, delta=δ + 1)
                    # print(f"label1={label}, label2={label_next}")
                    if a > 1:
                        arrays[label][1] = (arrays[label][1] + arrays[label_next][1 * num]) % q                        
                        count_ops+=1
                    else:
                        idx = k_list[h-eps_idx - 1] - k_eps + 1
                        arrays[label][1] = (arrays[label][idx]+ arrays[label_next][idx * num]) % q
                        count_ops+=1
                # ------------------------------
                # δ = 0
                # ------------------------------
                else:
                    a = x_keps
                    label1 = build_label(k_list, x_vec, eps_idx, delta=0 )
                    label2 = build_label(k_list, x_vec, eps_idx, delta=1)
                    # print(f"label1={label1}, label2={label2}")
                    num1 = 0 if eps_idx==1 else 1
                    idx = k_eps - k_list[(h-eps_idx+1)*num1]*num1 + 1
                    if a > 1:
                        arrays[label1][idx] = (arrays[label1][idx]+ arrays[label2][1 * num]) % q
                        count_ops+=1
                    else:
                        num1 = 0 if eps_idx==1 else 1
                        idx1 = k_list[h-eps_idx-1]  - k_list[(h-eps_idx+1)*num1]*num1  + 1
                        arrays[label1][idx] = (arrays[label1][idx1]+ arrays[label2][(idx1-idx+1) * num]) % q
                        count_ops+=1
    # --------------------------------------------------
    # Return
    # --------------------------------------------------
    return arrays[(0,)][k_list[-1] + 1] % q, arrays,count_ops

def ENU_Gen(K, H, I, W, N, h, end, q, 
current_vector,H_eps, newarrays
                            ):
    t_oprs=0
    for i in range(1, end + 1):
        num = I[i]
        
        for repeat in range(1, q):
            
            K[h + 1] = i
            current_vector[len(current_vector)-i] = repeat
            H[num] += 1
            H_eps = compute_H_eps(K, h, current_vector)
            
            val, newarrays, oprs = EVA_Gen(
                arrays=newarrays,
                x_vec=current_vector,
                H_n=H_eps[-1],
                H_eps=H_eps,
                d=d,
                k_list=K[1:h+2],
                q=q
                )
            # print("curr_vector and val",current_vector, val)
            our_eva_res[tuple(current_vector)] = val
            t_oprs+=oprs
            
            # print(f"{current_vector} {K[1:h+2]} {H_eps[:h+2]}")
            
            if H[num] < W[num]:
                _, child_oprs = ENU_Gen(
                    K, H, I, W, N,
                    h + 1, K[h + 1] - 1,
                    q, current_vector, H_eps, newarrays
                )
                t_oprs += child_oprs
            elif H[num] == W[num] and num > 0:
                _, child_oprs = ENU_Gen(
                    K, H, I, W, N,
                    h + 1, N[num - 1],
                    q, current_vector, H_eps, newarrays
                )
                t_oprs += child_oprs 
            # Backtrack: Reset vector and counter
            current_vector[len(current_vector)-i] = 0
            H[num] -= 1
    return our_eva_res, t_oprs





n =8
q = 3
degree = 3

'''
If want to fix N and W then 
write N = [n1, n2,...] such that sum(n_i)==n
and   W = [w1, w2,...] such that 0<w_i<=n_i
'''
N, W = random_W_N(n)                               # generate random ni and wi
poly = gen_random_poly(n, q, degree, sparsity=0.5) # generate random polynomial
f, d = build_polynomial(poly, n, q)                # build a function same as generated polynomial for direct evaluation.
N_org=N[:]
W_org=W[:]
print("Variables(n):", n, ", Degree(d):", d, ", Order of field(q):",q, "\n","polynomial:", poly)
print(f"|F_q^n| (full space size) : {q**n}")
print("\nStructured partial space parameters:")
print("  Blocks (N_i, W_i):")
print("N:", N)
print("W:", W)

arrays, _ = prepare_all_memory(n, d)
derivatives_at_0 = compute_all_derivatives_at_0(arrays, f, n, q)

#Initialization Phase: Initiate the first entry of each array
for label, value in derivatives_at_0.items():
    if label not in arrays:
        raise KeyError(f"Label {label} not found in arrays")

    if label == (0,):
        arrays[label][1] = value   # special rule for (0,)
    else:
        arrays[label][0] = value

our_eva_res = {}
newarrays=arrays

length = len(N)


for i in range(1, length):
    N[i] = N[i] + N[i - 1]

n = N[length - 1]
I = [0] * (n + 1)
for i in range(1, n + 1):
    for j in range(1, length):
        if N[j-1] < i <= N[j]:
            I[i] = j

K = [0] * (n + 1)
K[0] = n + 1
H_eps=[0]*(n+1)
H = [0] * length
initial_vector = [0] * n

our_eva_res, t_oprs =ENU_Gen(K, H, I, W, N, 0, K[0]-1, q, initial_vector,H_eps, newarrays)
our_eva_res[tuple(initial_vector)] = arrays[(0,)][1]

print("Evaluation part done")

print("\nFunctional Correctness Check:")
exhaustive_eva_res = {x: f(x) for x in product(range(q), repeat=n)}
card_part_space=0
for x in product(range(q), repeat=n):
    if not in_structured_space(list(x), N_org, W_org):
        continue
    
    card_part_space+=1
    eva_val  = our_eva_res[x]
    poly_val = exhaustive_eva_res[x]

    if eva_val != poly_val:
        raise ValueError(
            f"Mismatch detected!\n"
            f"  x        = {x}\n"
            f"  EVAGen(x)   = {eva_val}\n"
            f"  f(x)     = {poly_val}\n"
        )
else:
    print("All values match")

print("\nTime Complexity Check:")
if (t_oprs>d*card_part_space):
    raise ValueError(
            f"time error\n"
        )
else:
    print(f"|S| (partial space size)               : {card_part_space}")
    print(f"Total '+' operations in the eval phase : {t_oprs}")
    print(f"Theoretical bound (d · |S|)            : {d * card_part_space}")

print("\nMemory Complexity Check:")
alloted_mem = sum(len(v) for v in arrays.values())
if (alloted_mem-1 !=2 * comb(n+d, d) - 1):
    raise ValueError(
            f"time error\n"
    )
else:
    
    print(f"  Allocated memory  : arrays A_(.)    ={alloted_mem - 1}")  # alloted_mem-1 because only A[(0,)][1:] is used. That is A[(0,)][0] is nowhere used.
    print(f"  Theoretical bound : 2·C(n+d, d) − 1 ={2 * comb(n+d, d) - 1},\n\n\n")


    
