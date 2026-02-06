load("Poly_helper.sage")

#Quadratic evaluation
def EVAQuadratic(arrays, x_vec, H_n, H_eps, d, k_list, q):
    count_ops=0
    if H_n == 1:
        
        k1 = k_list[-1]
        
        arrays[(0,)][k1 + 1] = (arrays[(0,)][1] + arrays[(k1,)][0]) % q
        count_ops+=1

    if H_n==2: 
        k1 = k_list[-1]
        x_k1 = x_vec[len(x_vec) - k1]
        if x_k1>1:
            arrays[(k1,)][1]=(arrays[(k1,)][0] + arrays[(k1,k1,)][0]) % q
            arrays[(0,)][k1+1]=(arrays[(0,)][k1+1] + arrays[(k1,)][1]) % q
            count_ops+=2
        else: 
            k2 = k_list[-2]
            arrays[(k1,)][k2-k1+1]=(arrays[(k1,)][0] + arrays[(k1,k2,)][0]) % q
            arrays[(0,)][k1+1]=(arrays[(0,)][k2+1] + arrays[(k1,)][k2-k1+1]) % q
            count_ops+=2
    if H_n>2: 
        k1 = k_list[-1]
        x_k1 = x_vec[len(x_vec) - k1]
        if x_k1>1:
            if x_k1>2:
                arrays[(k1,)][1]=(arrays[(k1,)][1] + arrays[(k1,k1,)][0]) % q
                count_ops+=1
            else: 
                k2 = k_list[-2]
                arrays[(k1,)][1]=(arrays[(k1,)][k2-k1+1] + arrays[(k1,k1,)][0]) % q
                count_ops+=1
            arrays[(0,)][k1+1]=(arrays[(0,)][k1+1] + arrays[(k1,)][1]) % q 
            count_ops+=1
        else: 
            k2 = k_list[-2]
            x_k2 = x_vec[len(x_vec) - k2]
            if x_k2> 1:
                arrays[(k1,)][k2-k1+1]=(arrays[(k1,)][k2-k1+1] + arrays[(k1,k2)][0]) % q
                count_ops+=1
            else:
                k3 = k_list[-3]
                arrays[(k1,)][k2-k1+1]=(arrays[(k1,)][k3-k1+1] + arrays[(k1,k2)][0]) % q
                count_ops+=1
            arrays[(0,)][k1+1]=(arrays[(0,)][k2+1] + arrays[(k1,)][k2-k1+1]) % q
            count_ops+=1
            
    
    return arrays[(0,)][k_list[-1] + 1] % q, arrays,count_ops




def ENUQuadratic(K, H, I, W, N, h, end, q, 
current_vector,H_eps, newarrays
                            ):
    # our_eva_res={}
    t_oprs=0
    for i in range(1, end + 1):
        num = I[i]
        
        for repeat in range(1, q):
            
            K[h + 1] = i
            current_vector[len(current_vector)-i] = repeat
            H[num] += 1
            H_eps = compute_H_eps(K, h, current_vector)
            
            val, newarrays, oprs = EVAQuadratic(
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
            
            # print(f"{current_vector} {K[1:h+2]},{H_eps[:h+2]}, {H_eps[:]}")
            
            if H[num] < W[num]:
                _, child_oprs = ENUQuadratic(
                    K, H, I, W, N,
                    h + 1, K[h + 1] - 1,
                    q, current_vector, H_eps, newarrays
                )
                t_oprs += child_oprs
            elif H[num] == W[num] and num > 0:
                _, child_oprs = ENUQuadratic(
                    K, H, I, W, N,
                    h + 1, N[num - 1],
                    q, current_vector, H_eps, newarrays
                )
                t_oprs += child_oprs 
            # Backtrack: Reset vector and counter
            current_vector[len(current_vector)-i] = 0
            H[num] -= 1
    
    return our_eva_res, t_oprs



n = 7
q = 5
degree = 2

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

our_eva_res, t_oprs =ENUQuadratic(K, H, I, W, N, 0, K[0]-1, q, initial_vector,H_eps, newarrays)
our_eva_res[tuple(initial_vector)] = arrays[(0,)][1]

print("Evaluation part done")

print("\nFunctional Correctness Check:")
exhaustive_eva_res = {x: f(x) for x in product(range(q), repeat=n)}
card_part_space=0
''' 
The part below can be ignored by inserting this code in the ENUQuadratic function, as done in the Poly_evaQuadratic code. 
if val!=f(current_vector):
               
    raise ValueError(
    print("wrong result")
    
    )
But then we need to calculate the cardinality of the structured space separately, to verify the complexity units.
Here is the code for calculating the cardinality of a structured space. 
def structured_space_cardinality(n_list, w_list, q):
    from math import comb
    card = 1
    for n_j, w_j in zip(n_list, w_list):
        card *= sum(comb(n_j, i) * (q - 1)**i for i in range(w_j + 1))
    return card
print(structured_space_cardinality(N_org,W_org,q))

'''
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
            f"  EVA(x)   = {eva_val}\n"
            f"  f(x)     = {poly_val}"
        )
else:
    print("All values match")
print("\nTime Complexity Check:")
if (t_oprs>d*card_part_space):
    raise ValueError(
            f"time error\n"
            f"Total '+' operations in the eval phase : {t_oprs}"
            f"Theoretical bound (d · |S|)            : {d * card_part_space}"
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


    
