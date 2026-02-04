load("Poly_helper.sage")
#Quadratic evaluation
def EVAQuadratic(arrays, x_vec, H_n, H_eps, d, k_list, q):
    count_ops=0
    l = len(k_list)
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
            print(f"x: {current_vector}, f(x): {val}" )
            if val!=f(current_vector):
                
                raise ValueError(
                print("wrong result")
                
                )
            # print(f"{current_vector} {K[1:h+2]},{H_eps[:h+2]}, {H_eps[:]}")
            
            if H[num] < W[num]:
                ENUQuadratic(
                    K, H, I, W, N,
                    h + 1, K[h + 1] - 1,
                    q, current_vector, H_eps, newarrays
                )
            elif H[num] == W[num] and num > 0:
                ENUQuadratic(
                    K, H, I, W, N,
                    h + 1, N[num - 1],
                    q, current_vector, H_eps, newarrays
                )
            # Backtrack: Reset vector and counter
            current_vector[len(current_vector)-i] = 0
            H[num] -= 1
    



n = 5
q = 3
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
print(f"x: {initial_vector}, f(x): {arrays[(0,)][1]}" ) #print first vector=0^n and the f(first vector)= constant part of f.

ENUQuadratic(K, H, I, W, N, 0, K[0]-1, q, initial_vector,H_eps, newarrays)

print("Evaluation part done")

    