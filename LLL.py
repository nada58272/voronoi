import numpy as np

def gram_schmidt(basis):
    n = len(basis)
    b_star = [np.array(basis[i], dtype=float) for i in range(n)]
    mu = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i):
            mu[i,j] = np.dot(basis[i], b_star[j]) / np.dot(b_star[j], b_star[j])
            b_star[i] -= mu[i,j] * b_star[j]
        mu[i,i] = np.dot(b_star[i], b_star[i])  # β_i = |b*_i|²
    
    return b_star, mu


def size_reduction(basis, b_star, mu, k, l):
    """Выполняет Size-Reduction шаг."""
    q = round(mu[k, l])
    basis[k] -= q * basis[l].copy()
    mu[k, l] -= q
    for j in range(l):
        mu[k, j] -= q * mu[l, j]


def interchange(basis, mu, b_star, beta, k):
    """Выполняет Interchange шаг."""
    # Меняем местами векторы
    basis[k], basis[k - 1] = basis[k - 1].copy(), basis[k].copy()
    
    # Обновляем mu
    for j in range(k - 1):
        #print('mu[k, j]', mu[k, j], mu[k-1, j])
        mu[k, j], mu[k - 1, j] = mu[k - 1, j], mu[k, j]
        #print('mu[k, j]', mu[k, j])
    
    mu_k_km1 = mu[k, k - 1]
    beta_k = beta[k] - mu_k_km1 ** 2 * beta[k - 1]
    beta_km1 = beta[k - 1]
    
    
    # Обновляем mu
    mu[k, k - 1] = mu_k_km1 * beta_km1 / beta_k
    b = b_star[k-1]
    b_star[k-1] = b_star[k] + mu_k_km1 * b
    b_star[k] = (beta[k]/beta_k) * b - mu_k_km1 * b_star[k]
    
    # Обновляем beta
    beta[k] = beta[k - 1] * beta[k] / beta_k
    beta[k - 1] = beta_k
    
    
    for i in range(k + 1, len(basis)):
        temp = mu[i, k]
        mu[i, k] = mu[i, k - 1] - mu_k_km1 * temp
        mu[i, k - 1] = temp + mu[k, k-1] * mu[i, k]


def LLL(basis, delta=0.75):
    n = len(basis)
    basis = [np.array(vec, dtype=float) for vec in basis]
    b_star, mu = gram_schmidt(basis)
    beta = [np.dot(b_star[i], b_star[i]) for i in range(n)]
    
    k = 1
    while k < n:
        # Size Reduction
        for l in range(k-1, -1, -1):
            if abs(mu[k,l]) > 0.5:
                q = round(mu[k,l])
                basis[k] -= q * basis[l]
                # Обновляем mu
                mu[k,l] -= q
                for j in range(l):
                    mu[k,j] -= q * mu[l,j]
                # Пересчитываем Грама-Шмидта после изменения базиса
                b_star, mu = gram_schmidt(basis)
                beta = [np.dot(b_star[i], b_star[i]) for i in range(n)]
        
        # Проверка условия Ловаса
        if beta[k] < (delta - mu[k,k-1]**2) * beta[k-1]:
            # Interchange
            basis[k], basis[k-1] = basis[k-1].copy(), basis[k].copy()
            # Пересчитываем всё после обмена
            b_star, mu = gram_schmidt(basis)
            beta = [np.dot(b_star[i], b_star[i]) for i in range(n)]
            k = max(1, k-1)
        else:
            k += 1
    
    return np.array(basis)