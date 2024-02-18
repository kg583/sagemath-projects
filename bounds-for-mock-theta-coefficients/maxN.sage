# SageMath code used to calculate maxN(r,2;n) and generate Table 4
# Kevin Gomez -- 8/7/2020

# Rank of a partition p
def rank(p):
  return max(p) - len(p)

# N(r,t;n), the r (mod t) rank counting function
def N(r,t,n):
  return len(list(filter(lambda p: rank(p) % t == r, Partitions(n))))

# Multiplicative extension of N(r,t;n) for evaluation at partitions
def pN(r,t,p):
  return prod(N(r,t,part) for part in p)

# maxN(r,t;n), as defined in our paper, and the list of all partitions which achieve this maximum
def maxN(r,t,n):
  pNs = {tuple(p): pN(r,t,p) for p in Partitions(n)}
  m = max(pNs.values())
  return m, list(filter(lambda p: pNs[p] == m, pNs.keys()))

# Generates values displayed in Table 4
# Equivalent maximal partitions are reduced to all 2's
def table4(m=24):
  for n in range(1, m):
    even = maxN(0,2,n)
    odd = maxN(1,2,n)
    print(n, even[0], even[1], odd[0], list(filter(lambda p: p.count(4) == 0 and p.count(6) == 0, odd[1])))

if __name__ == "__main__":
  table4()
