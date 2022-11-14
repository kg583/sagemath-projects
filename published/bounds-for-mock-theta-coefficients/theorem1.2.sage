# SageMath code used to verify Theorem 1.2 for a finite number of small cases outside the range of our asymptotics
# Kevin Gomez -- 8/7/2020

def l(n):
  return pi * sqrt(24 * n - 1) / 6

# List of values of N(r,2;n) for 1 <= n <= 9999
# Requires https://oeis.org/A000025/b000025.txt and https://oeis.org/A000041/b000041.txt
def N(r):
  with open("b000025.txt") as file:
    alpha = list(map(lambda line: int(line.split()[1]), file.readlines()))

  with open("b000041.txt") as file:
    partition = list(map(lambda line: int(line.split()[1]), file.readlines()))
  
  Ne = [(p + a) / 2 for p, a in zip(partition, alpha)]
  No = [(p - a) / 2 for p, a in zip(partition, alpha)]
  
  return Ne if r == 0 else No

# Verification of small n in Lemma 6.1
# Requires https://oeis.org/A000025/b000025.txt and https://oeis.org/A000041/b000041.txt
def verify61(m=4543):
  Ne, No = N(0), N(1)
  constant = pi ** 2 * sqrt(3) / 36
  for n in range(2, m):
      lower = constant / l(n) ** 2 * (1 - 1/l(n)) * (1 - 1/sqrt(n)) * exp(l(n))
      upper = constant / l(n) ** 2 * (1 - 1/l(n)) * (1 + 1/sqrt(n)) * exp(l(n))
      print(lower, Ne[n], No[n], upper)
      if lower < min(Ne[n], No[n]) and max(Ne[n], No[n]) < upper:
          print("TRUE -- {}".format(n))
      else:
          print("FALSE -- {}".format(n))


# Verification of small cases in Theorem 1.2
# Requires https://oeis.org/A000025/b000025.txt and https://oeis.org/A000041/b000041.txt
def check12():
  Ne, No = N(0), N(1)

  # Constants C_a as displayed in Table 3
  C = {11: 2.20,
       12: 1.86,
       13: 1.62,
       14: 1.43,
       15: 1.27,
       16: 1.15,
       17: 1.05}

  # Note that (a, b) = (11, 11) fails for No (as our theorem only supports a,b >= 12 for No)
  for a in range(11, 18):
    b = a
    while b/a <= C[a]:
      if Ne[a] * Ne[b] > Ne[a + b] and No[a] * No[b] > No[a + b]:
        print("TRUE -- {}, {}".format(a, b))
      else:
        print("FALSE -- {}, {}".format(a, b))
      b += 1

if __name__ == "__main__":
    verify61()
    check12()
