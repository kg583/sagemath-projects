i = j = 0
a, b = 8, 0

while True:
    if j % 100000 == 0:
        print(i, b)
        
    d = a % 2
    e = [2, -1][d]
    
    a -= d
    while a & 1 == 0:
        a >>= 1
        a *= 3
        b += e
        i += 1
        
    a += d
    j += 1
