var("c")

P.<z> = QQ[]

def A(a, c, e, i):
    if a is None or c is None:
        return (None, None, e)
    
    a, c, f = [
        (a - 2 * c - 3, 3 * c / 2 + 2, 2 * c + 1 <= a - 2),
        (a - 2 * c - 5, 3 * c / 2 + 7 / 2, 2 * c + 1 <= a - 4),
        (None, None, 2 * c + 1 == a - 3),
        (6 * c + 11, 0, 2 * c + 1 == a - 2),
        (None, None, 2 * c + 1 == a - 1),
        (1, a / 4 + c + 1, 2 * c + 1 >= a + 1),
        (3, a / 4 + c + 5 / 2, 2 * c + 1 >= a + 1),
        (3 * a + 2, c - a / 2 + 1 / 2, c >= a)
    ][i]
    
    return a, c, [*e, f]
    
def apply(y, ns):
    x = 1
    for n in ns:
        x, y, _ = A(x, y, [], int(n))
        print(x, y)
        
    return (x, y)


def sol(ex):
    d, e = list(ex.polynomial(P))
    dd, ee = QQ(d).denom(), QQ(e).denom()
    return ee > dd or ee == dd and QQ(d).numer() % dd == QQ(e).numer() % dd

def pos(e):
    sol = solve(e, [c])
    if not sol:
        return
    
    sol = sol[0][0]
    return ">" in str(sol) or sol.rhs() > 0
    
apples = [("", 1, c, [])]
for j in range(5):
    apples = [(s + str(i), *A(x, y, e, i)) for s, x, y, e in apples for i in range(8)]
    apples = [(s, x, y, e) for s, x, y, e in apples if y and sol(y) and pos(e)]
    
[(s, y, solve(e, [c])) for s, x, y, e in apples]
