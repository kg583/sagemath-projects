def S(a, b, X=set()):
    seq = [a, b]
    while True:
        try:
            seq.append(next(d for d in divisors(seq[-1] + seq[-2]) if d not in seq + [*X]))
            
        except StopIteration:
            return seq
		
def is_replete(X):
	checked = set()
	for x in sorted(X, reverse=True):
		if x not in checked:
			for d in divisors(x):
				if d in X:
					checked.add(d)
				else:
					return False
	
	return True
