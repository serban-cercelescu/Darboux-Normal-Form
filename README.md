# Darboux-Normal-Form
The input format is:
- the first line consists of two numbers: `mod` and `n`
- the next `n` lines consist of a list of `n` numbers between `0` and `mod-1` each, describing an `n x n` skew-symmetric matrix

Example interaction:

stdin:
```
6 4
0 1 2 3
5 0 1 2
4 5 0 1
3 4 5 0
```

stdout:
``` 
Darboux form: 
0 5 0 0 
1 0 0 0 
0 0 0 0 
0 0 0 0 
Change of basis matrix: 
1 0 5 0 
0 0 0 1 
0 1 3 4 
0 5 4 1
```