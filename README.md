# FluxCouplingPy
Python Implementation of Cached Flux Coupling

```python
fname = "ecoli.lp"
model = read(fname)
mtx = cachedFCF(model)
sets = sets_from_coupling_mtx(mtx)
print(sets)
```
