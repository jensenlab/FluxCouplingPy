#!/bin/env python
import random
import numpy as np
import pandas as pd
from gurobipy import *

def optimize(model, maximize, obj):
    v = model.getVars()[obj]
    if maximize:
        model.setObjective(v, GRB.MAXIMIZE)
    if not maximize:
        model.setObjective(v, GRB.MINIMIZE)

    model.update()
    model.optimize()
    return(np.array(model.X))

def pmax(a, b):
    replace = np.where(b > a)
    a[replace] = b[replace]
    return(a)

def pmin(a, b):
    replace = np.where(b < a)
    a[replace] = b[replace]
    return(a)

def cachedFCF(model, min_fva_cor=0.9, fix_frac=0.1, fix_tol_frac=0.01,
      bnd_tol = 0.1, stored_obs = 4000, cor_iter = 5, cor_check = True, idxs = []):
      model.setParam('OutputFlag', False)

      # previous values
      prev_obj = model.Obj

      vars = model.getVars()
      n = len(vars)

      # reset model
      #model.Obj = np.zeros(n)

      n_idxs = len(idxs)

      if (n_idxs < 1):
          idxs = range(0,n)
          n_idxs = n

      coupled = np.zeros((n, n), dtype=bool)

      global_min = np.zeros(n)
      global_max = np.zeros(n)

      blocked = np.zeros(n, dtype=bool)
      active = np.ones(n, dtype=bool)

      flux = np.zeros([stored_obs, n])
      lp_calls = 0

      def not_fixed(x,y):
          if ((not np.isinf(x)) and (not np.isinf(y)) and (abs(x-y) > fix_tol_frac*max(x,y))):
              return(True)
          return(False)

      def fix_flux(max_, min_):
          avg = min_ + fix_frac*(max_ - min_)
          ct = 0
          while (np.isclose(avg, 0) & (ct < 5)):
              avg = avg  + fix_frac*(max_ - min_)
              ct = ct + 1
              print(ct)
          return(avg)

      def correlation_check(flux, i, j):
          def not_near(arr, val):
              close = np.isclose(arr, val)
              not_close = np.where(close == False)
              return(not_close)

          entries_i = not_near(flux[:,i], 0)
          entries_j = not_near(flux[:,j], 0)
          entries_ij = np.concatenate((entries_i, entries_j), axis=1)
          n_entries = len(np.unique(entries_ij))
          if (n_entries > cor_iter):
              C = np.corrcoef(flux[:,i], flux[:,j])[1,0]
              if (abs(C) < min_fva_cor):
                  return(False)
          return(True)

      def update_flux(flux_, idx, sol):
          if (stored_obs > 0):
              flux_[idx,] = sol
          return(flux_)

      for i in range(0,n_idxs):
          idx_i = idxs[i]

          sub_max = np.full(n, np.inf)
          sub_min = np.full(n, -1*np.inf)
          prev_lb = model.LB[idx_i]
          prev_ub = model.UB[idx_i]

          if ((not np.isclose(global_max[i], 0)) or (not np.isclose(global_min[i], 0))):
              fixed_val = fix_flux(global_max[i], global_min[i])
          else:
              if (not np.isclose(model.UB[i], 0)):
                  sol = optimize(model, True, idx_i)
                  lp_calls = lp_calls + 1
                  global_max = pmax(global_max, sol)
                  global_min = pmin(global_min, sol)
                  flux[lp_calls%stored_obs,] = sol
                  if ((not np.isclose(global_max[i], 0)) or (not np.isclose(global_min[i], 0))):
                      fixed_val = fix_flux(global_max[i], global_min[i])
              if (not np.isclose(model.LB[i], 0)):
                  sol = optimize(model, False, idx_i)
                  lp_calls = lp_calls + 1
                  global_max = pmax(global_max, sol)
                  global_min = pmin(global_min, sol)
                  flux[lp_calls%stored_obs,] = sol
                  if ((not np.isclose(global_max[i], 0)) or (not np.isclose(global_min[i], 0))):
                      fixed_val = fix_flux(global_max[i], global_min[i])
              if (np.isclose(global_max[i], 0) or np.isclose(global_min[i], 0)):
                  blocked[idx_i] = True
                  active[idx_i] = False
                  continue

          for j in range(i,n_idxs):
              idx_j = idxs[j]
              # fill this in

              if (blocked[idx_j] | (not active[idx_j])):
                  continue

              if (not_fixed(sub_max[idx_j], sub_min[idx_j])):
                  continue
              if (cor_check):
                  if (not correlation_check(flux, idx_i, idx_j)):
                      continue
              skip = False
              max = 0
              min = 0

              if (idx_i == idx_j):
                  #print('default coupled')
                  coupled[idx_i, idx_j] = True
                  continue

              if (not skip):
                  sol = optimize(model, True, idx_j)
                  lp_calls = lp_calls + 1
                  global_max = pmax(global_max, sol)
                  global_min = pmin(global_min, sol)
                  sub_max = pmax(sub_max, sol)
                  sub_min = pmin(sub_min, sol)
                  flux[lp_calls%stored_obs,] = sol
                  max = sol[idx_j]
                  skip = not_fixed(sub_max[idx_j], sub_min[idx_j])
              if (not skip):
                  sol = optimize(model, False, idx_j)
                  lp_calls = lp_calls + 1
                  global_max = pmax(global_max, sol)
                  global_min = pmin(global_min, sol)
                  sub_max = pmax(sub_max, sol)
                  sub_min = pmin(sub_min, sol)
                  flux[lp_calls%stored_obs,] = sol
                  min = sol[idx_j]
                  skip = not_fixed(sub_max[idx_j], sub_min[idx_j])
              if (np.isclose(max, 0) and np.isclose(min, 0)):
                  skip = True
              if (not skip):
                  coupled[idx_i, idx_i] = True
                  coupled[idx_j, idx_j] = True
                  active[idx_j] = False
                  coupled[idx_i, idx_j] = True
      #np.savetxt("foo.csv", coupled, delimiter=",")
      df = pd.DataFrame(data=coupled[0:,0:], columns=vars, index=vars)
      #df.to_csv('python_output.csv', sep='\t')
      #print(coupled)
      return(df)

def sets_from_coupling_mtx(mtx):
    elems_list = np.array(mtx.columns.values)
    coupled = mtx.values
    n = len(elems_list)
    active = np.ones(n, dtype=bool)
    sets = []

    for i in range(0, n):
        if (not active[i]):
            continue
        if (not coupled[i,i]):
            active[i] = False
            continue
        coupled_idxs = np.where(coupled[i,])
        sets.append(list(elems_list[coupled_idxs]))
        active[coupled_idxs] = False
    return(sets)

fname = "/home/dikshant/GitHub/PathwayMining/ecoli.lp"
model = read(fname)
mtx = cachedFCF(model)
sets = sets_from_coupling_mtx(mtx)
print(len(sets))
print(sets)
