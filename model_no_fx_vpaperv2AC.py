from gurobipy import *
from openpyxl import load_workbook
import pandas as pd
import os
import timeit
import time

pd.set_option('display.max_rows', 15)
pd.set_option('display.width', 500)

scenario_file = 'SWISS.xlsx'
t_start = 3600
t_end = 85200
t_rez = 600


# importing initial data
os.chdir("/Users/moukako/Dropbox/Academic/Research Paper/Code/Scenarios")
Sch = pd.read_excel(scenario_file, 'Schedule', convert_float=True)
Sch.index.name = 'f'
AC = pd.read_excel(scenario_file, 'Aircrafts', index_col='a', convert_float=True)

set_start = time.time()  # used for benchmarking

# creating sets
Sch_win = Sch[(Sch['d_STA'] > t_start) & (Sch['d_STD'] < t_end)]
F = Sch_win.index
A = Sch_win['a'].unique()
AC_win = AC[AC.index.isin(A)]
N = Sch_win[['Origin', 'Dest']].unstack().unique()
T = [t for t in xrange(t_start, t_end + t_rez, t_rez)]
D = range(0, len(T))
M = AC_win.loc[AC_win.index.isin(A), 'Type'].unique()

# determining entering, exiting and inner flights
Fenter = Sch_win[(Sch_win['d_STD'] < t_start) & (Sch_win['d_STA'] > t_start)].index
Fexit = Sch_win[(Sch_win['d_STD'] < t_end) & (Sch_win['d_STA'] > t_end)].index
Fins1 = Sch_win[(Sch_win['d_STD'] >= t_start) & (Sch_win['d_STA'] <= t_end)].index
Fins2 = Sch_win[(Sch_win['d_STD'] < t_start) & (Sch_win['d_STD'] > t_start)].index

# determining start & end positions of aircrafts and passenger groups
AC_win.sort_index(inplace=True)
AC_win['Source'] = Sch_win.loc[Sch_win.groupby('a')['d_STD'].idxmin(), 'Origin'].values
AC_win['Sink'] = Sch_win.loc[Sch_win.groupby('a')['d_STA'].idxmax(), 'Dest'].values
Mend = AC_win.groupby(['Sink', 'Type']).size()

# creating all possible arcs
ArcAC = [(a, f, d, Sch_win.loc[f, 'Origin'], Sch_win.loc[f, 'd_STD'] + d*t_rez, Sch_win.loc[f, 'Dest'], d_STA)
         for a in A for f, d, d_STA in ((f, d, AC_win.loc[a, 'TAT'] + Sch_win.loc[f, 'STA'] + d*t_rez)
         for f in Fins1 for d in D) if d_STA <= t_end]
ArcAC += [(a, f, 0, Sch_win.loc[f, 'Origin'], t_start, Sch_win.loc[f, 'Dest'], AC.loc[a, 'TAT'] + Sch_win.loc[f, 'STA'])
          for a in A for f in Fenter]
ArcAC += [(a, f, 0, Sch_win.loc[f, 'Origin'], Sch_win.loc[f, 'STD'], Sch_win.loc[f, 'Dest'], t_end)
          for a in A for f in Fexit]
ArcAC += [(a, f, d, Sch_win.loc[f, 'Origin'], d_STD, Sch_win.loc[f, 'Dest'], d_STA)
          for a in A for f, d, d_STD, d_STA in
          ((f, d, Sch_win.loc[f, 'STD'] + d*t_rez, AC_win.loc[a, 'TAT'] + Sch_win.loc[f, 'STA'] + d*t_rez)
          for f in Fins2 for d in D) if d_STD >= t_start and d_STA <= t_end]

ArcAC = tuplelist(ArcAC)

set_inter = time.time()  # used for benchmarking

# creating variables

m = Model('AC+Pax')
dFCanx = {f: m.addVar(vtype=GRB.BINARY, name="dFCanx[{}]".format(f)) for f in F}
dAF = {(a, f, d): m.addVar(vtype=GRB.BINARY, name="dAF[{},{},{}]".format(a, f, d)) for a, f, d, no, to, nd, td in
       ArcAC}
dAG = {(a, n, t): m.addVar(vtype=GRB.BINARY, name="dAG[{},{},{}]".format(a, n, t)) for a in A for n in N for t in T}
m.update()

# importing cost parameters
COP = pd.read_excel(scenario_file, 'COPaf', index_col='a', convert_float=True)
CD = pd.read_excel(scenario_file, 'CDd', index_col='d', convert_float=True)
CAG = pd.read_excel(scenario_file, 'CAGan', index_col='a', convert_float=True)

# setting objective function
m.setObjective(quicksum(dAF[a, f, d] * (COP.loc[a, f] + CD.loc[d, 'CDd']) for a, f, d in dAF.keys()) +
               quicksum(dAG[a, n, t] * CAG.loc[a, n] for a, n, t in dAG.keys()))

# constraint 1: flight coverage
for f0 in F:
    m.addConstr(
        quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td in ArcAC.select('*', f0, '*', '*', '*', '*', '*')) +
        dFCanx[f0] == 1, name="F_coverage_{}".format(f0))

# constraint 5: continuity at source and sink nodes
for a0 in A:
    m.addConstr(dAG[a0, AC_win.loc[a0, 'Source'], t_start] +
                quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td
                         in ArcAC.select(a0, '*', '*', AC_win.loc[a0, 'Source'], t_start, '*', '*')) == 1,
                name="SourceAC_{}".format(a0))
for n0 in Mend.index.get_level_values(0).unique():
    for m0 in M:
        m.addConstr(quicksum(dAG[a0, n0, t_end - t_rez] +
                             quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td
                                      in ArcAC.select(a0, '*', '*', '*', '*', n0, t_end))
                             for a0 in A if AC_win.loc[a0, 'Type'] == m0) >= Mend.loc[n0, m0],
                    name="SinkAC_{}_{}".format(n0, m0))

# constraint 6: continuity inside the graph
for a0 in A:
    for (n0, t0) in ((n, t) for n in N for t in T if (n, t) not in
            ((AC_win.loc[a0, 'Source'], t_start), (AC_win.loc[a0, 'Sink'], t_end))):
        if t0 == t_start:
            m.addConstr(quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td
                                 in ArcAC.select(a0, '*', '*', '*', '*', n0, t0)) -
                        dAG[a0, n0, t0] -
                        quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td
                                 in ArcAC.select(a0, '*', '*', n0, t0, '*', '*')) == 0,
                        name="A_continuity_{}_{}_{}".format(a0, n0, t0))
        elif t0 == t_end:
            continue
        else:
            m.addConstr(dAG[a0, n0, t0 - t_rez] +
                        quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td
                                 in ArcAC.select(a0, '*', '*', '*', '*', n0, t0)) -
                        dAG[a0, n0, t0] -
                        quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td
                                 in ArcAC.select(a0, '*', '*', n0, t0, '*', '*')) == 0,
                        name="A_continuity_{}_{}_{}".format(a0, n0, t0))

# constraint 7: flight delay, canx or forced time
Fdisr = pd.read_excel(scenario_file, 'Flight disruptions', index_col='f', convert_float=True)
for f0 in F:
    if Fdisr.loc[f0, 'Delay'] > 0:
        m.addConstr(quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td in ArcAC.select('*', f0)
                             if to < Fdisr.loc[f0, 'Delay'] + Sch_win.loc[f0, 'STD']) == 0
                    , name="F_delay_{}".format(f0))
    if pd.notnull(Fdisr.loc[f0, 'Forced_t']):
        m.addConstr(quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td in ArcAC.select('*', f0)
                             if to == Fdisr.loc[f0, 'Forced_t']) == 1
                    , name="F_forced_t_{}".format(f0))
    if Fdisr.loc[f0, 'Forced_canx'] == True:
        m.addConstr(dFCanx[f0] == 1, name="F_forced_canx_{}".format(f0))

# Aircraft unavailable
Aunav = pd.read_excel(scenario_file, 'Aircrafts unavailable', index_col='a', convert_float=True)
for a0 in A:
    if pd.notnull(Aunav.loc[a0, 'From']):
        m.addConstr(quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td in ArcAC.select(a0)
                             if Aunav.loc[a0, 'From'] <= to < Aunav.loc[a0, 'Until'] or
                             Aunav.loc[a0, 'From'] <= td - AC_win.loc[a0, 'TAT'] < Aunav.loc[a0, 'Until']) == 0,
                    name="P_unavailable_{}".format(a0))

# Airport unavailable
Nunav = pd.read_excel(scenario_file, 'Airports unavailable', index_col='n', convert_float=True)
for n0 in N:
    if pd.notnull(Nunav.loc[n0, 'From']):
        m.addConstr(quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td in ArcAC.select('*', '*', '*', n0)
                             if Aunav.loc[n0, 'From'] <= to < Aunav.loc[n0, 'Until'] or
                             Aunav.loc[n0, 'From'] <= td < Aunav.loc[n0, 'Until']) == 0, name="P_unavailable_{}".format(a0))

set_end = time.time()  # used for benchmarking

# updating the model and solving using Gurobi
# m.setParam('OutputFlag', False)
m.update()
m.optimize()
m.write('AC+Pax.lp')

solve_end = time.time()  # used for benchmarking

for v in m._Model__vars:  # prints results if need be
    if v.x > 0 and "dAG" not in v.VarName:
        print v.VarName, "=", v.x
print set_start - import_time, set_inter - set_start, set_end - set_start, solve_end - set_end, write_end - solve_end