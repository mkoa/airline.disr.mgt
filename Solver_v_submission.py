from gurobipy import *
from openpyxl import load_workbook
import pandas as pd
import os
import timeit
import time


def add_t(t1, t2):
    m1, m2 = t1 % 100, t2 % 100
    return 100 * (t1 // 100 + t2 // 100 + (m1 + m2) // 60) + (m1 + m2) % 60

def sub_t(t1, s1):
    m1, m2 = t1 % 100, abs(s1) % 100
    return 100 * (t1 // 100 + s1 // 100 + (m1 - m2) // 60) + (m1 - m2) % 60

def add_d(delay, t0, time_rez):
    delay = (t0 % 100 + delay * time_rez)
    return 100 * (t0 // 100 + delay // 60) + delay % 60


def add_dt(delay, t1, t0, time_rez):
    return add_t(t1, add_d(delay, t0, time_rez))


def wrapper(func, *args, **kwargs): #used for benchmarking
    def wrapped():
        return func(*args, **kwargs)
    return wrapped


def time_function(func, *args, **kwargs): #used for benchmarking
    wrapped = wrapper(func, *args, **kwargs)
    print timeit.timeit(wrapped, number=1)


def solve(scenario_file, t_start, t_end, t_rez):
    # importing initial data
    import_time = time.time() #used for benchmarking
    os.chdir("/Users/moukako/Dropbox/Academic/ISMM/Research Project/Code/Gurobi/Scenarios_Integrated")
    Sch = pd.read_excel(scenario_file, 'Schedule', index_col='f', convert_float=True)
    AC = pd.read_excel(scenario_file, 'Aircrafts', index_col='a', convert_float=True)
    Pax = pd.read_excel(scenario_file, 'Pax', index_col='p', convert_float=True)

    set_start = time.time()  #used for benchmarking

    # creating sets
    F = [f for f in Sch.index if Sch.d_STA[f] > t_start and Sch.d_STD[f] < t_end]
    A = pd.unique(AC.index.ravel())
    P = pd.unique(Pax.index.ravel())
    N = pd.unique(Sch[['Origin', 'Dest']].values.ravel())
    T = [t for t in xrange(t_start, t_end + t_rez, t_rez) if t % 100 < 60]
    D = range(0, len(T))
    M = pd.unique(AC.Type.values.ravel())

    # determining entering, exiting and inner flights
    Fenter = [f for f in F if Sch.d_STD[f] < t_start < Sch.d_STA[f]]
    Fexit = [f for f in F if Sch.STD[f] < t_end < Sch.STA[f]]
    Fins1 = [f for f in F if Sch.STD[f] >= t_start and Sch.d_STA[f] <= t_end]
    Fins2 = [f for f in F if Sch.STD[f] < t_start < Sch.d_STD[f]]

    # determining start & end positions of aircrafts and passenger groups
    Asource = {a: (Sch.Origin[Sch.d_STD[(f for f in (AC.Ro[a].split(',')) if f in F and ~Sch.is_Canx[f])].argmin()])
                for a in A}
    Asink = {a: (Sch.Dest[Sch.d_STD[(f for f in (AC.Ro[a].split(',')) if f in F and ~Sch.is_Canx[f])].argmax()])
              for a in A}
    Psource = {p: (Sch.Origin[Sch.d_STD[(f for f in (Pax.It[p].split(',')) if f in F and ~Sch.is_Canx[f])].argmin()])
                for p in P}
    Psink = {p: (Sch.Dest[Sch.d_STD[(f for f in (Pax.It[p].split(',')) if f in F and ~Sch.is_Canx[f])].argmax()])
              for p in P}
    P_STA = {p: Sch.STA[(f for f in (Pax.It[p].split(',')) if f in F)].max() for p in P}
    P_firstF = {p: Sch.STD[(f for f in (Pax.It[p].split(',')) if f in F)].argmin() for p in P }
    Mend_n = {(n, m): sum(1 if Asink[a] == n and AC.Type[a] == m else 0 for a in A) for n in N for m in M}

    # determining recovery options for last flights
    R = {p: [f for f in F if Sch.Dest[f] == Psink[p] and
             Sch.STD[f] >= Sch.STD[P_firstF[p]]] for p in P}

    # creating all possible arcs
    ArcAC = [(a, f, d, Sch.Origin[f], add_d(d, Sch.STD[f], t_rez), Sch.Dest[f], d_STA)
             for a in A for f, d, d_STA in
             ((f, d, add_dt(d, AC.TAT[a], Sch.STA[f], t_rez)) for f in Fins1 for d in D)
             if d_STA <= t_end]
    ArcAC += [(a, f, 0, Sch.Origin[f], t_start, Sch.Dest[f], add_t(AC.TAT[a], Sch.STA[f])) for a in A for f in Fenter]
    ArcAC += [(a, f, 0, Sch.Origin[f], Sch.STD[f], Sch.Dest[f], t_end)for a in A for f in Fexit]
    ArcAC += [(a, f, d, Sch.Origin[f], d_STD, Sch.Dest[f], d_STA)
              for a in A for f, d, d_STD, d_STA in
              ((f, d, add_d(d, Sch.STD[f], t_rez), add_dt(d, AC.TAT[a], Sch.STA[f], t_rez)) for f in Fins2 for d in D)
              if d_STD >= t_start and d_STA <= t_end]
    ArcPax = [(p, f, d, Sch.Origin[f], add_d(d, Sch.STD[f], t_rez), Sch.Dest[f], d_STA)
              for p in P for f, d, d_STA in
              ((f, d, add_dt(d, Pax.CXT[p], Sch.STA[f], t_rez)) for f in Fins1 for d in D
               if Sch.STD[f] >= Sch.STD[P_firstF[p]])
              if d_STA <= t_end]
    ArcPax += [(p, f, 0, Sch.Origin[f], t_start, Sch.Dest[f], add_t(Pax.CXT[p], Sch.STA[f])) for p in P for f in Fenter
               if Sch.STD[f] >= Sch.STD[P_firstF[p]]]
    ArcPax += [(p, f, 0, Sch.Origin[f], Sch.STD[f], Sch.Dest[f], t_end)for p in P for f in Fexit
               if Sch.STD[f] >= Sch.STD[P_firstF[p]]]
    ArcPax += [(p, f, d, Sch.Origin[f], d_STD, Sch.Dest[f], d_STA)
               for p in P for f, d, d_STD, d_STA in
               ((f, d, add_d(d, Sch.STD[f], t_rez), add_dt(d, Pax.CXT[p], Sch.STA[f], t_rez)) for f in Fins2 for d in D
                if Sch.STD[f] >= Sch.STD[P_firstF[p]])
               if d_STD >= t_start and d_STA <= t_end]
    ArcAC, ArcPax = tuplelist(ArcAC), tuplelist(ArcPax)

    # creating variables
    m = Model('AC+Pax')
    dFCanx = {f: m.addVar(vtype=GRB.BINARY, name="dFCanx[{}]".format(f)) for f in F}
    dAF = {(a, f, d): m.addVar(vtype=GRB.BINARY, name="dAF[{},{},{}]".format(a, f, d)) for a, f, d, no, to, nd, td in ArcAC}
    PCanx = {p: m.addVar(vtype=GRB.INTEGER, name="PCanx[{}]".format(p)) for p in P}
    DTot = {p: m.addVar(vtype=GRB.INTEGER, name="DTot[{}]".format(p)) for p in P}
    PF = {(p, f, d): m.addVar(vtype=GRB.INTEGER, name="PF[{},{},{}]".format(p, f, d)) for p, f, d, no, to, nd, td in ArcPax}
    dAG = {(a, n, t): m.addVar(vtype=GRB.BINARY, name="dAG[{},{},{}]".format(a, n, t)) for a in A for n in N for t in T}
    PG = {(p, n, t): m.addVar(vtype=GRB.INTEGER, name="PG[{},{},{}]".format(p, n, t)) for p in P for n in N for t in T}
    m.update()

    # importing cost parameters
    COP = pd.read_excel(scenario_file, 'COP', index_col='a', convert_float=True)
    CDp = pd.read_excel(scenario_file, 'CDp', index_col='p', convert_float=True)
    CDf = pd.read_excel(scenario_file, 'CDf', index_col='d', convert_float=True)
    CAG = pd.read_excel(scenario_file, 'CAG', index_col='a', convert_float=True)
    CPG = pd.read_excel(scenario_file, 'CPG', index_col='p', convert_float=True)
    CC = pd.read_excel(scenario_file, 'CC', index_col='p', convert_float=True)

    # setting objective function
    m.setObjective(quicksum(dAF[a, f, d] * (COP[f][a] + CDf['CDf'][d]) for a, f, d in dAF.keys()) +
                   quicksum(dAG[a, n, t] * CAG[n][a] for a, n, t in dAG.keys()) +
                   quicksum(DTot[p] * CDp['CDp'][p] for p in DTot.keys()) +
                   quicksum(PCanx[p] * CC['CCp'][p] for p in PCanx.keys()) +
                   quicksum(PG[p, n, t] * CPG[n][p] for p, n, t in PG.keys()))

    # constraint 1: flight coverage
    for f0 in F:
        m.addConstr(quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td in ArcAC.select('*', f0, '*', '*', '*', '*', '*')) +
                    dFCanx[f0] == 1, name="F_coverage_{}".format(f0))

    # constraint 2: passenger coverage
    for p0 in P:
        m.addConstr(quicksum(PF[p, f, d] for p, f, d, no, to, nd, td in ArcPax.select(p0, '*', '*', '*', '*', '*', '*')) +
                    PCanx[p0] >= Pax.Npax[p0], name="Pax_coverage_{}".format(p0))

    # constraint 3: PF & AF binding
    for f0 in F:
        for d0 in D:
            try:
                m.addConstr(quicksum(PF[p, f0, d0] for p, f, d, no, to, nd, td in ArcPax.select('*', f0, d0)) <=
                            quicksum(AC.Cap[a] * dAF[a, f0, d0] for a in A), name="AF&PF_binding_{}_{}".format(f0, d0))
            except:
                continue

    # constraint 4: DTot & PF binding
    for p0 in P:
        m.addConstr(quicksum(PF[p0, r, d] * (add_d(d, Sch.STA[r], t_rez) - P_STA[p0])
                             for r in R[p0] for p, f, d, no, to, nd, td in ArcPax.select(p0, r, '*', '*', '*', '*', '*')) <=
                    DTot[p0], name="DTot&PF_binding_{}".format(p0))

    # constraint 5: continuity at source and sink nodes
    for a0 in A:
        m.addConstr(dAG[a0, Asource[a0], t_start] +
                    quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td
                             in ArcAC.select(a0, '*', '*', Asource[a0], t_start, '*', '*')) == 1,
                    name="SourceAC_{}".format(a0))
    for n0 in N:
        for m0 in M:
            m.addConstr(quicksum(dAG[a0, n0, add_d(-1, t_end, t_rez)] +
                                 quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td
                                          in ArcAC.select(a0, '*', '*', '*', '*', n0, t_end))
                                 for a0 in A if AC.Type[a0] == m0) >= Mend_n[n0, m0], name="SinkAC_{}_{}".format(n0, m0))
    for p0 in P:
        m.addConstr(PG[p0, Psource[p0], t_start] +
                    quicksum(PF[p, f, d] for p, f, d, no, to, nd, td
                             in ArcPax.select(p0, '*', '*', Psource[p0], t_start, '*', '*')) +
                    PCanx[p0] == Pax.Npax[p0], name="SourcePax_{}".format(p0))
        m.addConstr(PG[p0, Psink[p0], add_d(-1, t_end, t_rez)] +
                    quicksum(PF[p, f, d] for p, f, d, no, to, nd, td
                             in ArcPax.select(p0, '*', '*', '*', '*', Psink[p0], t_end)) +
                    PCanx[p0] == Pax.Npax[p0], name="SinkPax_{}".format(p0))

    # constraint 6: continuity inside the graph
    for a0 in A:
        for (n0, t0) in ((n, t) for n in N for t in T if (n, t) not in ((Asource[a0], t_start), (Asink[a0], t_end))):
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
                m.addConstr(dAG[a0, n0, add_d(-1, t0, t_rez)] +
                            quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td
                                     in ArcAC.select(a0, '*', '*', '*', '*', n0, t0)) -
                            dAG[a0, n0, t0] -
                            quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td
                                     in ArcAC.select(a0, '*', '*', n0, t0, '*', '*')) == 0,
                            name="A_continuity_{}_{}_{}".format(a0, n0, t0))
    for p0 in P:
        for (n0, t0) in ((n, t) for n in N for t in T if (n, t) not in ((Psource[p0], t_start), (Psink[p0], t_end))):
            if t0 == t_start:
                m.addConstr(quicksum(PF[p, f, d] for p, f, d, no, to, nd, td
                                     in ArcPax.select(p0, '*', '*', '*', '*', n0, t0)) -
                            PG[p0, n0, t0] -
                            quicksum(PF[p, f, d] for p, f, d, no, to, nd, td
                                     in ArcPax.select(p0, '*', '*', n0, t0, '*', '*')) == 0,
                            name="Pax_continuity_{}_{}_{}".format(p0, n0, t0))
            elif t0 == t_end:
                m.addConstr(PG[p0, n0, add_d(-1, t0, t_rez)] +
                            quicksum(PF[p, f, d] for p, f, d, no, to, nd, td
                                     in ArcPax.select(p0, '*', '*', '*', '*', n0, t0)) -
                            quicksum(PF[p, f, d] for p, f, d, no, to, nd, td
                                     in ArcPax.select(p0, '*', '*', n0, t0, '*', '*')) == 0,
                            name="Pax_continuity_{}_{}_{}".format(p0, n0, t0))
            else:
                m.addConstr(PG[p0, n0, add_d(-1, t0, t_rez)] +
                            quicksum(PF[p, f, d] for p, f, d, no, to, nd, td
                                     in ArcPax.select(p0, '*', '*', '*', '*', n0, t0)) -
                            PG[p0, n0, t0] -
                            quicksum(PF[p, f, d] for p, f, d, no, to, nd, td
                                     in ArcPax.select(p0, '*', '*', n0, t0, '*', '*')) == 0,
                            name="Pax_continuity_{}_{}_{}".format(p0, n0, t0))

    # constraint 7: flight delay, canx or forced time
    Fdisr = pd.read_excel(scenario_file, 'Flight disruptions', index_col='f', convert_float=True)
    for f0 in F:
        if pd.notnull(Fdisr.Delay[f0]):
            m.addConstr(quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td in ArcAC.select('*', f0)
                                 if to < add_t(Fdisr.Delay[f0], Sch.STD[f0])) == 0, name="F_delay_{}".format(f0))
        if pd.notnull(Fdisr.Forced_t[f0]):
            m.addConstr(quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td in ArcAC.select('*', f0)
                                 if to == Fdisr.Forced_t[f0]) == 1, name="F_forced_t_{}".format(f0))
        if pd.notnull(Fdisr.Forced_canx[f0]):
            m.addConstr(dFCanx[f0] == 1, name="F_forced_canx_{}".format(f0))

    # Aircraft unavailable
    Aunav = pd.read_excel(scenario_file, 'Aircrafts unavailable', index_col='a', convert_float=True)
    for a0 in A:
        if pd.notnull(Aunav.From[a0]):
            m.addConstr(quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td in ArcAC.select(a0)
                                 if Aunav.From[a0] <= to < Aunav.Until[a0] or
                                 Aunav.From[a0] <= sub_t(td, AC.TAT[a0]) < Aunav.Until[a0]) == 0,
                        name="P_unavailable_{}".format(a0))

    # Airport unavailable
    Nunav = pd.read_excel(scenario_file, 'Airports unavailable', index_col='n', convert_float=True)
    for n0 in N:
        if pd.notnull(Nunav.From[n0]):
            m.addConstr(quicksum(dAF[a, f, d] for a, f, d, no, to, nd, td in ArcAC.select('*', '*', '*', n0)
                                 if Aunav.From[n0] <= to < Aunav.Until[n0] or
                                 Aunav.From[n0] <= td < Aunav.Until[n0]) == 0, name="P_unavailable_{}".format(a0))

    set_end = time.time()  #used for benchmarking

    # updating the model and solving using Gurobi
    m.setParam('OutputFlag', False)
    m.update()
    m.optimize()
    m.write('AC+Pax.lp')

    solve_end = time.time()  #used for benchmarking

    for f in F:
        if dFCanx[f].x == 1:
            #print dFCanx[f].VarName, '=', dFCanx[f].x
            Sch.loc[f, 'a'] = a
            Sch.loc[f, 'd_STD'] = add_t(Sch['STD'][f], 600) #max delay for canx flights
            Sch.loc[f, 'd_STA'] = add_t(Sch['STA'][f], 600)
            Sch.loc[f, 'is_Canx'] = True
    for a0 in A:
        Ro = ''
        for a, f, d, no, to, nd, td in ArcAC.select(a0, '*', '*', '*', '*', '*', '*'):
            if dAF[a, f, d].x == 1:
                #print dAF[a, f, d].VarName, '=', dAF[a, f, d].x
                Sch.loc[f, 'a'] = a
                Sch.loc[f, 'd_STD'] = to
                Sch.loc[f, 'd_STA'] = sub_t(td, AC.TAT[a])
                Sch.loc[f, 'is_Canx'] = False
                if Ro == '':
                    Ro += f
                else:
                    Ro += ','+ f
                AC.Ro = Ro
    for p0 in P:
        It = ''
        for p, f, d, no, to, nd, td in ArcPax.select(p0, '*', '*', '*', '*', '*', '*'):
            if PF[p, f, d].x > 0:
                # print PF[p, f, d].VarName, '=', PF[p, f, d].x
                Pax.loc[p, 'Npax'] = PF[p, f, d].x
                if It == '':
                    It += f
                else:
                    It += ',' + f

    writer = pd.ExcelWriter('Results.xlsx', engine='openpyxl')
    book = load_workbook('Results.xlsx')
    writer.book = book
    writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
    Sch.to_excel(writer, sheet_name='Schedule')
    AC.to_excel(writer, sheet_name='AC')
    Pax.to_excel(writer, sheet_name='Pax')
    writer.save()

    write_end = time.time()  #used for benchmarking
    for v in m._Model__vars: # prints results if need be
        if v.x > 0 and "dAG" not in v.VarName and "PG" not in v.VarName:
            print v.VarName, "=", v.x
    print set_start-import_time, set_end - set_start, solve_end - set_end, write_end - solve_end