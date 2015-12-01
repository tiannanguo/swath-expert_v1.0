__author__ = 'Tiannan Guo, ETH Zurich 2015'


from lifelines.datasets import load_waltons
df = load_waltons()
print df.head()

T = df['T']
E = df['E']


from lifelines import KaplanMeierFitter
kmf = KaplanMeierFitter()
kmf.fit(T, event_observed=E)

kmf.survival_function_
kmf.median_
kmf.plot()