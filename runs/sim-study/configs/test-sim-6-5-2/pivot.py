import pandas as pd

metrics = pd.read_csv('metrics.csv')
p = metrics.pivot_table(index=['dataseed', 'mcmcseed', 'scale', 'DIC', 'LPML'],
                        columns='i')
p.reset_index().round(2).to_csv('metrics_pivot_table.csv', index=False)
