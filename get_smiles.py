import pandas
import pubchempy as pcp
ivo = pandas.read_csv("~/yaupon/metabolomics/metabolite_names.csv")


results = pcp.get_compounds('Saikosaponin BK1', 'name')
print(results)
for compound in results:
    print(compound.canonical_smiles)

def func(x):
    #print(x)
    results = pcp.get_compounds(x,'name')
    for compound in results:
        return compound.canonical_smiles


ivo['Compound'] = ivo['Name'].apply(func)
ivo.to_csv("smiles_raw.csv")
