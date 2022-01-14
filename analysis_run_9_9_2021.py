import pandas as pd
import numpy as np
import scipy.optimize
from matplotlib import pyplot as plt

data_dir = 'E:\\slamseq_run_9_9_2021\\count\\'
results_dir = 'E:\\slamseq_run_9_9_2021\\jons_analysis\\'
data_dir2 = 'E:\\slamseq_run_9_9_2021\\test2\\'

################################################
# original data cleaning for slamdunk output ###
################################################

A3_1_0 = pd.read_csv(data_dir + 'Set1-A3-0_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
A3_1_1 = pd.read_csv(data_dir + 'Set1-A3-1_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
A3_1_2 = pd.read_csv(data_dir + 'Set1-A3-2_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
A3_1_4 = pd.read_csv(data_dir + 'Set1-A3-4_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)

A3_2_0 = pd.read_csv(data_dir + 'Set2-A3-0_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
A3_2_1 = pd.read_csv(data_dir + 'Set2-A3-1_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
A3_2_2 = pd.read_csv(data_dir + 'Set2-A3-2_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
A3_2_4 = pd.read_csv(data_dir + 'Set2-A3-4_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)

A3_3_0 = pd.read_csv(data_dir + 'Set3-A3-0_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
A3_3_1 = pd.read_csv(data_dir + 'Set3-A3-1_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
A3_3_2 = pd.read_csv(data_dir + 'Set3-A3-2_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
A3_3_4 = pd.read_csv(data_dir + 'Set3-A3-4_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)

B2_1_0 = pd.read_csv(data_dir + 'Set1-B2-0_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
B2_1_1 = pd.read_csv(data_dir + 'Set1-B2-1_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
B2_1_2 = pd.read_csv(data_dir + 'Set1-B2-2_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
B2_1_4 = pd.read_csv(data_dir + 'Set1-B2-4_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)

B2_2_0 = pd.read_csv(data_dir + 'Set2-B2-0_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
B2_2_1 = pd.read_csv(data_dir + 'Set2-B2-1_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
B2_2_2 = pd.read_csv(data_dir + 'Set2-B2-2_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
B2_2_4 = pd.read_csv(data_dir + 'Set2-B2-4_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)

B2_3_0 = pd.read_csv(data_dir + 'Set3-B2-0_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
B2_3_1 = pd.read_csv(data_dir + 'Set3-B2-1_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
B2_3_2 = pd.read_csv(data_dir + 'Set3-B2-2_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
B2_3_4 = pd.read_csv(data_dir + 'Set3-B2-4_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)

M17_1_0 = pd.read_csv(data_dir + 'Set1-M17-0_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
M17_1_1 = pd.read_csv(data_dir + 'Set1-M17-1_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
M17_1_2 = pd.read_csv(data_dir + 'Set1-M17-2_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
M17_1_4 = pd.read_csv(data_dir + 'Set1-M17-4_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)

M17_2_0 = pd.read_csv(data_dir + 'Set2-M17-0_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
M17_2_1 = pd.read_csv(data_dir + 'Set2-M17-1_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
M17_2_2 = pd.read_csv(data_dir + 'Set2-M17-2_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
M17_2_4 = pd.read_csv(data_dir + 'Set2-M17-4_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)

M17_3_0 = pd.read_csv(data_dir + 'Set3-M17-0_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
M17_3_1 = pd.read_csv(data_dir + 'Set3-M17-1_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
M17_3_2 = pd.read_csv(data_dir + 'Set3-M17-2_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)
M17_3_4 = pd.read_csv(data_dir + 'Set3-M17-4_slamdunk_mapped_filtered_tcount.tsv', sep='\t', header=2).drop(['multimapCount', 'ConversionRateLower', 'ConversionRateUpper'], axis=1)

A3_a = A3_1_0.merge(A3_1_1, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_0', '_1'))
A3_b = A3_1_2.merge(A3_1_4, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_2', '_4'))
A3_1 = A3_a.merge(A3_b, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer')
A3_1.to_csv(results_dir + 'A3_1_raw.csv', index=False)

B2_a = B2_1_0.merge(B2_1_1, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_0', '_1'))
B2_b = B2_1_2.merge(B2_1_4, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_2', '_4'))
B2_1 = B2_a.merge(B2_b, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer')
B2_1.to_csv(results_dir + 'B2_1_raw.csv', index=False)

M17_a = M17_1_0.merge(M17_1_1, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_0', '_1'))
M17_b = M17_1_2.merge(M17_1_4, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_2', '_4'))
M17_1 = M17_a.merge(M17_b, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer')
M17_1.to_csv(results_dir + 'M17_1_raw.csv', index=False)

A3_a = A3_2_0.merge(A3_2_1, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_0', '_1'))
A3_b = A3_2_2.merge(A3_2_4, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_2', '_4'))
A3_2 = A3_a.merge(A3_b, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer')
A3_2.to_csv(results_dir + 'A3_2_raw.csv', index=False)

B2_a = B2_2_0.merge(B2_2_1, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_0', '_1'))
B2_b = B2_2_2.merge(B2_2_4, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_2', '_4'))
B2_2 = B2_a.merge(B2_b, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer')
B2_2.to_csv(results_dir + 'B2_2_raw.csv', index=False)

M17_a = M17_2_0.merge(M17_2_1, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_0', '_1'))
M17_b = M17_2_2.merge(M17_2_4, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_2', '_4'))
M17_2 = M17_a.merge(M17_b, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer')
M17_2.to_csv(results_dir + 'M17_2_raw.csv', index=False)

A3_a = A3_3_0.merge(A3_3_1, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_0', '_1'))
A3_b = A3_3_2.merge(A3_3_4, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_2', '_4'))
A3_3 = A3_a.merge(A3_b, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer')
A3_3.to_csv(results_dir + 'A3_3_raw.csv', index=False)

B2_a = B2_3_0.merge(B2_3_1, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_0', '_1'))
B2_b = B2_3_2.merge(B2_3_4, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_2', '_4'))
B2_3 = B2_a.merge(B2_b, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer')
B2_3.to_csv(results_dir + 'B2_3_raw.csv', index=False)

M17_a = M17_3_0.merge(M17_3_1, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_0', '_1'))
M17_b = M17_3_2.merge(M17_3_4, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer', suffixes=('_2', '_4'))
M17_3 = M17_a.merge(M17_b, on=['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand'], how='outer')
M17_3.to_csv(results_dir + 'M17_3_raw.csv', index=False)

#############################################
# new data cleaning for grand slam output ###
#############################################

df = pd.read_csv(data_dir2 + '24h.tsv', sep='\t')
df2 = df.filter(regex=r'(Gene|Symbol|Length|Conversion|Coverage)', axis=1)
df2.drop(list(df2.filter(regex = 'Double')), axis = 1, inplace = True)
df2.columns = df2.columns.str.replace('/mnt/e/slamseq_run_9_9_2021/filter/','')

for s in ['Set1','Set2','Set3']:
    for c in ['A3','B2','M17']:
        for t in [0,1,2,4]:
            df2['{0}-{1}-{2}-Conversion-Rate'.format(s,c,t)] = df2['{0}-{1}-{2} Conversions'.format(s,c,t)] / df2['{0}-{1}-{2} Coverage'.format(s,c,t)]

df2.drop(list(df2.filter(regex = 'Conversions')), axis = 1, inplace = True)
df2.drop(list(df2.filter(regex = 'Coverage')), axis = 1, inplace = True)

df2.to_csv(data_dir2 + 'conversion_rates.csv', index=False)

A3_1 = df2.filter(regex=r'(Gene|Symbol|Length|Set1-A3)', axis=1)
A3_2 = df2.filter(regex=r'(Gene|Symbol|Length|Set2-A3)', axis=1)
A3_3 = df2.filter(regex=r'(Gene|Symbol|Length|Set3-A3)', axis=1)

B2_1 = df2.filter(regex=r'(Gene|Symbol|Length|Set1-B2)', axis=1)
B2_2 = df2.filter(regex=r'(Gene|Symbol|Length|Set2-B2)', axis=1)
B2_3 = df2.filter(regex=r'(Gene|Symbol|Length|Set3-B2)', axis=1)

M17_1 = df2.filter(regex=r'(Gene|Symbol|Length|Set1-M17)', axis=1)
M17_2 = df2.filter(regex=r'(Gene|Symbol|Length|Set2-M17)', axis=1)
M17_3 = df2.filter(regex=r'(Gene|Symbol|Length|Set3-M17)', axis=1)

A3_1.columns = ['Gene','Symbol','Length','ConversionRate_0','ConversionRate_1','ConversionRate_2','ConversionRate_4']
A3_2.columns = ['Gene','Symbol','Length','ConversionRate_0','ConversionRate_1','ConversionRate_2','ConversionRate_4']
A3_3.columns = ['Gene','Symbol','Length','ConversionRate_0','ConversionRate_1','ConversionRate_2','ConversionRate_4']

B2_1.columns = ['Gene','Symbol','Length','ConversionRate_0','ConversionRate_1','ConversionRate_2','ConversionRate_4']
B2_2.columns = ['Gene','Symbol','Length','ConversionRate_0','ConversionRate_1','ConversionRate_2','ConversionRate_4']
B2_3.columns = ['Gene','Symbol','Length','ConversionRate_0','ConversionRate_1','ConversionRate_2','ConversionRate_4']

M17_1.columns = ['Gene','Symbol','Length','ConversionRate_0','ConversionRate_1','ConversionRate_2','ConversionRate_4']
M17_2.columns = ['Gene','Symbol','Length','ConversionRate_0','ConversionRate_1','ConversionRate_2','ConversionRate_4']
M17_3.columns = ['Gene','Symbol','Length','ConversionRate_0','ConversionRate_1','ConversionRate_2','ConversionRate_4']

A3_1['ID'] = A3_1['Gene'] +'_'+ A3_1['Symbol'].map(str) +'_'+ A3_1['Length'].map(str)
A3_2['ID'] = A3_2['Gene'] +'_'+ A3_2['Symbol'].map(str) +'_'+ A3_2['Length'].map(str)
A3_3['ID'] = A3_3['Gene'] +'_'+ A3_3['Symbol'].map(str) +'_'+ A3_3['Length'].map(str)

B2_1['ID'] = B2_1['Gene'] +'_'+ B2_1['Symbol'].map(str) +'_'+ B2_1['Length'].map(str)
B2_2['ID'] = B2_2['Gene'] +'_'+ B2_2['Symbol'].map(str) +'_'+ B2_2['Length'].map(str)
B2_3['ID'] = B2_3['Gene'] +'_'+ B2_3['Symbol'].map(str) +'_'+ B2_3['Length'].map(str)

M17_1['ID'] = M17_1['Gene'] +'_'+ M17_1['Symbol'].map(str) +'_'+ M17_1['Length'].map(str)
M17_2['ID'] = M17_2['Gene'] +'_'+ M17_2['Symbol'].map(str) +'_'+ M17_2['Length'].map(str)
M17_3['ID'] = M17_3['Gene'] +'_'+ M17_3['Symbol'].map(str) +'_'+ M17_3['Length'].map(str)

#############################################
# fit then average
#############################################

A3_1 = pd.read_csv(results_dir + 'A3_1_raw.csv')
B2_1 = pd.read_csv(results_dir + 'B2_1_raw.csv')
M17_1 = pd.read_csv(results_dir + 'M17_1_raw.csv')

A3_2 = pd.read_csv(results_dir + 'A3_2_raw.csv')
B2_2 = pd.read_csv(results_dir + 'B2_2_raw.csv')
M17_2 = pd.read_csv(results_dir + 'M17_2_raw.csv')

A3_3 = pd.read_csv(results_dir + 'A3_3_raw.csv')
B2_3 = pd.read_csv(results_dir + 'B2_3_raw.csv')
M17_3 = pd.read_csv(results_dir + 'M17_3_raw.csv')

waterfall = {}

waterfall['A3_1_init'] = A3_1.shape[0]
waterfall['B2_1_init'] = B2_1.shape[0]
waterfall['M17_1_init'] = M17_1.shape[0]

waterfall['A3_2_init'] = A3_2.shape[0]
waterfall['B2_2_init'] = B2_2.shape[0]
waterfall['M17_2_init'] = M17_2.shape[0]

waterfall['A3_3_init'] = A3_3.shape[0]
waterfall['B2_3_init'] = B2_3.shape[0]
waterfall['M17_3_init'] = M17_3.shape[0]

# remove fragments where 4 timepoint is 0

A3_1_f = A3_1.loc[A3_1['ConversionRate_4'] > 0]
B2_1_f = B2_1.loc[B2_1['ConversionRate_4'] > 0]
M17_1_f = M17_1.loc[M17_1['ConversionRate_4'] > 0]

A3_2_f = A3_2.loc[A3_2['ConversionRate_4'] > 0]
B2_2_f = B2_2.loc[B2_2['ConversionRate_4'] > 0]
M17_2_f = M17_2.loc[M17_2['ConversionRate_4'] > 0]

A3_3_f = A3_3.loc[A3_1['ConversionRate_4'] > 0]
B2_3_f = B2_3.loc[B2_1['ConversionRate_4'] > 0]
M17_3_f = M17_3.loc[M17_1['ConversionRate_4'] > 0]

waterfall['A3_1_conv_4>0'] = A3_1_f.shape[0]
waterfall['B2_1_conv_4>0'] = B2_1_f.shape[0]
waterfall['M17_1_conv_4>0'] = M17_1_f.shape[0]

waterfall['A3_2_conv_4>0'] = A3_2_f.shape[0]
waterfall['B2_2_conv_4>0'] = B2_2_f.shape[0]
waterfall['M17_2_conv_4>0'] = M17_2_f.shape[0]

waterfall['A3_3_conv_4>0'] = A3_3_f.shape[0]
waterfall['B2_3_conv_4>0'] = B2_3_f.shape[0]
waterfall['M17_3_conv_4>0'] = M17_3_f.shape[0]

# remove fragments where 4 timepoint is the only positive conversion rate

A3_1_ff = A3_1_f.loc[A3_1_f['ConversionRate_0']+A3_1_f['ConversionRate_1']+A3_1_f['ConversionRate_2'] != 0]
B2_1_ff = B2_1_f.loc[B2_1_f['ConversionRate_0']+B2_1_f['ConversionRate_1']+B2_1_f['ConversionRate_2'] != 0]
M17_1_ff = M17_1_f.loc[M17_1_f['ConversionRate_0']+M17_1_f['ConversionRate_1']+M17_1_f['ConversionRate_2'] != 0]

A3_2_ff = A3_2_f.loc[A3_2_f['ConversionRate_0']+A3_2_f['ConversionRate_1']+A3_2_f['ConversionRate_2'] != 0]
B2_2_ff = B2_2_f.loc[B2_2_f['ConversionRate_0']+B2_2_f['ConversionRate_1']+B2_2_f['ConversionRate_2'] != 0]
M17_2_ff = M17_2_f.loc[M17_2_f['ConversionRate_0']+M17_2_f['ConversionRate_1']+M17_2_f['ConversionRate_2'] != 0]

A3_3_ff = A3_3_f.loc[A3_3_f['ConversionRate_0']+A3_3_f['ConversionRate_1']+A3_3_f['ConversionRate_2'] != 0]
B2_3_ff = B2_3_f.loc[B2_3_f['ConversionRate_0']+B2_3_f['ConversionRate_1']+B2_3_f['ConversionRate_2'] != 0]
M17_3_ff = M17_3_f.loc[M17_3_f['ConversionRate_0']+M17_3_f['ConversionRate_1']+M17_3_f['ConversionRate_2'] != 0]

waterfall['A3_1_conv_x>0'] = A3_1_ff.shape[0]
waterfall['B2_1_conv_x>0'] = B2_1_ff.shape[0]
waterfall['M17_1_conv_x>0'] = M17_1_ff.shape[0]

waterfall['A3_2_conv_x>0'] = A3_2_ff.shape[0]
waterfall['B2_2_conv_x>0'] = B2_2_ff.shape[0]
waterfall['M17_2_conv_x>0'] = M17_2_ff.shape[0]

waterfall['A3_3_conv_x>0'] = A3_3_ff.shape[0]
waterfall['B2_3_conv_x>0'] = B2_3_ff.shape[0]
waterfall['M17_3_conv_x>0'] = M17_3_ff.shape[0]

# normalize to timepoint 4

A3_1_ff['ConversionRate_4_n'] = 1
A3_1_ff['ConversionRate_2_n'] = A3_1_ff['ConversionRate_2'] / A3_1_ff['ConversionRate_4']
A3_1_ff['ConversionRate_1_n'] = A3_1_ff['ConversionRate_1'] / A3_1_ff['ConversionRate_4']
A3_1_ff['ConversionRate_0_n'] = A3_1_ff['ConversionRate_0'] / A3_1_ff['ConversionRate_4']

B2_1_ff['ConversionRate_4_n'] = 1
B2_1_ff['ConversionRate_2_n'] = B2_1_ff['ConversionRate_2'] / B2_1_ff['ConversionRate_4']
B2_1_ff['ConversionRate_1_n'] = B2_1_ff['ConversionRate_1'] / B2_1_ff['ConversionRate_4']
B2_1_ff['ConversionRate_0_n'] = B2_1_ff['ConversionRate_0'] / B2_1_ff['ConversionRate_4']

M17_1_ff['ConversionRate_4_n'] = 1
M17_1_ff['ConversionRate_2_n'] = M17_1_ff['ConversionRate_2'] / M17_1_ff['ConversionRate_4']
M17_1_ff['ConversionRate_1_n'] = M17_1_ff['ConversionRate_1'] / M17_1_ff['ConversionRate_4']
M17_1_ff['ConversionRate_0_n'] = M17_1_ff['ConversionRate_0'] / M17_1_ff['ConversionRate_4']

A3_2_ff['ConversionRate_4_n'] = 1
A3_2_ff['ConversionRate_2_n'] = A3_2_ff['ConversionRate_2'] / A3_2_ff['ConversionRate_4']
A3_2_ff['ConversionRate_1_n'] = A3_2_ff['ConversionRate_1'] / A3_2_ff['ConversionRate_4']
A3_2_ff['ConversionRate_0_n'] = A3_2_ff['ConversionRate_0'] / A3_2_ff['ConversionRate_4']

B2_2_ff['ConversionRate_4_n'] = 1
B2_2_ff['ConversionRate_2_n'] = B2_2_ff['ConversionRate_2'] / B2_2_ff['ConversionRate_4']
B2_2_ff['ConversionRate_1_n'] = B2_2_ff['ConversionRate_1'] / B2_2_ff['ConversionRate_4']
B2_2_ff['ConversionRate_0_n'] = B2_2_ff['ConversionRate_0'] / B2_2_ff['ConversionRate_4']

M17_2_ff['ConversionRate_4_n'] = 1
M17_2_ff['ConversionRate_2_n'] = M17_2_ff['ConversionRate_2'] / M17_2_ff['ConversionRate_4']
M17_2_ff['ConversionRate_1_n'] = M17_2_ff['ConversionRate_1'] / M17_2_ff['ConversionRate_4']
M17_2_ff['ConversionRate_0_n'] = M17_2_ff['ConversionRate_0'] / M17_2_ff['ConversionRate_4']

A3_3_ff['ConversionRate_4_n'] = 1
A3_3_ff['ConversionRate_2_n'] = A3_3_ff['ConversionRate_2'] / A3_3_ff['ConversionRate_4']
A3_3_ff['ConversionRate_1_n'] = A3_3_ff['ConversionRate_1'] / A3_3_ff['ConversionRate_4']
A3_3_ff['ConversionRate_0_n'] = A3_3_ff['ConversionRate_0'] / A3_3_ff['ConversionRate_4']

B2_3_ff['ConversionRate_4_n'] = 1
B2_3_ff['ConversionRate_2_n'] = B2_3_ff['ConversionRate_2'] / B2_3_ff['ConversionRate_4']
B2_3_ff['ConversionRate_1_n'] = B2_3_ff['ConversionRate_1'] / B2_3_ff['ConversionRate_4']
B2_3_ff['ConversionRate_0_n'] = B2_3_ff['ConversionRate_0'] / B2_3_ff['ConversionRate_4']

M17_3_ff['ConversionRate_4_n'] = 1
M17_3_ff['ConversionRate_2_n'] = M17_3_ff['ConversionRate_2'] / M17_3_ff['ConversionRate_4']
M17_3_ff['ConversionRate_1_n'] = M17_3_ff['ConversionRate_1'] / M17_3_ff['ConversionRate_4']
M17_3_ff['ConversionRate_0_n'] = M17_3_ff['ConversionRate_0'] / M17_3_ff['ConversionRate_4']

'''
A3_1_ff['ID'] = A3_1_ff['Chromosome'] +'_'+ A3_1_ff['Start'].map(str) +'_'+ A3_1_ff['End'].map(str) +'_'+ A3_1_ff['Name']
A3_2_ff['ID'] = A3_2_ff['Chromosome'] +'_'+ A3_2_ff['Start'].map(str) +'_'+ A3_2_ff['End'].map(str) +'_'+ A3_2_ff['Name']
A3_3_ff['ID'] = A3_3_ff['Chromosome'] +'_'+ A3_3_ff['Start'].map(str) +'_'+ A3_3_ff['End'].map(str) +'_'+ A3_3_ff['Name']

B2_1_ff['ID'] = B2_1_ff['Chromosome'] +'_'+ B2_1_ff['Start'].map(str) +'_'+ B2_1_ff['End'].map(str) +'_'+ B2_1_ff['Name']
B2_2_ff['ID'] = B2_2_ff['Chromosome'] +'_'+ B2_2_ff['Start'].map(str) +'_'+ B2_2_ff['End'].map(str) +'_'+ B2_2_ff['Name']
B2_3_ff['ID'] = B2_3_ff['Chromosome'] +'_'+ B2_3_ff['Start'].map(str) +'_'+ B2_3_ff['End'].map(str) +'_'+ B2_3_ff['Name']

M17_1_ff['ID'] = M17_1_ff['Chromosome'] +'_'+ M17_1_ff['Start'].map(str) +'_'+ M17_1_ff['End'].map(str) +'_'+ M17_1_ff['Name']
M17_2_ff['ID'] = M17_2_ff['Chromosome'] +'_'+ M17_2_ff['Start'].map(str) +'_'+ M17_2_ff['End'].map(str) +'_'+ M17_2_ff['Name']
M17_3_ff['ID'] = M17_3_ff['Chromosome'] +'_'+ M17_3_ff['Start'].map(str) +'_'+ M17_3_ff['End'].map(str) +'_'+ M17_3_ff['Name']
'''

A3_1['ID'] = A3_1['Chromosome'] +'_'+ A3_1['Start'].map(str) +'_'+ A3_1['End'].map(str) +'_'+ A3_1['Name']
A3_2['ID'] = A3_2['Chromosome'] +'_'+ A3_2['Start'].map(str) +'_'+ A3_2['End'].map(str) +'_'+ A3_2['Name']
A3_3['ID'] = A3_3['Chromosome'] +'_'+ A3_3['Start'].map(str) +'_'+ A3_3['End'].map(str) +'_'+ A3_3['Name']

B2_1['ID'] = B2_1['Chromosome'] +'_'+ B2_1['Start'].map(str) +'_'+ B2_1['End'].map(str) +'_'+ B2_1['Name']
B2_2['ID'] = B2_2['Chromosome'] +'_'+ B2_2['Start'].map(str) +'_'+ B2_2['End'].map(str) +'_'+ B2_2['Name']
B2_3['ID'] = B2_3['Chromosome'] +'_'+ B2_3['Start'].map(str) +'_'+ B2_3['End'].map(str) +'_'+ B2_3['Name']

M17_1['ID'] = M17_1['Chromosome'] +'_'+ M17_1['Start'].map(str) +'_'+ M17_1['End'].map(str) +'_'+ M17_1['Name']
M17_2['ID'] = M17_2['Chromosome'] +'_'+ M17_2['Start'].map(str) +'_'+ M17_2['End'].map(str) +'_'+ M17_2['Name']
M17_3['ID'] = M17_3['Chromosome'] +'_'+ M17_3['Start'].map(str) +'_'+ M17_3['End'].map(str) +'_'+ M17_3['Name']
'''
a3_1 = A3_1_ff['ID'].tolist()
b2_1 = B2_1_ff['ID'].tolist()
m17_1 = M17_1_ff['ID'].tolist()

a3_2 = A3_2_ff['ID'].tolist()
b2_2 = B2_2_ff['ID'].tolist()
m17_2 = M17_2_ff['ID'].tolist()

a3_3 = A3_3_ff['ID'].tolist()
b2_3 = B2_3_ff['ID'].tolist()
m17_3 = M17_3_ff['ID'].tolist()
'''

a3_1 = A3_1['ID'].tolist()
b2_1 = B2_1['ID'].tolist()
m17_1 = M17_1['ID'].tolist()

a3_2 = A3_2['ID'].tolist()
b2_2 = B2_2['ID'].tolist()
m17_2 = M17_2['ID'].tolist()

a3_3 = A3_3['ID'].tolist()
b2_3 = B2_3['ID'].tolist()
m17_3 = M17_3['ID'].tolist()

# calculate half life using the new equation

#A3_1_ff['half_life_0_1'] = -np.log(2) / np.log(1-((A3_1_ff['ConversionRate_1'] - A3_1_ff['ConversionRate_0']) / (A3_1_ff['CoverageOnTs_1'] - A3_1_ff['CoverageOnTs_0'])))
#A3_1_ff.to_csv(results_dir+'test.csv')

# fit exponential curve to each fragment and calculate half life, keep only positive half lives
# EDIT 11/11/2021: fit exponential curve to each GROUP of fragments, new code is marked

#x = np.array([0,1,2,4]).astype(float)
x = np.array([0,1,2,4,0,1,2,4,0,1,2,4]).astype(float)

def exp_fit(x, a, b, c, d):
    return a * np.exp(b * (x + c)) + d

def sigmoid_fit(x, L, x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)

'''
loop = [a3_1,
        b2_1,
        m17_1,
        a3_2,
        b2_2,
        m17_2,
        a3_3,
        b2_3,
        m17_3
        ]

loop = [a3_1,
        b2_1,
        m17_1
        ]

for fuck in loop:
    if fuck == a3_1:
        shit = 'A3_1'
        dick = A3_1_ff
    elif fuck == b2_1:
        shit = 'B2_1'
        dick = B2_1_ff
    elif fuck == m17_1:
        shit = 'M17_1'
        dick = M17_1_ff
    elif fuck == a3_2:
        shit = 'A3_2'
        dick = A3_2_ff
    elif fuck == b2_2:
        shit = 'B2_2'
        dick = B2_2_ff
    elif fuck == m17_2:
        shit = 'M17_2'
        dick = M17_2_ff
    elif fuck == a3_3:
        shit = 'A3_3'
        dick = A3_3_ff
    elif fuck == b2_3:
        shit = 'B2_3'
        dick = B2_3_ff
    elif fuck == m17_3:
        shit = 'M17_3'
        dick = M17_3_ff
    if fuck == a3_1:
        shit1 = 'A3_1'
        dick1 = A3_1
        shit2 = 'A3_2'
        dick2 = A3_2
        shit3 = 'A3_3'
        dick3 = A3_3
    elif fuck == b2_1:
        shit1 = 'B2_1'
        dick1 = B2_1
        shit2 = 'B2_2'
        dick2 = B2_2
        shit3 = 'B2_3'
        dick3 = B2_3
    elif fuck == m17_1:
        shit1 = 'M17_1'
        dick1 = M17_1
        shit2 = 'M17_2'
        dick2 = M17_2
        shit3 = 'M17_3'
        dick3 = M17_3'''

i=0
half_lifes = {}
k_exp = {}
k_sig = {}
w1 = a3_1[0] # for testing only
z_thresh = 1.05
for w1 in a3_1:
    i += 1
    #dick_0 = float(dick.loc[dick['ID'] == w1]['ConversionRate_0_n'])
    #dick_1 = float(dick.loc[dick['ID'] == w1]['ConversionRate_1_n'])
    #dick_2 = float(dick.loc[dick['ID'] == w1]['ConversionRate_2_n'])
    #dick_4 = float(dick.loc[dick['ID'] == w1]['ConversionRate_4_n'])
    
    a3_dick_10 = float(A3_1.loc[A3_1['ID'] == w1]['ConversionRate_0'])
    a3_dick_11 = float(A3_1.loc[A3_1['ID'] == w1]['ConversionRate_1'])
    a3_dick_12 = float(A3_1.loc[A3_1['ID'] == w1]['ConversionRate_2'])
    a3_dick_14 = float(A3_1.loc[A3_1['ID'] == w1]['ConversionRate_4'])
    a3_dick_20 = float(A3_2.loc[A3_2['ID'] == w1]['ConversionRate_0'])
    a3_dick_21 = float(A3_2.loc[A3_2['ID'] == w1]['ConversionRate_1'])
    a3_dick_22 = float(A3_2.loc[A3_2['ID'] == w1]['ConversionRate_2'])
    a3_dick_24 = float(A3_2.loc[A3_2['ID'] == w1]['ConversionRate_4'])
    a3_dick_30 = float(A3_3.loc[A3_3['ID'] == w1]['ConversionRate_0'])
    a3_dick_31 = float(A3_3.loc[A3_3['ID'] == w1]['ConversionRate_1'])
    a3_dick_32 = float(A3_3.loc[A3_3['ID'] == w1]['ConversionRate_2'])
    a3_dick_34 = float(A3_3.loc[A3_3['ID'] == w1]['ConversionRate_4'])
    
    b2_dick_10 = float(B2_1.loc[B2_1['ID'] == w1]['ConversionRate_0'])
    b2_dick_11 = float(B2_1.loc[B2_1['ID'] == w1]['ConversionRate_1'])
    b2_dick_12 = float(B2_1.loc[B2_1['ID'] == w1]['ConversionRate_2'])
    b2_dick_14 = float(B2_1.loc[B2_1['ID'] == w1]['ConversionRate_4'])
    b2_dick_20 = float(B2_2.loc[B2_2['ID'] == w1]['ConversionRate_0'])
    b2_dick_21 = float(B2_2.loc[B2_2['ID'] == w1]['ConversionRate_1'])
    b2_dick_22 = float(B2_2.loc[B2_2['ID'] == w1]['ConversionRate_2'])
    b2_dick_24 = float(B2_2.loc[B2_2['ID'] == w1]['ConversionRate_4'])
    b2_dick_30 = float(B2_3.loc[B2_3['ID'] == w1]['ConversionRate_0'])
    b2_dick_31 = float(B2_3.loc[B2_3['ID'] == w1]['ConversionRate_1'])
    b2_dick_32 = float(B2_3.loc[B2_3['ID'] == w1]['ConversionRate_2'])
    b2_dick_34 = float(B2_3.loc[B2_3['ID'] == w1]['ConversionRate_4'])
    
    m17_dick_10 = float(M17_1.loc[M17_1['ID'] == w1]['ConversionRate_0'])
    m17_dick_11 = float(M17_1.loc[M17_1['ID'] == w1]['ConversionRate_1'])
    m17_dick_12 = float(M17_1.loc[M17_1['ID'] == w1]['ConversionRate_2'])
    m17_dick_14 = float(M17_1.loc[M17_1['ID'] == w1]['ConversionRate_4'])
    m17_dick_20 = float(M17_2.loc[M17_2['ID'] == w1]['ConversionRate_0'])
    m17_dick_21 = float(M17_2.loc[M17_2['ID'] == w1]['ConversionRate_1'])
    m17_dick_22 = float(M17_2.loc[M17_2['ID'] == w1]['ConversionRate_2'])
    m17_dick_24 = float(M17_2.loc[M17_2['ID'] == w1]['ConversionRate_4'])
    m17_dick_30 = float(M17_3.loc[M17_3['ID'] == w1]['ConversionRate_0'])
    m17_dick_31 = float(M17_3.loc[M17_3['ID'] == w1]['ConversionRate_1'])
    m17_dick_32 = float(M17_3.loc[M17_3['ID'] == w1]['ConversionRate_2'])
    m17_dick_34 = float(M17_3.loc[M17_3['ID'] == w1]['ConversionRate_4'])
    
    try:
        #print('trying a3')
        #y = np.array([dick_0, dick_1, dick_2, dick_4]).astype(float)
        y_a3 = np.array([a3_dick_10,
                      a3_dick_11,
                      a3_dick_12,
                      a3_dick_14,
                      a3_dick_20,
                      a3_dick_21,
                      a3_dick_22,
                      a3_dick_24,
                      a3_dick_30,
                      a3_dick_31,
                      a3_dick_32,
                      a3_dick_34
                      ]).astype(float)

        y2 = y_a3
        y2[y2==0] = np.nan
        
        ############
        # added code to run z score outlier test on each timepoint
        ############
        
        a3_0_avg = (y2[0] + y2[4] + y2[8]) / 3
        if ~np.isnan(a3_0_avg):
            a3_0_std = np.std([y2[0], y2[4], y2[8]])
            a3_10_z = (y2[0] - a3_0_avg) / a3_0_std
            a3_20_z = (y2[4] - a3_0_avg) / a3_0_std
            a3_30_z = (y2[8] - a3_0_avg) / a3_0_std
            
            cnt = 0
            if abs(a3_10_z) > z_thresh: cnt += 1
            if abs(a3_20_z) > z_thresh: cnt += 1
            if abs(a3_30_z) > z_thresh: cnt += 1
            
            if cnt < 2:
                if abs(a3_10_z) > z_thresh: y2[0] = np.nan
                if abs(a3_20_z) > z_thresh: y2[4] = np.nan
                if abs(a3_30_z) > z_thresh: y2[8] = np.nan
        
        a3_1_avg = (y2[1] + y2[5] + y2[9]) / 3
        if ~np.isnan(a3_1_avg):
            a3_1_std = np.std([y2[1], y2[5], y2[9]])
            a3_11_z = (y2[1] - a3_1_avg) / a3_1_std
            a3_21_z = (y2[5] - a3_1_avg) / a3_1_std
            a3_31_z = (y2[9] - a3_1_avg) / a3_1_std
            
            cnt = 0
            if abs(a3_11_z) > z_thresh: cnt += 1
            if abs(a3_21_z) > z_thresh: cnt += 1
            if abs(a3_31_z) > z_thresh: cnt += 1
            
            if cnt < 2:
                if abs(a3_11_z) > z_thresh: y2[1] = np.nan
                if abs(a3_21_z) > z_thresh: y2[5] = np.nan
                if abs(a3_31_z) > z_thresh: y2[9] = np.nan
        
        a3_2_avg = (y2[2] + y2[6] + y2[10]) / 3
        if ~np.isnan(a3_2_avg):
            a3_2_std = np.std([y2[2], y2[6], y2[10]])
            a3_12_z = (y2[2] - a3_2_avg) / a3_2_std
            a3_22_z = (y2[6] - a3_2_avg) / a3_2_std
            a3_32_z = (y2[10] - a3_2_avg) / a3_2_std
            
            cnt = 0
            if abs(a3_12_z) > z_thresh: cnt += 1
            if abs(a3_22_z) > z_thresh: cnt += 1
            if abs(a3_32_z) > z_thresh: cnt += 1
            
            if cnt < 2:
                if abs(a3_12_z) > z_thresh: y2[2] = np.nan
                if abs(a3_22_z) > z_thresh: y2[6] = np.nan
                if abs(a3_32_z) > z_thresh: y2[10] = np.nan
        
        a3_4_avg = (y2[3] + y2[7] + y2[11]) / 3
        if ~np.isnan(a3_4_avg):
            a3_4_std = np.std([y2[3], y2[7], y2[11]])
            a3_14_z = (y2[3] - a3_4_avg) / a3_4_std
            a3_24_z = (y2[7] - a3_4_avg) / a3_4_std
            a3_34_z = (y2[11] - a3_4_avg) / a3_4_std
            
            cnt = 0
            if abs(a3_14_z) > z_thresh: cnt += 1
            if abs(a3_24_z) > z_thresh: cnt += 1
            if abs(a3_34_z) > z_thresh: cnt += 1
            
            if cnt < 2:
                if abs(a3_14_z) > z_thresh: y2[3] = np.nan
                if abs(a3_24_z) > z_thresh: y2[7] = np.nan
                if abs(a3_34_z) > z_thresh: y2[11] = np.nan
        
        ############

        indexy = ~(np.isnan(y2))
        y1 = y2[indexy]
        index = ~(np.isnan(y2))
        x1 = x[index]

        popt_exp_a3, pcov_exp_a3 = scipy.optimize.curve_fit(exp_fit,
                                                            x1,
                                                            y1,
                                                            check_finite = False,
                                                            method='trf',
                                                            maxfev=10000,
                                                            bounds = ([-np.inf, 0, -np.inf, -np.inf],[np.inf, 10, np.inf, np.inf])
                                                            )

        hl = np.log(2 / popt_exp_a3[1])
        
        if hl > 0:
            half_lifes[w1] = hl
        
        residuals = y1 - exp_fit(x1, *popt_exp_a3)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y1-np.mean(y1))**2)
        r_squared_exp_a3 = 1 - (ss_res / ss_tot)
        
        p0 = [max(y1),np.median(y1),1,min(y1)] # this is an mandatory initial guess

        popt_sig_a3, pcov_sig_a3 = scipy.optimize.curve_fit(sigmoid_fit,
                                                            x1,
                                                            y1,
                                                            p0,
                                                            check_finite = False,
                                                            method='trf',
                                                            maxfev=10000,
                                                            bounds = ([-np.inf, -np.inf, 0, -np.inf],[np.inf, np.inf, np.inf, np.inf])
                                                            )
        
        residuals = y1 - sigmoid_fit(x1, *popt_sig_a3)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y1-np.mean(y1))**2)
        r_squared_sig_a3 = 1 - (ss_res / ss_tot)
        
        k_exp[w1] = (popt_exp_a3,r_squared_exp_a3)
        
    except:
        popt_sig_a3 = np.array([np.nan,np.nan, np.nan, np.nan])
        r_squared_sig_a3 = np.nan
        
        popt_exp_a3 = np.array([np.nan,np.nan,np.nan,np.nan])
        r_squared_exp_a3 = np.nan
        pass
    
    try:
        #print('trying b2')
        #y = np.array([dick_0, dick_1, dick_2, dick_4]).astype(float)
        y_b2 = np.array([b2_dick_10,
                      b2_dick_11,
                      b2_dick_12,
                      b2_dick_14,
                      b2_dick_20,
                      b2_dick_21,
                      b2_dick_22,
                      b2_dick_24,
                      b2_dick_30,
                      b2_dick_31,
                      b2_dick_32,
                      b2_dick_34
                      ]).astype(float)

        y2 = y_b2
        y2[y2==0] = np.nan
        
        ############
        # added code to run z score outlier test on each timepoint
        ############
        
        b2_0_avg = (y2[0] + y2[4] + y2[8]) / 3
        if ~np.isnan(b2_0_avg):
            b2_0_std = np.std([y2[0], y2[4], y2[8]])
            b2_10_z = (y2[0] - b2_0_avg) / b2_0_std
            b2_20_z = (y2[4] - b2_0_avg) / b2_0_std
            b2_30_z = (y2[8] - b2_0_avg) / b2_0_std
            
            cnt = 0
            if abs(b2_10_z) > z_thresh: cnt += 1
            if abs(b2_20_z) > z_thresh: cnt += 1
            if abs(b2_30_z) > z_thresh: cnt += 1
            
            if cnt < 2:
                if abs(b2_10_z) > z_thresh: y2[0] = np.nan
                if abs(b2_20_z) > z_thresh: y2[4] = np.nan
                if abs(b2_30_z) > z_thresh: y2[8] = np.nan
        
        b2_1_avg = (y2[1] + y2[5] + y2[9]) / 3
        if ~np.isnan(b2_1_avg):
            b2_1_std = np.std([y2[1], y2[5], y2[9]])
            b2_11_z = (y2[1] - b2_1_avg) / b2_1_std
            b2_21_z = (y2[5] - b2_1_avg) / b2_1_std
            b2_31_z = (y2[9] - b2_1_avg) / b2_1_std
            
            cnt = 0
            if abs(b2_11_z) > z_thresh: cnt += 1
            if abs(b2_21_z) > z_thresh: cnt += 1
            if abs(b2_31_z) > z_thresh: cnt += 1
            
            if cnt < 2:
                if abs(b2_11_z) > z_thresh: y2[1] = np.nan
                if abs(b2_21_z) > z_thresh: y2[5] = np.nan
                if abs(b2_31_z) > z_thresh: y2[9] = np.nan
        
        b2_2_avg = (y2[2] + y2[6] + y2[10]) / 3
        if ~np.isnan(b2_2_avg):
            b2_2_std = np.std([y2[2], y2[6], y2[10]])
            b2_12_z = (y2[2] - b2_2_avg) / b2_2_std
            b2_22_z = (y2[6] - b2_2_avg) / b2_2_std
            b2_32_z = (y2[10] - b2_2_avg) / b2_2_std
            if abs(b2_12_z) > z_thresh: y2[2] = np.nan
            if abs(b2_22_z) > z_thresh: y2[6] = np.nan
            if abs(b2_32_z) > z_thresh: y2[10] = np.nan
        
        b2_4_avg = (y2[3] + y2[7] + y2[11]) / 3
        if ~np.isnan(b2_4_avg):
            b2_4_std = np.std([y2[3], y2[7], y2[11]])
            b2_14_z = (y2[3] - b2_4_avg) / b2_4_std
            b2_24_z = (y2[7] - b2_4_avg) / b2_4_std
            b2_34_z = (y2[11] - b2_4_avg) / b2_4_std
            
            cnt = 0
            if abs(b2_12_z) > z_thresh: cnt += 1
            if abs(b2_22_z) > z_thresh: cnt += 1
            if abs(b2_32_z) > z_thresh: cnt += 1
            
            if cnt < 2:
                if abs(b2_14_z) > z_thresh: y2[3] = np.nan
                if abs(b2_24_z) > z_thresh: y2[7] = np.nan
                if abs(b2_34_z) > z_thresh: y2[11] = np.nan
        
        ############        
        
        indexy = ~(np.isnan(y2))
        y1 = y2[indexy]
        index = ~(np.isnan(y2))
        x1 = x[index]
    
        popt_exp_b2, pcov_exp_b2 = scipy.optimize.curve_fit(exp_fit,
                                                            x1,
                                                            y1,
                                                            check_finite = False,
                                                            method='trf',
                                                            maxfev=10000,
                                                            bounds = ([-np.inf, 0, -np.inf, -np.inf],[np.inf, 10, np.inf, np.inf])
                                                            )
        
        p0 = [max(y1),np.median(y1),1,min(y1)] # this is an mandatory initial guess

        popt_sig_b2, pcov_sig_b2 = scipy.optimize.curve_fit(sigmoid_fit,
                                                            x1,
                                                            y1,
                                                            p0,
                                                            check_finite = False,
                                                            method='trf',
                                                            maxfev=10000,
                                                            bounds = ([-np.inf, -np.inf, 0, -np.inf],[np.inf, np.inf, np.inf,np.inf])
                                                            )

        hl = np.log(2 / popt_exp_b2[1])

        if hl > 0:
            half_lifes[w1] = hl
        
        residuals = y1 - exp_fit(x1, *popt_exp_b2)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y1-np.mean(y1))**2)
        r_squared_exp_b2 = 1 - (ss_res / ss_tot)
        
        residuals = y1 - sigmoid_fit(x1, *popt_sig_b2)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y1-np.mean(y1))**2)
        r_squared_sig_b2 = 1 - (ss_res / ss_tot)
        
        k_exp[w1] = (popt_exp_b2,r_squared_exp_b2)
        
    except:
        popt_sig_b2 = np.array([np.nan,np.nan,np.nan,np.nan])
        r_squared_sig_b2 = np.nan
        
        popt_exp_b2 = np.array([np.nan,np.nan,np.nan,np.nan])
        r_squared_exp_b2 = np.nan
        pass
    
    try:
        #print('trying m17')
        #y = np.array([dick_0, dick_1, dick_2, dick_4]).astype(float)
        y_m17 = np.array([m17_dick_10,
                      m17_dick_11,
                      m17_dick_12,
                      m17_dick_14,
                      m17_dick_20,
                      m17_dick_21,
                      m17_dick_22,
                      m17_dick_24,
                      m17_dick_30,
                      m17_dick_31,
                      m17_dick_32,
                      m17_dick_34
                      ]).astype(float)

        y2 = y_m17
        y2[y2==0] = np.nan
        
        ############
        # added code to run z score outlier test on each timepoint
        ############
        
        m17_0_avg = (y2[0] + y2[4] + y2[8]) / 3
        if ~np.isnan(m17_0_avg):
            m17_0_std = np.std([y2[0], y2[4], y2[8]])
            m17_10_z = (y2[0] - m17_0_avg) / m17_0_std
            m17_20_z = (y2[4] - m17_0_avg) / m17_0_std
            m17_30_z = (y2[8] - m17_0_avg) / m17_0_std
            
            cnt = 0
            if abs(m17_10_z) > z_thresh: cnt += 1
            if abs(m17_20_z) > z_thresh: cnt += 1
            if abs(m17_30_z) > z_thresh: cnt += 1
            
            if cnt < 2:
                if abs(m17_10_z) > z_thresh: y2[0] = np.nan
                if abs(m17_20_z) > z_thresh: y2[4] = np.nan
                if abs(m17_30_z) > z_thresh: y2[8] = np.nan
        
        m17_1_avg = (y2[1] + y2[5] + y2[9]) / 3
        if ~np.isnan(m17_1_avg):
            m17_1_std = np.std([y2[1], y2[5], y2[9]])
            m17_11_z = (y2[1] - m17_1_avg) / m17_1_std
            m17_21_z = (y2[5] - m17_1_avg) / m17_1_std
            m17_31_z = (y2[9] - m17_1_avg) / m17_1_std
            
            cnt = 0
            if abs(m17_11_z) > z_thresh: cnt += 1
            if abs(m17_21_z) > z_thresh: cnt += 1
            if abs(m17_31_z) > z_thresh: cnt += 1
            
            if cnt < 2:
                if abs(m17_11_z) > z_thresh: y2[1] = np.nan
                if abs(m17_21_z) > z_thresh: y2[5] = np.nan
                if abs(m17_31_z) > z_thresh: y2[9] = np.nan
        
        m17_2_avg = (y2[2] + y2[6] + y2[10]) / 3
        if ~np.isnan(m17_2_avg):
            m17_2_std = np.std([y2[2], y2[6], y2[10]])
            m17_12_z = (y2[2] - m17_2_avg) / m17_2_std
            m17_22_z = (y2[6] - m17_2_avg) / m17_2_std
            m17_32_z = (y2[10] - m17_2_avg) / m17_2_std
            
            cnt = 0
            if abs(m17_12_z) > z_thresh: cnt += 1
            if abs(m17_22_z) > z_thresh: cnt += 1
            if abs(m17_32_z) > z_thresh: cnt += 1
            
            if cnt < 2:
                if abs(m17_12_z) > z_thresh: y2[2] = np.nan
                if abs(m17_22_z) > z_thresh: y2[6] = np.nan
                if abs(m17_32_z) > z_thresh: y2[10] = np.nan
        
        m17_4_avg = (y2[3] + y2[7] + y2[11]) / 3
        if ~np.isnan(m17_4_avg):
            m17_4_std = np.std([y2[3], y2[7], y2[11]])
            m17_14_z = (y2[3] - m17_4_avg) / m17_4_std
            m17_24_z = (y2[7] - m17_4_avg) / m17_4_std
            m17_34_z = (y2[11] - m17_4_avg) / m17_4_std
            
            cnt = 0
            if abs(m17_14_z) > z_thresh: cnt += 1
            if abs(m17_24_z) > z_thresh: cnt += 1
            if abs(m17_34_z) > z_thresh: cnt += 1
            
            if cnt < 2:
                if abs(m17_14_z) > z_thresh: y2[3] = np.nan
                if abs(m17_24_z) > z_thresh: y2[7] = np.nan
                if abs(m17_34_z) > z_thresh: y2[11] = np.nan
        
        ############
        
        indexy = ~(np.isnan(y2))
        y1 = y2[indexy]
        index = ~(np.isnan(y2))
        x1 = x[index]
    
        popt_exp_m17, pcov_exp_m17 = scipy.optimize.curve_fit(exp_fit,
                                                              x1,
                                                              y1,
                                                              check_finite = False,
                                                              method='trf',
                                                              maxfev=10000,
                                                              bounds = ([-np.inf, 0, -np.inf, -np.inf],[np.inf, 10, np.inf, np.inf])
                                                              )
        
        p0 = [max(y1),np.median(y1),1,min(y1)] # this is an mandatory initial guess

        popt_sig_m17, pcov_sig_m17 = scipy.optimize.curve_fit(sigmoid_fit,
                                                              x1,
                                                              y1,
                                                              p0,
                                                              check_finite = False,
                                                              method='trf',
                                                              maxfev=10000,
                                                              bounds = ([-np.inf, -np.inf, 0, -np.inf],[np.inf, np.inf,np.inf,np.inf])
                                                              )

        hl = np.log(2 / popt_exp_m17[1])

        if hl > 0:
            half_lifes[w1] = hl
        
        residuals = y1 - exp_fit(x1, *popt_exp_m17)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y1-np.mean(y1))**2)
        r_squared_exp_m17 = 1 - (ss_res / ss_tot)
        
        residuals = y1 - sigmoid_fit(x1, *popt_sig_m17)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y1-np.mean(y1))**2)
        r_squared_sig_m17 = 1 - (ss_res / ss_tot)
        
        k_exp[w1] = (popt_exp_m17,r_squared_exp_m17)
    except:
        popt_sig_m17 = np.array([np.nan,np.nan,np.nan,np.nan])
        r_squared_sig_m17 = np.nan
        
        popt_exp_m17 = np.array([np.nan,np.nan,np.nan,np.nan])
        r_squared_exp_m17 = np.nan
        pass
    
    k_sig[w1] = (popt_sig_a3[0],
                 popt_sig_a3[1],
                 popt_sig_a3[2],
                 popt_sig_a3[3],
                 r_squared_sig_a3,
                 popt_exp_a3[0],
                 popt_exp_a3[1],
                 popt_exp_a3[2],
                 popt_exp_a3[3],
                 r_squared_exp_a3,
                 a3_dick_10,
                 a3_dick_11,
                 a3_dick_12,
                 a3_dick_14,
                 a3_dick_20,
                 a3_dick_21,
                 a3_dick_22,
                 a3_dick_24,
                 a3_dick_30,
                 a3_dick_31,
                 a3_dick_32,
                 a3_dick_34,
                 popt_sig_b2[0],
                 popt_sig_b2[1],
                 popt_sig_b2[2],
                 popt_sig_b2[3],
                 r_squared_sig_b2,
                 popt_exp_b2[0],
                 popt_exp_b2[1],
                 popt_exp_b2[2],
                 popt_exp_b2[3],
                 r_squared_exp_b2,
                 b2_dick_10,
                 b2_dick_11,
                 b2_dick_12,
                 b2_dick_14,
                 b2_dick_20,
                 b2_dick_21,
                 b2_dick_22,
                 b2_dick_24,
                 b2_dick_30,
                 b2_dick_31,
                 b2_dick_32,
                 b2_dick_34,
                 popt_sig_m17[0],
                 popt_sig_m17[1],
                 popt_sig_m17[2],
                 popt_sig_m17[3],
                 r_squared_sig_m17,
                 popt_exp_m17[0],
                 popt_exp_m17[1],
                 popt_exp_m17[2],
                 popt_exp_m17[3],
                 r_squared_exp_m17,
                 m17_dick_10,
                 m17_dick_11,
                 m17_dick_12,
                 m17_dick_14,
                 m17_dick_20,
                 m17_dick_21,
                 m17_dick_22,
                 m17_dick_24,
                 m17_dick_30,
                 m17_dick_31,
                 m17_dick_32,
                 m17_dick_34)
    #print(i,len(a3_1))
    print("{0:.3%}".format(i/len(a3_1)))
#waterfall['{}_fit'.format(shit1)] = i-j

#halfs_df = pd.DataFrame.from_dict(half_lifes, orient='index')
#halfs_df.reset_index(inplace=True)
#halfs_df.columns = ['Name','half_life']
#halfs_df.loc[halfs_df['half_life'] > 0].describe()
#halfs_df.to_csv(results_dir + 'half_lifes_{}.csv'.format(shit), index=False)

#k_exp_df = pd.DataFrame.from_dict(k_exp, orient='index')
#k_exp_df.reset_index(inplace=True)
#k_exp_df.columns = ['Name','exp_a', 'exp_b', 'exp_c', 'exp_d', 'r2_exp']

k_sig_df = pd.DataFrame.from_dict(k_sig, orient='index')
k_sig_df.reset_index(inplace=True)
k_sig_df.columns = ['Name',
                    'A3_sig_L',
                    'A3_sig_x0',
                    'A3_sig_k',
                    'A3_sig_b',
                    'A3_r2_sig',
                    'A3_exp_a',
                    'A3_exp_b',
                    'A3_exp_c',
                    'A3_exp_d',
                    'A3_r2_exp',
                    'A3_rate_10',
                    'A3_rate_11',
                    'A3_rate_12',
                    'A3_rate_14',
                    'A3_rate_20',
                    'A3_rate_21',
                    'A3_rate_22',
                    'A3_rate_24',
                    'A3_rate_30',
                    'A3_rate_31',
                    'A3_rate_32',
                    'A3_rate_34',
                    'B2_sig_L',
                    'B2_sig_x0',
                    'B2_sig_k',
                    'B2_sig_b',
                    'B2_r2_sig',
                    'B2_exp_a',
                    'B2_exp_b',
                    'B2_exp_c',
                    'B2_exp_d',
                    'B2_r2_exp',
                    'B2_rate_10',
                    'B2_rate_11',
                    'B2_rate_12',
                    'B2_rate_14',
                    'B2_rate_20',
                    'B2_rate_21',
                    'B2_rate_22',
                    'B2_rate_24',
                    'B2_rate_30',
                    'B2_rate_31',
                    'B2_rate_32',
                    'B2_rate_34',
                    'M17_sig_L',
                    'M17_sig_x0',
                    'M17_sig_k',
                    'M17_sig_b',
                    'M17_r2_sig',
                    'M17_exp_a',
                    'M17_exp_b',
                    'M17_exp_c',
                    'M17_exp_d',
                    'M17_r2_exp',
                    'M17_rate_10',
                    'M17_rate_11',
                    'M17_rate_12',
                    'M17_rate_14',
                    'M17_rate_20',
                    'M17_rate_21',
                    'M17_rate_22',
                    'M17_rate_24',
                    'M17_rate_30',
                    'M17_rate_31',
                    'M17_rate_32',
                    'M17_rate_34',
                    ]
#k_sig_df.to_csv(results_dir + 'k_{}_AVG.csv'.format('all'), index=False)

k_sig_df.to_csv(data_dir2 + 'k_{}_AVG_trf_bounded_left_only.csv'.format('all'), index=False)


#k_df = k_exp_df.merge(k_sig_df, on='Name', how='outer')
#k_df.to_csv(results_dir + 'k_{}.csv'.format(shit), index=False)

#pd.DataFrame(waterfall.items()).to_csv(results_dir + "waterfall_{}_AVG.csv".format('all'))

k_sig_df = pd.read_csv(results_dir + 'k_{}_AVG.csv'.format('all'))
k_sig_df['ID'] = k_sig_df['Name'].str.split('_').str[3]

ivanas_csv = pd.read_csv(results_dir + 'transcription rates 11 15 2021.csv')

new_db = pd.read_csv(results_dir + 'mart_export.txt')
new_db['gene_stable_id'] = new_db['Gene stable ID version'].str.split('.').str[0]
new_db_gene = new_db[['gene_stable_id','Gene name']].drop_duplicates()

k_final_db = k_sig_df.merge(new_db_gene, left_on='ID', right_on='gene_stable_id', how='left')
k_final_db.to_csv(results_dir+'k_{}_AVG_db.csv'.format('all'), index=False)

ivanas_csv_db = ivanas_csv.merge(new_db_gene, left_on='Name', right_on='gene_stable_id', how='left')
ivanas_csv_db.to_csv(results_dir+'ivanas_csv_db.csv'.format('all'), index=False)

###############################################################################
# Tool for plotting selected gene fit lines

data = pd.read_csv(data_dir2+'k_all_AVG_trf_bounded_left_only.csv')
data['gene'] = data['Name'].str.split('_').str[0]
x_p = np.linspace(0, 4, num=50)
x = np.array([0,1,2,4,0,1,2,4,0,1,2,4]).astype(float)

def exp_fit(x, a, b, c, d):
    return a * np.exp(b * (x + c)) + d

def sigmoid_fit(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0)))+b
    return (y)

# grab the gene with the highest r2 in each bin of k values
    
bins = np.linspace(0, 20, num=41)
data = data.loc[data['A3_r2_sig'] >= 0.6]

genes = []

for bin1 in bins:
    filtered = data.loc[((data['A3_sig_k'] > bin1) & (data['A3_sig_k'] < bin1 + 0.5))]
    genes.append(filtered.sort_values(by='A3_r2_sig', ascending=True).iloc[0]['gene'])

for gene in genes:
    d = data[data['gene']==gene].to_dict('records')[0]
    
    A3_y = np.array([d['A3_rate_10'],
                     d['A3_rate_11'],
                     d['A3_rate_12'],
                     d['A3_rate_14'],
                     d['A3_rate_20'],
                     d['A3_rate_21'],
                     d['A3_rate_22'],
                     d['A3_rate_24'],
                     d['A3_rate_30'],
                     d['A3_rate_31'],
                     d['A3_rate_32'],
                     d['A3_rate_34'],])
    
    A3_popt_sig = np.array([d['A3_sig_L'],
                            d['A3_sig_x0'],
                            d['A3_sig_k'],
                            d['A3_sig_b'],])
    
    B2_y = np.array([d['B2_rate_10'],
                     d['B2_rate_11'],
                     d['B2_rate_12'],
                     d['B2_rate_14'],
                     d['B2_rate_20'],
                     d['B2_rate_21'],
                     d['B2_rate_22'],
                     d['B2_rate_24'],
                     d['B2_rate_30'],
                     d['B2_rate_31'],
                     d['B2_rate_32'],
                     d['B2_rate_34'],])
    
    B2_popt_sig = np.array([d['B2_sig_L'],
                            d['B2_sig_x0'],
                            d['B2_sig_k'],
                            d['B2_sig_b'],])
    
    M17_y = np.array([d['M17_rate_10'],
                     d['M17_rate_11'],
                     d['M17_rate_12'],
                     d['M17_rate_14'],
                     d['M17_rate_20'],
                     d['M17_rate_21'],
                     d['M17_rate_22'],
                     d['M17_rate_24'],
                     d['M17_rate_30'],
                     d['M17_rate_31'],
                     d['M17_rate_32'],
                     d['M17_rate_34'],])
    
    M17_popt_sig = np.array([d['M17_sig_L'],
                            d['M17_sig_x0'],
                            d['M17_sig_k'],
                            d['M17_sig_b'],])
    
    plt.plot(x, A3_y, 'ko', label="A3", c='r')
    plt.plot(x_p, sigmoid_fit(x_p, *A3_popt_sig), 'r-', label="Fitted Curve")
    plt.title(gene)
    #plt.plot(x, B2_y, 'ko', label="B2", c='g')
    #plt.plot(x_p, sigmoid_fit(x_p, *B2_popt_sig), 'g-', label="Fitted Curve")
    #plt.plot(x, M17_y, 'ko', label="M17", c='b')
    #plt.plot(x_p, sigmoid_fit(x_p, *M17_popt_sig), 'b-', label="Fitted Curve")
    plt.legend()
    plt.savefig(results_dir + str(d['A3_sig_k']) + '_' + str(d['A3_r2_sig']) + '.jpg')
    plt.clf()

################################################################################
# new post-processing analysis

r2_thresh = 0.6
var_thresh = 0.25
    
k_A3_1 = pd.read_csv(results_dir + 'k_A3_1.csv')
k_B2_1 = pd.read_csv(results_dir + 'k_B2_1.csv')
k_M17_1 = pd.read_csv(results_dir + 'k_M17_1.csv')

k_A3_2 = pd.read_csv(results_dir + 'k_A3_2.csv')
k_B2_2 = pd.read_csv(results_dir + 'k_B2_2.csv')
k_M17_2 = pd.read_csv(results_dir + 'k_M17_2.csv')

k_A3_3 = pd.read_csv(results_dir + 'k_A3_3.csv')
k_B2_3 = pd.read_csv(results_dir + 'k_B2_3.csv')
k_M17_3 = pd.read_csv(results_dir + 'k_M17_3.csv')

# some strings in my r2 column for some reason, remove them

k_A3_1 = k_A3_1[pd.to_numeric(k_A3_1['r2_sig'], errors='coerce').notnull()]
k_A3_2 = k_A3_2[pd.to_numeric(k_A3_2['r2_sig'], errors='coerce').notnull()]
k_A3_3 = k_A3_3[pd.to_numeric(k_A3_3['r2_sig'], errors='coerce').notnull()]

k_B2_1 = k_B2_1[pd.to_numeric(k_B2_1['r2_sig'], errors='coerce').notnull()]
k_B2_2 = k_B2_2[pd.to_numeric(k_B2_2['r2_sig'], errors='coerce').notnull()]
k_B2_3 = k_B2_3[pd.to_numeric(k_B2_3['r2_sig'], errors='coerce').notnull()]

k_M17_1 = k_M17_1[pd.to_numeric(k_M17_1['r2_sig'], errors='coerce').notnull()]
k_M17_2 = k_M17_2[pd.to_numeric(k_M17_2['r2_sig'], errors='coerce').notnull()]
k_M17_3 = k_M17_3[pd.to_numeric(k_M17_3['r2_sig'], errors='coerce').notnull()]

# filter by r2 threshold value

k_A3_1 = k_A3_1.loc[k_A3_1['r2_sig'].astype(float) > r2_thresh]
k_A3_2 = k_A3_2.loc[k_A3_2['r2_sig'].astype(float) > r2_thresh]
k_A3_3 = k_A3_3.loc[k_A3_3['r2_sig'].astype(float) > r2_thresh]

k_B2_1 = k_B2_1.loc[k_B2_1['r2_sig'].astype(float) > r2_thresh]
k_B2_2 = k_B2_2.loc[k_B2_2['r2_sig'].astype(float) > r2_thresh]
k_B2_3 = k_B2_3.loc[k_B2_3['r2_sig'].astype(float) > r2_thresh]

k_M17_1 = k_M17_1.loc[k_M17_1['r2_sig'].astype(float) > r2_thresh]
k_M17_2 = k_M17_2.loc[k_M17_2['r2_sig'].astype(float) > r2_thresh]
k_M17_3 = k_M17_3.loc[k_M17_3['r2_sig'].astype(float) > r2_thresh]

# filter by k > 0
k_A3_1 = k_A3_1.loc[k_A3_1['sig_k'].astype(float) > 0]
k_A3_2 = k_A3_2.loc[k_A3_2['sig_k'].astype(float) > 0]
k_A3_3 = k_A3_3.loc[k_A3_3['sig_k'].astype(float) > 0]

k_B2_1 = k_B2_1.loc[k_B2_1['sig_k'].astype(float) > 0]
k_B2_2 = k_B2_2.loc[k_B2_2['sig_k'].astype(float) > 0]
k_B2_3 = k_B2_3.loc[k_B2_3['sig_k'].astype(float) > 0]

k_M17_1 = k_M17_1.loc[k_M17_1['sig_k'].astype(float) > 0]
k_M17_2 = k_M17_2.loc[k_M17_2['sig_k'].astype(float) > 0]
k_M17_3 = k_M17_3.loc[k_M17_3['sig_k'].astype(float) > 0]

names = {'sig_L': 'sig_L_Set3',
         'sig_x0': 'sig_x0_Set3',
         'sig_k': 'sig_k_Set3',
         'sig_b': 'sig_b_Set3',
         'r2_sig': 'r2_sig_Set3',
         'rate_0': 'rate_0_Set3',
         'rate_1': 'rate_1_Set3',
         'rate_2': 'rate_2_Set3',
         'rate_4': 'rate_4_Set3'
         }

k_A3_a = k_A3_1.merge(k_A3_2, on='Name', how='outer', suffixes=('_Set1', '_Set2'))
k_A3 = k_A3_a.merge(k_A3_3, on='Name', how='outer').rename(columns = names)

k_A3['num_replicas'] = 3 - k_A3[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].isnull().sum(axis=1)
k_A3['k_avg'] = k_A3[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].mean(axis=1)
k_A3['k_std'] = k_A3[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].std(axis=1)
k_A3['k_var'] = k_A3['k_std'] / k_A3['k_avg']
k_A3['Set1_k_z'] = (k_A3['sig_k_Set1'] - k_A3['k_avg']) / k_A3['k_std']
k_A3['Set2_k_z'] = (k_A3['sig_k_Set2'] - k_A3['k_avg']) / k_A3['k_std']
k_A3['Set3_k_z'] = (k_A3['sig_k_Set3'] - k_A3['k_avg']) / k_A3['k_std']

k_A3['max_z'] = k_A3[['Set1_k_z','Set2_k_z','Set3_k_z']].abs().idxmax(axis = 1)

Set1_drops = k_A3.loc[((k_A3['k_var'] > var_thresh) & (k_A3['num_replicas'] == 3) & (k_A3['max_z'] == 'Set1_k_z'))]['Name'].tolist()
Set2_drops = k_A3.loc[((k_A3['k_var'] > var_thresh) & (k_A3['num_replicas'] == 3) & (k_A3['max_z'] == 'Set2_k_z'))]['Name'].tolist()
Set3_drops = k_A3.loc[((k_A3['k_var'] > var_thresh) & (k_A3['num_replicas'] == 3) & (k_A3['max_z'] == 'Set3_k_z'))]['Name'].tolist()

k_A3_1_f = k_A3_1[~k_A3_1['Name'].isin(Set1_drops)]
k_A3_2_f = k_A3_2[~k_A3_2['Name'].isin(Set2_drops)]
k_A3_3_f = k_A3_3[~k_A3_3['Name'].isin(Set3_drops)]

k_A3_a = k_A3_1_f.merge(k_A3_2_f, on='Name', how='outer', suffixes=('_Set1', '_Set2'))
k_A3_f = k_A3_a.merge(k_A3_3_f, on='Name', how='outer').rename(columns = names)

k_A3_f['num_replicas'] = 3 - k_A3_f[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].isnull().sum(axis=1)
k_A3_f['k_avg'] = k_A3_f[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].mean(axis=1)
k_A3_f['k_std'] = k_A3_f[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].std(axis=1)
k_A3_f['k_var'] = k_A3_f['k_std'] / k_A3_f['k_avg']
k_A3_f['avg_r2'] = k_A3_f[['r2_sig_Set1','r2_sig_Set2','r2_sig_Set3']].mean(axis=1)

k_A3_f = k_A3_f.loc[k_A3_f['k_var'] < var_thresh]
k_A3_final = k_A3_f[['Name','k_avg','k_std','k_var','avg_r2']]

k_B2_a = k_B2_1.merge(k_B2_2, on='Name', how='outer', suffixes=('_Set1', '_Set2'))
k_B2 = k_B2_a.merge(k_B2_3, on='Name', how='outer').rename(columns = names)

k_B2['num_replicas'] = 3 - k_B2[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].isnull().sum(axis=1)
k_B2['k_avg'] = k_B2[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].mean(axis=1)
k_B2['k_std'] = k_B2[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].std(axis=1)
k_B2['k_var'] = k_B2['k_std'] / k_B2['k_avg']
k_B2['Set1_k_z'] = (k_B2['sig_k_Set1'] - k_B2['k_avg']) / k_B2['k_std']
k_B2['Set2_k_z'] = (k_B2['sig_k_Set2'] - k_B2['k_avg']) / k_B2['k_std']
k_B2['Set3_k_z'] = (k_B2['sig_k_Set3'] - k_B2['k_avg']) / k_B2['k_std']

k_B2['max_z'] = k_B2[['Set1_k_z','Set2_k_z','Set3_k_z']].abs().idxmax(axis = 1)

Set1_drops = k_B2.loc[((k_B2['k_var'] > var_thresh) & (k_B2['num_replicas'] == 3) & (k_B2['max_z'] == 'Set1_k_z'))]['Name'].tolist()
Set2_drops = k_B2.loc[((k_B2['k_var'] > var_thresh) & (k_B2['num_replicas'] == 3) & (k_B2['max_z'] == 'Set2_k_z'))]['Name'].tolist()
Set3_drops = k_B2.loc[((k_B2['k_var'] > var_thresh) & (k_B2['num_replicas'] == 3) & (k_B2['max_z'] == 'Set3_k_z'))]['Name'].tolist()

k_B2_1_f = k_B2_1[~k_B2_1['Name'].isin(Set1_drops)]
k_B2_2_f = k_B2_2[~k_B2_2['Name'].isin(Set2_drops)]
k_B2_3_f = k_B2_3[~k_B2_3['Name'].isin(Set3_drops)]

k_B2_a = k_B2_1_f.merge(k_B2_2_f, on='Name', how='outer', suffixes=('_Set1', '_Set2'))
k_B2_f = k_B2_a.merge(k_B2_3_f, on='Name', how='outer').rename(columns = names)

k_B2_f['num_replicas'] = 3 - k_B2_f[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].isnull().sum(axis=1)
k_B2_f['k_avg'] = k_B2_f[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].mean(axis=1)
k_B2_f['k_std'] = k_B2_f[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].std(axis=1)
k_B2_f['k_var'] = k_B2_f['k_std'] / k_B2_f['k_avg']
k_B2_f['avg_r2'] = k_B2_f[['r2_sig_Set1','r2_sig_Set2','r2_sig_Set3']].mean(axis=1)

k_B2_f = k_B2_f.loc[k_B2_f['k_var'] < var_thresh]
k_B2_final = k_B2_f[['Name','k_avg','k_std','k_var','avg_r2']]


k_M17_a = k_M17_1.merge(k_M17_2, on='Name', how='outer', suffixes=('_Set1', '_Set2'))
k_M17 = k_M17_a.merge(k_M17_3, on='Name', how='outer').rename(columns = names)

k_M17['num_replicas'] = 3 - k_M17[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].isnull().sum(axis=1)
k_M17['k_avg'] = k_M17[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].mean(axis=1)
k_M17['k_std'] = k_M17[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].std(axis=1)
k_M17['k_var'] = k_M17['k_std'] / k_M17['k_avg']
k_M17['Set1_k_z'] = (k_M17['sig_k_Set1'] - k_M17['k_avg']) / k_M17['k_std']
k_M17['Set2_k_z'] = (k_M17['sig_k_Set2'] - k_M17['k_avg']) / k_M17['k_std']
k_M17['Set3_k_z'] = (k_M17['sig_k_Set3'] - k_M17['k_avg']) / k_M17['k_std']

k_M17['max_z'] = k_M17[['Set1_k_z','Set2_k_z','Set3_k_z']].abs().idxmax(axis = 1)

Set1_drops = k_M17.loc[((k_M17['k_var'] > var_thresh) & (k_M17['num_replicas'] == 3) & (k_M17['max_z'] == 'Set1_k_z'))]['Name'].tolist()
Set2_drops = k_M17.loc[((k_M17['k_var'] > var_thresh) & (k_M17['num_replicas'] == 3) & (k_M17['max_z'] == 'Set2_k_z'))]['Name'].tolist()
Set3_drops = k_M17.loc[((k_M17['k_var'] > var_thresh) & (k_M17['num_replicas'] == 3) & (k_M17['max_z'] == 'Set3_k_z'))]['Name'].tolist()

k_M17_1_f = k_M17_1[~k_M17_1['Name'].isin(Set1_drops)]
k_M17_2_f = k_M17_2[~k_M17_2['Name'].isin(Set2_drops)]
k_M17_3_f = k_M17_3[~k_M17_3['Name'].isin(Set3_drops)]

k_M17_a = k_M17_1_f.merge(k_M17_2_f, on='Name', how='outer', suffixes=('_Set1', '_Set2'))
k_M17_f = k_M17_a.merge(k_M17_3_f, on='Name', how='outer').rename(columns = names)

k_M17_f['num_replicas'] = 3 - k_M17_f[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].isnull().sum(axis=1)
k_M17_f['k_avg'] = k_M17_f[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].mean(axis=1)
k_M17_f['k_std'] = k_M17_f[['sig_k_Set1','sig_k_Set2','sig_k_Set3']].std(axis=1)
k_M17_f['k_var'] = k_M17_f['k_std'] / k_M17_f['k_avg']
k_M17_f['avg_r2'] = k_M17_f[['r2_sig_Set1','r2_sig_Set2','r2_sig_Set3']].mean(axis=1)

k_M17_f = k_M17_f.loc[k_M17_f['k_var'] < var_thresh]
k_M17_final = k_M17_f[['Name','k_avg','k_std','k_var','avg_r2']]

# merge into final output

names = {'k_avg':'k_avg_M17',
         'k_std':'k_std_M17',
         'k_var':'k_var_M17',
         'avg_r2':'avg_r2_M17'
         }

k_final_a = k_A3_final.merge(k_B2_final, on='Name', how='outer', suffixes=('_A3', 'B2'))
k_final = k_final_a.merge(k_M17_final, on='Name', how='outer').rename(columns = names)
k_final['ID'] = k_final['Name'].str.split('_').str[3]
k_final.to_csv(results_dir + 'targets_{0}_{1}.csv'.format(r2_thresh, var_thresh), index=False)

new_db = pd.read_csv(results_dir + 'mart_export.txt')
new_db['gene_stable_id'] = new_db['Gene stable ID version'].str.split('.').str[0]
new_db_gene = new_db[['gene_stable_id','Gene name']].drop_duplicates()
k_final_db = k_final.merge(new_db_gene, left_on='ID', right_on='gene_stable_id', how='left')
k_final_db.to_csv(results_dir+'targets_db_{0}_{1}.csv', index=False)

######################################################
# new analysis of grand slam data
######################################################

import pandas as pd
import numpy as np
import scipy.optimize
from matplotlib import pyplot as plt

data_dir = 'E:\\slamseq_run_9_9_2021\\count\\'
results_dir = 'E:\\slamseq_run_9_9_2021\\jons_analysis\\'

df = pd.read_csv('C:\\Users\\Jon\\Downloads\\test2\\24h.tsv', sep = '\t')

cells = ['A3','B2','M17']
times = [1,2,4]
runs = ['1','2','3']
z_thresh = 1.05
max_HL = 3

for cell in cells:
    for time in times:
        for run in runs:
            df['Set{0}_{1}_T{2}_R'.format(run, cell, time)] = df['/mnt/e/slamseq_run_9_9_2021/filter/Set{0}-{1}-{2} MAP'.format(run, cell, time)] - df['/mnt/e/slamseq_run_9_9_2021/filter/Set{0}-{1}-0 MAP'.format(run, cell)]
            df['Set{0}_{1}_T{2}_HL'.format(run, cell, time)] = -time*np.log(2)/np.log(1-df['Set{0}_{1}_T{2}_R'.format(run, cell, time)])
            
for cell in cells:
    for time in times:
        df['{0}_T{1}_POS'.format(cell, time)] = (df[['Set1_{0}_T{1}_HL'.format(cell, time), 'Set2_{0}_T{1}_HL'.format(cell, time),'Set3_{0}_T{1}_HL'.format(cell, time)]] >0).sum(1)
        df['{0}_T{1}_AVG'.format(cell, time)] = df[['Set1_{0}_T{1}_HL'.format(cell, time),'Set2_{0}_T{1}_HL'.format(cell, time),'Set3_{0}_T{1}_HL'.format(cell, time)]].mean(axis=1)
        df['{0}_T{1}_STD'.format(cell, time)] = df[['Set1_{0}_T{1}_HL'.format(cell, time),'Set2_{0}_T{1}_HL'.format(cell, time),'Set3_{0}_T{1}_HL'.format(cell, time)]].std(axis=1)
        df['{0}_T{1}_Set1_Z'.format(cell, time)] = abs((df['Set1_{0}_T{1}_HL'.format(cell, time)] - df['{0}_T{1}_AVG'.format(cell, time)]) / df['{0}_T{1}_STD'.format(cell, time)])
        df['{0}_T{1}_Set2_Z'.format(cell, time)] = abs((df['Set2_{0}_T{1}_HL'.format(cell, time)] - df['{0}_T{1}_AVG'.format(cell, time)]) / df['{0}_T{1}_STD'.format(cell, time)])
        df['{0}_T{1}_Set3_Z'.format(cell, time)] = abs((df['Set3_{0}_T{1}_HL'.format(cell, time)] - df['{0}_T{1}_AVG'.format(cell, time)]) / df['{0}_T{1}_STD'.format(cell, time)])
        
        df['Set1_{0}_T{1}_REALHL'.format(cell, time)] = df['Set1_{0}_T{1}_HL'.format(cell, time)]
        df['Set2_{0}_T{1}_REALHL'.format(cell, time)] = df['Set2_{0}_T{1}_HL'.format(cell, time)]
        df['Set3_{0}_T{1}_REALHL'.format(cell, time)] = df['Set3_{0}_T{1}_HL'.format(cell, time)]
        
        df[['Set1_{0}_T{1}_REALHL'.format(cell, time),'Set2_{0}_T{1}_REALHL'.format(cell, time),'Set3_{0}_T{1}_REALHL'.format(cell, time)]] = df[['Set1_{0}_T{1}_REALHL'.format(cell, time),'Set2_{0}_T{1}_REALHL'.format(cell, time),'Set3_{0}_T{1}_REALHL'.format(cell, time)]].replace([np.inf, -np.inf], np.nan)
        
        df.loc[df['{0}_T{1}_Set1_Z'.format(cell, time)] > z_thresh, 'Set1_{0}_T{1}_REALHL'.format(cell, time)] = np.nan
        df.loc[df['{0}_T{1}_Set2_Z'.format(cell, time)] > z_thresh, 'Set2_{0}_T{1}_REALHL'.format(cell, time)] = np.nan
        df.loc[df['{0}_T{1}_Set3_Z'.format(cell, time)] > z_thresh, 'Set3_{0}_T{1}_REALHL'.format(cell, time)] = np.nan
        
        df.loc[df['Set1_{0}_T1_HL'.format(cell, time)] > 5, 'Set1_{0}_T1_REALHL'.format(cell, time)] = np.nan
        df.loc[df['Set2_{0}_T1_HL'.format(cell, time)] > 5, 'Set2_{0}_T1_REALHL'.format(cell, time)] = np.nan
        df.loc[df['Set3_{0}_T1_HL'.format(cell, time)] > 5, 'Set3_{0}_T1_REALHL'.format(cell, time)] = np.nan
        
        df.loc[df['Set1_{0}_T2_HL'.format(cell, time)] > 10, 'Set1_{0}_T2_REALHL'.format(cell, time)] = np.nan
        df.loc[df['Set2_{0}_T2_HL'.format(cell, time)] > 10, 'Set2_{0}_T2_REALHL'.format(cell, time)] = np.nan
        df.loc[df['Set3_{0}_T2_HL'.format(cell, time)] > 10, 'Set3_{0}_T2_REALHL'.format(cell, time)] = np.nan
        
        df.loc[df['Set1_{0}_T4_HL'.format(cell, time)] > 15, 'Set1_{0}_T4_REALHL'.format(cell, time)] = np.nan
        df.loc[df['Set2_{0}_T4_HL'.format(cell, time)] > 15, 'Set2_{0}_T4_REALHL'.format(cell, time)] = np.nan
        df.loc[df['Set3_{0}_T4_HL'.format(cell, time)] > 15, 'Set3_{0}_T4_REALHL'.format(cell, time)] = np.nan
        
        df.loc[df['Set1_{0}_T{1}_HL'.format(cell, time)] <= 0, 'Set1_{0}_T{1}_REALHL'.format(cell, time)] = np.nan
        df.loc[df['Set2_{0}_T{1}_HL'.format(cell, time)] <= 0, 'Set2_{0}_T{1}_REALHL'.format(cell, time)] = np.nan
        df.loc[df['Set3_{0}_T{1}_HL'.format(cell, time)] <= 0, 'Set3_{0}_T{1}_REALHL'.format(cell, time)] = np.nan
        
        df['{0}_T{1}_HL_AVG'.format(cell, time)] = (df[['Set1_{0}_T{1}_REALHL'.format(cell, time), 'Set2_{0}_T{1}_REALHL'.format(cell, time),'Set3_{0}_T{1}_REALHL'.format(cell, time)]]).mean(axis=1)
        df['{0}_T{1}_HL_STD'.format(cell, time)] = (df[['Set1_{0}_T{1}_REALHL'.format(cell, time), 'Set2_{0}_T{1}_REALHL'.format(cell, time),'Set3_{0}_T{1}_REALHL'.format(cell, time)]]).std(axis=1)
        df['{0}_T{1}_HL_NUM'.format(cell, time)] = (df[['Set1_{0}_T{1}_REALHL'.format(cell, time), 'Set2_{0}_T{1}_REALHL'.format(cell, time),'Set3_{0}_T{1}_REALHL'.format(cell, time)]] >0).sum(1)

# make flags for each requirement and a final flag for "valid half lifes that pass all the requirements"

new_df = df[['Gene',
             'Symbol',
             'Length',
             'Set1_A3_T1_HL',
             'Set2_A3_T1_HL',
             'Set3_A3_T1_HL',
             'A3_T1_HL_NUM',
             'A3_T1_HL_AVG',
             'A3_T1_HL_STD',
             'Set1_A3_T2_HL',
             'Set2_A3_T2_HL',
             'Set3_A3_T2_HL',
             'A3_T2_HL_NUM',
             'A3_T2_HL_AVG',
             'A3_T2_HL_STD',
             'Set1_A3_T4_HL',
             'Set2_A3_T4_HL',
             'Set3_A3_T4_HL',
             'A3_T4_HL_NUM',
             'A3_T4_HL_AVG',
             'A3_T4_HL_STD',
             'Set1_B2_T1_HL',
             'Set2_B2_T1_HL',
             'Set3_B2_T1_HL',
             'B2_T1_HL_NUM',
             'B2_T1_HL_AVG',
             'B2_T1_HL_STD',
             'Set1_B2_T2_HL',
             'Set2_B2_T2_HL',
             'Set3_B2_T2_HL',
             'B2_T2_HL_NUM',
             'B2_T2_HL_AVG',
             'B2_T2_HL_STD',
             'Set1_B2_T4_HL',
             'Set2_B2_T4_HL',
             'Set3_B2_T4_HL',
             'B2_T4_HL_NUM',
             'B2_T4_HL_AVG',
             'B2_T4_HL_STD',
             'Set1_M17_T1_HL',
             'Set2_M17_T1_HL',
             'Set3_M17_T1_HL',
             'M17_T1_HL_NUM',
             'M17_T1_HL_AVG',
             'M17_T1_HL_STD',
             'Set1_M17_T2_HL',
             'Set2_M17_T2_HL',
             'Set3_M17_T2_HL',
             'M17_T2_HL_NUM',
             'M17_T2_HL_AVG',
             'M17_T2_HL_STD',
             'Set1_M17_T4_HL',
             'Set2_M17_T4_HL',
             'Set3_M17_T4_HL',
             'M17_T4_HL_NUM',
             'M17_T4_HL_AVG',
             'M17_T4_HL_STD'
             ]]

new_df.to_csv('C:\\Users\\Jon\\Downloads\\test2\\test.csv', index=False)


###################################################
# previous post-processing analysis
# group and average half lifes by group
###################################################

A3_HL_n = pd.read_csv(results_dir + 'half_lifes_A3.csv')
B2_HL_n = pd.read_csv(results_dir + 'half_lifes_B2.csv')
HCT_HL_n = pd.read_csv(results_dir + 'half_lifes_HCT.csv')

A3_HL_n['ID'] = A3_HL_n['Name'].str.split("_").str[0]
B2_HL_n['ID'] = B2_HL_n['Name'].str.split("_").str[0]
HCT_HL_n['ID'] = HCT_HL_n['Name'].str.split("_").str[0]

merged1 = A3_HL_n[['ID','half_life']].merge(B2_HL_n[['ID','half_life']], on='ID', suffixes=('_A3','_B2'), how='outer')
merged2 = merged1.merge(HCT_HL_n[['ID','half_life']], on='ID', how='outer')
merged2.columns = ['ID','half_life_A3','half_life_B2','half_life_HCT']
merged3 = merged2.merge(new_db, left_on='ID', right_on='Transcript stable ID version', how='left')
merged3.to_csv(results_dir + 'half_lifes_merged.csv', index=False)

new_db.loc[new_db['Transcript stable ID version'] == 'ENST00000666741.1']
HCT_HL_n.loc[HCT_HL_n['ID'] == 'ENST00000666741.1']

merged3_grouped = merged3.groupby(['ID']).agg({'half_life_A3':[np.mean,np.std],\
                                                                         'half_life_B2':[np.mean,np.std],\
                                                                         'half_life_HCT':[np.mean,np.std]})
merged3_grouped.reset_index(inplace=True)
merged3_grouped.to_csv(results_dir + 'half_lifes_averaged_by_ID.csv', index=False)

# get transcripts in each gene/transcript type group

transcripts_in_groups = pd.DataFrame(merged3.groupby(['Gene name', "Transcript type"])['Transcript name'].apply(lambda x: "{%s}" % ', '.join(x)))
transcripts_in_groups.reset_index(inplace=True)
transcripts_in_groups.columns = ['Gene name','Transcript type','Transcripts in group']
merged4 = merged3_grouped.merge(transcripts_in_groups, on=['Gene name','Transcript type'], how='left')
merged4.to_csv(data_dir + 'grouped_by_gene_name_and_transcript_type.csv', index=False)

# calculate variant lengths

len_dic = {}
transcripts = list(set(A3['Transcript name'].dropna().tolist()))
for t in transcripts:
    print(t)
    if t != "":
        df = A3.loc[A3['Transcript name'] == t]
        tot = []
        for index, row in df.iterrows():
            leng = [*range(int(row['Start']),int(row['End']))]
            tot.append(leng)
        tot1 = len(set(tot[0]))
        len_dic[t] = tot1

lengths_df = pd.DataFrame.from_dict(len_dic, orient='index')
lengths_df.reset_index(inplace=True)
lengths_df.columns = ['Transcript name','true_length']
merged3 = pd.read_csv(data_dir + 'half_life_differences.csv')
merged4 = merged3.merge(lengths_df, on='Transcript name', how='left')
merged4.to_csv(data_dir + 'half_lifes_with_true_length.csv', index=False)
