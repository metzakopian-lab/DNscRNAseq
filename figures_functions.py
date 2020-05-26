import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc

def readRtable(filename):
    scores = pd.read_csv(filename)
    scores = scores.set_index('Unnamed: 0')
    scores.index.name = 'index'
    return scores

def export_to_mtx(adata, key):
    from scipy.io import mmwrite
    mmwrite('paper-data/modelling/' + key + '.mtx', adata.raw.X)
    adata.obs.to_csv('paper-data/modelling/' + key + '-cell_metadata.csv')
    adata.var.to_csv('paper-data/modelling/' + key + '-gene_metadata.csv')

def split_ann(ada):
    return {
        'wt-all': ada[ada.obs['genotype'].isin(['WT'])].copy(),
        'wt-unt': ada[ada.obs['sample'] == 'WT_UNT'].copy(),
        'unt-rot': ada[ada.obs['sample'].isin(['WT_UNT','WT_ROT'])].copy(),
        'unt-rot-dan1-dan2':ada[ada.obs['sample'].isin([ 'WT_ROT', 'WT_UNT']) & ada.obs['Cell Type Classification'].isin(['DAn1','DAn2'])].copy(),
        'unt-rot-dan1': ada[ada.obs['sample'].isin(['WT_UNT','WT_ROT']) & (ada.obs['Cell Type Classification'] == 'DAn1')].copy(),
        'unt-rot-dan2':ada[ada.obs['sample'].isin([ 'WT_ROT', 'WT_UNT']) & ada.obs['Cell Type Classification'].isin(['DAn2'])].copy(),
        'unt-tun' : ada[ada.obs['treatment'].isin(['TUN','UNT']) & ada.obs['genotype'].isin(['WT'])].copy(),
        'unt-tun-dan1' : ada[ada.obs['treatment'].isin(['TUN','UNT']) & ada.obs['genotype'].isin(['WT']) & ada.obs['Cell Type Classification'].isin(['DAn1'])].copy(),
        'unt-tun-dan1-dan2' : ada[ada.obs['treatment'].isin(['TUN','UNT']) & ada.obs['genotype'].isin(['WT']) & ada.obs['Cell Type Classification'].isin(['DAn1','DAn2'])].copy(),
        'wt-het-unt-dan1-dan2': ada[ada.obs['genotype'].isin(['WT','A53T']) &  
                               ada.obs['Cell Type Classification'].isin(['DAn1','DAn2']) &
                                   ada.obs['treatment'].isin(['UNT'])].copy(),
        'unt-tun-dan2' : ada[ada.obs['treatment'].isin(['TUN','UNT']) & ada.obs['genotype'].isin(['WT']) & ada.obs['Cell Type Classification'].isin(['DAn2'])].copy(),
        'unt-het': ada[ada.obs['genotype'].isin(['WT','A53T'])& ada.obs['treatment'].isin(['UNT'])].copy(),
        'het-tun': ada[ada.obs['genotype'].isin(['WT','A53T'])& ada.obs['treatment'].isin(['UNT','TUN'])].copy(),
        'het-rot': ada[ada.obs['genotype'].isin(['WT','A53T'])& ada.obs['treatment'].isin(['UNT','ROT'])].copy(),

        'unt-het-dan1': ada[ada.obs['genotype'].isin(['WT','A53T'])& 
                                ada.obs['treatment'].isin(['UNT']) & 
                                ada.obs['Cell Type Classification'].isin(['DAn1'])].copy(),
        'rot-het-dan1': ada[ada.obs['genotype'].isin(['WT','A53T'])& 
                                ada.obs['treatment'].isin(['ROT']) & 
                                ada.obs['Cell Type Classification'].isin(['DAn1'])].copy(),
        'rot-het-dan2': ada[ada.obs['genotype'].isin(['WT','A53T'])& 
                                ada.obs['treatment'].isin(['ROT']) & 
                                ada.obs['Cell Type Classification'].isin(['DAn2'])].copy(),
        'tun-het-dan1': ada[ada.obs['genotype'].isin(['WT','A53T'])& 
                                ada.obs['treatment'].isin(['TUN']) & 
                                ada.obs['Cell Type Classification'].isin(['DAn1'])].copy(),
        'tun-het-dan2': ada[ada.obs['genotype'].isin(['WT','A53T'])& 
                                ada.obs['treatment'].isin(['TUN']) & 
                                ada.obs['Cell Type Classification'].isin(['DAn2'])].copy(),
        'wt-het-dan1': ada[ada.obs['genotype'].isin(['WT','A53T']) &  
                               ada.obs['Cell Type Classification'].isin(['DAn1'])].copy(),
        'wt-het-dan1-rot': ada[ada.obs['genotype'].isin(['WT','A53T']) &  
                               ada.obs['Cell Type Classification'].isin(['DAn1']) &
                               ada.obs['treatment'].isin(['UNT','ROT'])].copy(),
        'wt-het-dan1-dan2': ada[ada.obs['genotype'].isin(['WT','A53T']) &  
                               ada.obs['Cell Type Classification'].isin(['DAn1','DAn2'])].copy(),
        'wt-het-tun-dan1-dan2': ada[ada.obs['genotype'].isin(['WT','A53T']) &
                               ada.obs['Cell Type Classification'].isin(['DAn1','DAn2']) & ada.obs['treatment'].isin(['UNT','TUN'])].copy(),
        'all-samples': ada}



def plotheatmap(ann_heatmap, genes, groupby, filename,  normalise_per_cell = True, vmin = -1.5, vmax = 1.5,figsize=(23,2)):    
    ann_heatmap = ann_heatmap.copy()
    ann_heatmap.X = ann_heatmap.raw.X
    if normalise_per_cell:
        sc.pp.normalize_per_cell(ann_heatmap, counts_per_cell_after = 15000, copy = False)
        
    ann_heatmap = ann_heatmap[:,ann_heatmap.var_names.isin(genes)].copy()

    records = {}
    htmap_df = pd.DataFrame(columns=ann_heatmap.var_names)
    #nn_heatmap.
    for grp in ann_heatmap.obs[groupby].unique():
        ctrl_group = ann_heatmap[ann_heatmap.obs[groupby] != grp].copy()
        test_group = ann_heatmap[ann_heatmap.obs[groupby] == grp].copy()
        records[grp] = (np.log2(test_group.X.mean(axis = 0) + 1)- np.log2(ctrl_group.X.mean(axis = 0) + 1)).tolist()[0]
    htmap_df = pd.DataFrame.from_dict(records, orient = 'index', columns=list(ann_heatmap.var_names))
    htmap_df = htmap_df.sort_index()
    htmap_df = htmap_df[genes]
    f,ax = plt.subplots(1,1, figsize = figsize)
    sns.heatmap(htmap_df, 
                square= False, 
                cmap = "coolwarm", 
                vmin = vmin, 
                vmax = vmax, 
                linewidths=.5, 
                cbar_kws=dict(
                              use_gridspec = False,
                              aspect = 8,
                              anchor = (-0.3, 0.0)),
                ax = ax)
    fn = filename
    if normalise_per_cell:
        fn += '-seqdepth-norm'
    f.savefig(fn+'.png', bbox_inches='tight', dpi = 300)
    plt.clf()

    
def plotheatmap_hypertest(ann_heatmap, genes, groupby, filename, vmin = 0, vmax = 1.5,figsize=(23,2)):
    from scipy.stats import hypergeom
    ann_heatmap = ann_heatmap.copy()
    ann_heatmap.X = ann_heatmap.raw.X
    ann_heatmap = ann_heatmap[:, ann_heatmap.var_names.isin(genes)].copy()
    ann_heatmap.var['ncells'] = (ann_heatmap.X > 0).sum(axis = 0).tolist()[0]
    records = {}
    htmap_df = pd.DataFrame(columns=ann_heatmap.var_names)
    #nn_heatmap.
    M = len(ann_heatmap.obs)
    for grp in ann_heatmap.obs[groupby].unique():        
        test_group = ann_heatmap[ann_heatmap.obs[groupby] == grp].copy()
        test_group.var['ncells'] = (test_group.X > 0).sum(axis = 0).tolist()[0]
        
        N = len(test_group.obs)
        records[grp] = []
        for g in genes:
            n = ann_heatmap.var.loc[g]['ncells']
            x = test_group.var.loc[g]['ncells']
            # add a small value to avoid 0 as a pvalue
            hyper_test_pval = 1 - hypergeom.cdf(x,M,n,N) + 10e-30 
            records[grp].append(hyper_test_pval)
    
    htmap_df = pd.DataFrame.from_dict(records, orient = 'index', columns = genes)
    htmap_df = htmap_df.sort_index()
    f,ax = plt.subplots(1,1, figsize = figsize)
    sns.heatmap(htmap_df.applymap(lambda x : -np.log10(x)), 
                square= False, 
                cmap = sns.light_palette('red', as_cmap=True), 
                vmin = vmin, 
                vmax = vmax, 
                linewidths=.5, 
                cbar_kws=dict(
                              use_gridspec = False,
                              aspect = 8,label = '$-log(p_{val})$',
                              anchor = (-0.3, 0.0)),
                ax = ax, )
    fn = filename
    f.savefig(fn+'.pdf', bbox_inches='tight', dpi = 300)
    plt.close(f)
    return htmap_df


def violin_from_dict(ann_violin, 
                     dict_list, 
                     category_label, 
                     prefix,
                     taskid,
                     colormap = None,
                     figsize = (1.2,1.2),
                    lfc = False):
    from scipy.stats import ranksums
    ann_violin = ann_violin.copy()
    ann_violin.X = ann_violin.raw.X
    sc.pp.normalize_per_cell(ann_violin, copy = False, counts_per_cell_after = 15000)
    
    for t, _genes in dict_list.items():
        if len(_genes) > 1:
            fig_all, ax_sub = plt.subplots(1,len(_genes), figsize = ((figsize[0]*len(_genes)), figsize[1] + 0.3))
        for _, g in enumerate(_genes):
            
            if not g in ann_violin.var_names:
                print(g,'not found')
                continue
            ann_violin.obs['exp'] = ann_violin.X[:,ann_violin.var_names == g].A.reshape(-1)

            ann_violin.obs['l2fc'] = np.log2(ann_violin.X[:,ann_violin.var_names == g].A.reshape(-1) + 1)
            ann_violin.obs[g] = ann_violin.obs['l2fc']
            # Figure properties
            fig, ax = plt.subplots(figsize=figsize)
                #                rc={'font.size': 32, 'axes.labelsize': 18, 'legend.fontsize': 18, 
                #    'axes.titlesize': 20, 'xtick.labelsize': 20, 'ytick.labelsize': 20}
                #plt.rcParams.update(**rc)
                # Violin Plot
            
            mask = ann_violin.obs[category_label] == ann_violin.obs[category_label].cat.categories[0]
            if lfc:
                rep = np.log2(ann_violin.obs['exp'][~mask].mean() + 1) - np.log2(ann_violin.obs['exp'][mask].mean() + 1)
            else:
                stat,pval = ranksums(ann_violin.obs['l2fc'][mask], ann_violin.obs['l2fc'][~mask])
                rep = pval
                    
            from statannot import add_stat_annotation
            if ann_violin.obs[category_label].dtype.name =='category':
                ann_violin.obs[category_label] =  ann_violin.obs[category_label].cat.remove_unused_categories()
            axs = [ax]
            if len(_genes) > 1:
                axs.append(ax_sub[_])
            for _ax in axs:
                sns.violinplot(
                    data = ann_violin.obs,
                    palette = colormap, 
                    y = g, 
                    x = category_label,
                    linewidth=1,
                    ax = _ax)
                if lfc:
                    add_stat_annotation(_ax,
                                        data = ann_violin.obs, 
                                        y = g, 
                                        x = category_label,
                                        box_pairs=[(ann_violin.obs[category_label].unique())],
                                        perform_stat_test=False, 
                                        pvalues=[rep],
                                        text_format = 'custom',
                                        line_offset_to_box=0.2, line_offset=0.1, line_height=0.05, linewidth=0.6, text_offset = 0.5
                    )
                else:
                    add_stat_annotation(_ax,
                                        data = ann_violin.obs, 
                                        y = g, 
                                        x = category_label,
                                        box_pairs=[(ann_violin.obs[category_label].unique())],
                                        perform_stat_test=False, pvalues=[rep],
                                        line_offset_to_box=0.2, line_offset=0.1, line_height=0.05, linewidth=0.6, text_offset = 0.5
                    )

                    
                
                    
                    
                _ax.set_xlabel('')
                _ax.set_title(g)
                _ax.set_ylabel('')
            
            ax.set_ylabel(r'$ \log_{2}( expression) $')
            if _ == 0 and len(_genes) > 1:
                ax_sub[_].set_ylabel(r'$ \log_{2}( expression) $', fontsize = 7 )
            if len(_genes) > 1:
                ax_sub[_].tick_params(axis = 'y', pad = -3)
            path = prefix + t + '-' + g +'-'+taskid+'-'+ ".pdf"
            fig.savefig(path, dpi = 300, bbox_inches='tight')
            plt.close(fig)
        
        fig_all.tight_layout(pad = 0.3)
        fig_all.savefig(prefix + t +'-'+taskid+'-'+ ".pdf", bbox_inches = 'tight')
        plt.close('all')
            

#violin('unt-rot-dan1', 'test', 'treatment', colormap = 'ROT treatment',  dirname='test/')
                
def rank_genes(ann_test, gene_names, level_names, compare, reference, groups = 'all', normalize = True, stat = False):
    import itertools
    import scanpy as sc
    ann_test = ann_test.copy()
    ann_test.X = ann_test.raw.X
    if normalize:
        sc.pp.normalize_per_cell(ann_test, counts_per_cell_after= 15000, copy=False)
        
        
    ann_test = ann_test[:,ann_test.var_names.isin(gene_names)].copy()
    
    levels = {}
    for l in level_names:
        levels[l] = list(ann_test.obs[l].unique())
    
    level_keys = list(levels.keys())
    for cvar in itertools.product(*list(levels.values())):
        ann_tmp = ann_test
        for key, value in zip(level_keys, cvar):
            ann_tmp = ann_tmp[ann_tmp.obs[key] == value]
        group_cardinality = len(ann_tmp.obs[compare].unique())

        for grp in ann_tmp.obs[compare].unique():
            if grp == reference:
                continue
            if groups != 'all' and not (grp in groups):
                continue

            if group_cardinality == 2:
                label_str = '_'.join(cvar)
            else:
                label_str = '_'.join(cvar) + '_' + grp
            if stat:
                sc.tl.rank_genes_groups(ann_tmp, compare, n_genes = ann_tmp.X.shape[1], reference = reference, groups = [grp], use_raw = False)
                gene_groups = ann_tmp.uns['rank_genes_groups']
                group_name = gene_groups['names'].dtype.names[0]
                res_df = pd.DataFrame({fname: gene_groups[fname][group_name] for fname in ['pvals']}, index = gene_groups['names'][group_name])
                yield(label_str, res_df)
            else:
                yield(label_str, ann_tmp[ann_tmp.obs[compare] == reference].copy(), ann_tmp[ann_tmp.obs[compare] == grp].copy())

        
        

def plotheatmap_pairs(gct, genes, filename, title ,vmin = -1.5, vmax = 1.5, figsize = (23,2)):    
    records = {}
    gl = None
    gct_size = 0
    gene_size = max(map(len, genes))
    groups_size = 0
    for grp, ctrl_group, test_group in gct:
        gct_size += 1
        if gl is None:
            gl = list(ctrl_group.var_names)
            groups_size = 1 + grp.count("_")
        records[grp] = (np.log2(test_group.X.mean(axis = 0) + 1)- np.log2(ctrl_group.X.mean(axis = 0) + 1)).reshape(-1).tolist()[0]
    htmap_df = pd.DataFrame.from_dict(records, orient = 'index', columns = gl)
    htmap_df = htmap_df.sort_index(axis = 0)
    # sort gene values
    htmap_df = htmap_df[genes]
    

    #est_figsize = (3+len(genes)*0.3 + 0.5*groups_size - (gene_size - gct_size)*0.05,  (gct_size*0.3) + (gene_size - gct_size)*0.1 + (0.05*gene_size))
    #print(est_figsize)
    f,ax = plt.subplots(1,1, figsize = figsize)
    sns.heatmap(htmap_df, 
                square= True, 
                cmap = "coolwarm", 
                vmin = vmin, 
                vmax = vmax, 
                linewidths=.5, 
                cbar_kws=dict(label = '$log_{2}(foldchange)$',
                              use_gridspec = False,
                              aspect = 8,
                              anchor = (-0.3, 0.0)),
                ax = ax, ).set_title(title)
    filename += '.pdf'
    print('Saving to',filename)
    f.savefig(filename, bbox_inches='tight', dpi = 300)
    plt.clf()
    return htmap_df


def plotheatmap_pairs_stat(gct, genes, filename, title ,vmin = 0, vmax = 30, figsize = (23,2)):    
    records = pd.DataFrame( index = genes)
    for grp, res_df in gct:
        records[grp] = res_df['pvals']
    htmap_df = records
    # sort gene values
    htmap_df = htmap_df.T
    htmap_df = htmap_df.sort_index()
    htmap_df = htmap_df[genes]
    
    #est_figsize = (3+len(genes)*0.3 + 0.5*groups_size - (gene_size - gct_size)*0.05,  (gct_size*0.3) + (gene_size - gct_size)*0.1 + (0.05*gene_size))
    #print(est_figsize)
    f,ax = plt.subplots(1,1, figsize = figsize)
    sns.heatmap(htmap_df.applymap(lambda x : -np.log10(x)), 
                square= True, 
                cmap = sns.light_palette('red', as_cmap=True), 
                vmin = vmin, 
                vmax = vmax, 
                linewidths=.5, 
                cbar_kws=dict(label = '$-log(p_{val})$',
                              use_gridspec = False,
                              aspect = 8,
                              anchor = (-0.3, 0.0)),
                ax = ax, ).set_title(title)
    filename += '.pdf'
    print('Saving to',filename)
    f.savefig(filename, bbox_inches='tight', dpi = 300)
    plt.clf()
    return htmap_df



