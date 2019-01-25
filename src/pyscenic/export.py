# coding=utf-8

import os
import numpy as np
import pandas as pd
import loompy as lp
#from sklearn.manifold.t_sne import TSNE
import umap
from typing import List, Mapping, Sequence, Optional
from operator import attrgetter
from multiprocessing import cpu_count
from itertools import chain, repeat, islice
import networkx as nx
import zlib
import base64
import json
import warnings
from .aucell import aucell
from .genesig import Regulon
from .binarization import binarize

def compress_encode(value):
    return base64.b64encode(zlib.compress(json.dumps(value).encode('ascii'))).decode('ascii')

def decompress_meta(meta):
    try:
        meta = meta.decode('ascii')
        return json.loads(zlib.decompress(base64.b64decode(meta)))
    except AttributeError:
        return json.loads(zlib.decompress(base64.b64decode(meta.encode('ascii'))).decode('ascii'))

def export2loom(ex_mtx: pd.DataFrame,
                regulons: List[Regulon],
                out_fname: str,
                cell_clusters: Optional[Mapping[str,str]]=None,
                cell_annotations: Optional[Mapping[str,str]]=None,
                tree_structure: Sequence[str] = (),
                title: Optional[str] = None,
                nomenclature: str = "Unknown",
                embeddings = None,
                auc_mtx = None,
                auc_thresholds = None,
                num_workers = cpu_count(),
                compress = True,
                replaceSpacesInRegulons=True):
    """
    Create a loom file for a single cell experiment to be used in SCope.

    :param ex_mtx: The expression matrix (n_cells x n_genes).
    :param regulons: A list of Regulons.
    :param out_fname: The name of the loom file to create.
    :param cell_clusters: A dataframe with cell cluster/type annotation (in SCOPE it can be searched in the 'gene' tab). Rows: cells; columns: each clustering or annotation (e.g. different resolutions...)).
    :param cell_annotations: A dataframe with cell annotations (in SCOPE it can be used in the 'compare' tab; they should be categorical). Rows: cells; columns: each cannotation (e.g. phenotype/condition, run date, ...)).
    :param tree_structure: A sequence of strings that defines the category tree structure.
    :param title: The title for this loom file. If None, the basename of the filename is used as the title.
    :param nomenclature: The name of the genome.
    :param embeddings: A dataframe with different embeddings: cells (rows) x coordinates for each embedding (columns). Each embedding (e.g. tSNE, UMAP, PCA, difussion map, ...) should consist of two columns with sufix '_X' and '_Y'. The prefix will be used as name for the embedding (e.g.: tsne_X, tsne_Y, umap_X, umap_Y). If 'None' an UMAP will be calculated with the default settings.
    :param auc_mtx: AUCell regulon activity matrix. If 'None' it will be calculated with the default settings.
    :param auc_thresholds: Binarized AUC thresholds for each regulon. If 'None' it will be calculated with the default settings.
    :param num_workers: The number of cores to use for AUCell regulon enrichment.
    :param compress: compress metadata (required for SCope).
    :param replaceSpacesInRegulons: Whether to rename regulons to remove spaces (required for current version of SCOPE).
    """
    # Information on the general loom file format: http://linnarssonlab.org/loompy/format/index.html
    # Information on the SCope specific alterations: https://github.com/aertslab/SCope/wiki/Data-Format

    ################################################
    ##### Check matching cell IDs
    cellIds_expMx = ex_mtx.index.values

    if not cell_clusters is None:
        if( not all(cellIds_expMx == cell_clusters.index.values) ):
             raise Exception('Cell IDs in the expression and the cell_clusters dataframe do not match.')

    if not cell_annotations is None:
        if( not all(cellIds_expMx == cell_annotations.index.values) ):
             raise Exception('Cell IDs in the expression and the cell_annotations dataframe do not match.')

    if not auc_mtx is None:
        if( not all(cellIds_expMx == auc_mtx.index.values) ):
             raise Exception('Cell IDs in the expression and the AUC matrices do not match.')

    if not embeddings is None:
        if( not all(cellIds_expMx == embeddings.index.values) ):
             raise Exception('Cell IDs in the expression matrix and the embeddings do not match.')



    ################################################
    #### Clusters / Cell type annotation
    if cell_clusters is None:
        cell_clusters = pd.DataFrame(data=['-']*ex_mtx.shape[0],
                               index=ex_mtx.index,
                               columns=['All cells'])

    # Encode cell type clusters
    clusterings=cell_clusters
    clusteringIds_dict=[]
    for i in range(0,clusterings.shape[1]):
        colName=clusterings.columns.values[i]
        name2idx = dict(map(reversed, enumerate(sorted(set(clusterings.to_dict()[colName].values())))))
        clusterings[colName] = clusterings[colName].replace(clusterings).replace(name2idx)
        clusteringIds_dict=clusteringIds_dict+[{
                     "id": i,
                     "group": colName,
                     "name": colName,
                     "clusters": [{"id": idx, "description": name} for name, idx in name2idx.items()]
                 }]
    clusterings.columns=[str(i) for i in range(0, clusterings.shape[1], 1)]

    # Encode annotations
    annots_dict=[{"name": "", "values": []}]

    if not cell_annotations is None:
        annots_dict=[]
        for i in range(0, cell_annotations.shape[1]):
            colName=cell_annotations.columns.values[i]
            annots_dict=annots_dict+[{
                         "name": colName,
                         "values": list(set(cell_annotations.to_dict()[colName].values()))
                     }]

    ################################################
    #### Regulons and AUC (calculate AUC if needed)
    # Check if regulon name is compatible with SCOPE
    regulonNames_regulons = [r.name for r in regulons]
    regulonNames_AucMx = [r for r in auc_mtx.columns.values]
    if( not set(regulonNames_regulons) == set(regulonNames_AucMx) ):
        raise Exception('Regulon names in the regulons and the AUC matrix do not match.')

    if(any([' ' in r for r in regulonNames_regulons]) & replaceSpacesInRegulons):
        regulons = [r.rename(r.name.replace(' ','_')) for r in regulons]
        auc_mtx=auc_mtx.rename(dict(zip(auc_mtx.columns.values, [r.replace(' ','_') for r in  auc_mtx.columns.values])), axis='columns')
        warnings.warn('Regulon names contain spaces, they have been replaced by "_".')
    del(regulonNames_regulons)
    del(regulonNames_AucMx)

    # Temporary until it is solved in SCOPE:
    if( not all(['_(' in r.name for r in regulons]) ):
        warnings.warn('Regulon names do not seem to be compatible with SCOPE. They should include "_(" to allow selection of the TF motif.')


    # Calculate regulon enrichment per cell using AUCell.
    if auc_mtx is None:
        auc_mtx = aucell(ex_mtx, regulons, num_workers=num_workers) # (n_cells x n_regulons)
        auc_mtx = auc_mtx.loc[ex_mtx.index]

    # Binarize matrix for AUC thresholds.
    if auc_thresholds is None:
        _, auc_thresholds = binarize(auc_mtx)

    # Calculate the number of genes per cell.
    binary_mtx = ex_mtx.copy()
    binary_mtx[binary_mtx != 0] = 1.0
    ngenes = binary_mtx.sum(axis=1).astype(int)

    # Encode genes in regulons as "binary" membership matrix.
    genes = np.array(ex_mtx.columns)
    n_genes = len(genes)
    n_regulons = len(regulons)
    data = np.zeros(shape=(n_genes, n_regulons), dtype=int)
    for idx, regulon in enumerate(regulons):
        data[:, idx] = np.isin(genes, regulon.genes).astype(int)
    regulon_assignment = pd.DataFrame(data=data,
                                      index=ex_mtx.columns,
                                      columns=list(map(attrgetter('name'), regulons)))

    # Encode motif logo
    def fetch_logo(context):
        for elem in context:
            if elem.endswith('.png'):
                return elem
        return ""

    name2logo = {reg.name: fetch_logo(reg.context) for reg in regulons}
    regulon_thresholds = [{"motifData": name2logo.get(name, ""),
                            "regulon": name,
                            "defaultThresholdValue":(threshold if isinstance(threshold, float) else threshold[0]),
                            "defaultThresholdName": "guassian_mixture_split",
                            "allThresholds": {"guassian_mixture_split": (threshold if isinstance(threshold, float) else threshold[0])}
                            } for name, threshold in auc_thresholds.iteritems()]

    ################################################
    ######## EMBEDDINGS
    if embeddings is None:
        # embeddings = pd.DataFrame(data=TSNE().fit_transform(auc_mtx),
        #                               index=ex_mtx.index, columns=['tSNE_X', 'tSNE_Y'])
        runUmap = umap.UMAP(n_neighbors=5, min_dist=0.5, metric='correlation').fit_transform
        embeddings = pd.DataFrame(data=runUmap(auc_mtx), index=auc_mtx.index, columns=['UMAP_X', 'UMAP_Y'])

    else:
        if(len(embeddings.columns)<2):
             raise Exception('The embedding should have at least two columns.')

    # Get default embedding (columns 0 and 1):
    if( (not '_X' in embeddings.columns[0]) | (not '_Y' in embeddings.columns[1]) ):
        raise Exception('The embedding columns should have the sufix _X and _Y.')
    default_embedding=embeddings.iloc[:,0:2]
    default_embedding_name=set([default_embedding.columns[0].replace('_X', ''), default_embedding.columns[0].replace('_X', '')])
    if(len(default_embedding_name)>1):
        raise Exception('The embedding columns should have the same prefix (e.g. tSNE_X , tSNE_Y.')
    default_embedding_name=default_embedding_name.pop()
    default_embedding.columns=['_X', '_Y']
    embeddings_names_dict=[{"id": -1, "name": default_embedding_name}]

    # Store remaining embeddings:
    embeddings=embeddings.iloc[:,2:]
    embeddingsX=None
    if(embeddings.shape[1]>0):
        embeddingsX = embeddings.iloc[:,embeddings.columns.str.contains('_X')]
        embeddingsY = embeddings.iloc[:,embeddings.columns.str.contains('_Y')]
        embeddings_names = embeddingsX.columns.str.replace('_X','')
        if(not (embeddings_names ==embeddingsY.columns.str.replace('_Y','')).all()):
            raise Exception('Embeddings prefixes for _X and _Y do not match')
        embeddings_names=embeddings_names.tolist()
        embeddingsX.columns = [str(i) for i in range(0, embeddingsX.shape[1], 1)]
        embeddingsY.columns = embeddingsX.columns
        embeddings_names_all=embeddings_names.insert(0, default_embedding_name)
        embeddings_names_dict=pd.DataFrame({'id':[i for i in range(-1, embeddingsX.shape[1], 1)], 'name':embeddings_names}).to_dict('records')

    ################################################
    ######## Save into loom file
    # Convert to numpy structured array
    def create_structure_array(df):
        return np.array([tuple(row) for row in df.values],
                        dtype=np.dtype(list(zip(df.columns, df.dtypes))))

    column_attrs = {
        "CellID": ex_mtx.index.values.astype('str'),
        "nGene": ngenes.values,
        "Embedding": create_structure_array(default_embedding),
        "RegulonsAUC": create_structure_array(auc_mtx),
        "Clusterings": create_structure_array(clusterings),
        "ClusterID": clusterings.values
        }

    if not embeddingsX is None:
        column_attrs['Embeddings_X'] = create_structure_array(embeddingsX)
        column_attrs['Embeddings_Y'] = create_structure_array(embeddingsY)

    if not cell_annotations is None:
        for i in range(0,cell_annotations.shape[1]):
            column_attrs[cell_annotations.columns[i]] = column_attrs[cell_annotations.columns[i]] = np.array(cell_annotations.iloc[:,i])

    row_attrs = {
        "Gene": ex_mtx.columns.values.astype('str'),
        "Regulons": create_structure_array(regulon_assignment),
        }

    general_attrs = {
            "title": os.path.splitext(os.path.basename(out_fname))[0] if title is None else title,
            "MetaData": {
                    "embeddings": embeddings_names_dict,
                    "annotations": annots_dict,
                    "clusterings": clusteringIds_dict,
                    "regulonThresholds": regulon_thresholds
                },
            "Genome": nomenclature}

    # Add tree structure.
    # All three levels need to be supplied
    assert len(tree_structure) <= 3, ""
    general_attrs.update(("SCopeTreeL{}".format(idx+1), category)
                         for idx, category in enumerate(list(islice(chain(tree_structure, repeat("")), 3))))

    # Compress MetaData global attribute (required for SCOPE)
    if compress:
        general_attrs["MetaData"] = compress_encode(value=general_attrs["MetaData"])

    # Create loom file for use with the SCope tool.
    # The loom file format opted for rows as genes to facilitate growth along the column axis (i.e add more cells)
    # PySCENIC chose a different orientation because of limitation set by the feather format: selectively reading
    # information from disk can only be achieved via column selection. For the ranking databases this is of utmost
    # importance.
    lp.create(filename=out_fname,
              layers=ex_mtx.T.values,
              row_attrs=row_attrs,
              col_attrs=column_attrs,
              file_attrs=general_attrs)

def export_asGmt(geneSet, fileName):
    """

    Export the regulons or co-expression modules. Format compatible to import from SCENIC in R.

    :param geneSet: Co-expression modules or regulons.
    :param fileName: The name of the file to create.
    """

    fileName=os.path.splitext(fileName)[0]+'.gmt'
    with open(fileName, 'w') as f:
        for reg in geneSet:
            f.write(reg.name + "\t.\t" + '\t'.join(list(reg.gene2weight.keys())) + '\n')

def export_motifs2txt(motifEnr, fileName):
    """

    Export the results from motif enrichment as text. Format compatible to import from SCENIC in R.

    :param motifEnr: Results from the motif enrichment
    :param fileName: The name of the file to create.
    """
    met=motifEnr['Enrichment']
    met.Context = [list(dbn)[2] for dbn in met.Context]
    met.TargetGenes=["; ".join(sorted([gs[0] for gs in row])) for row in met.TargetGenes]
    met=met.drop(columns='Annotation')
    met.columns.values[met.columns.values=='Context']='DB'
    met.to_csv(fileName, sep='\t')

def export_regulons(regulons: Sequence[Regulon], fname: str) -> None:
    """

    Export regulons as GraphML.

    :param regulons: The sequence of regulons to export.
    :param fname: The name of the file to create.
    """
    graph = nx.DiGraph()
    for regulon in regulons:
        src_name = regulon.transcription_factor
        graph.add_node(src_name, group='transcription_factor')
        edge_type = 'activating' if 'activating' in regulon.context else 'inhibiting'
        node_type = 'activated_target' if 'activating' in regulon.context else 'inhibited_target'
        for dst_name, edge_strength in regulon.gene2weight.items():
            graph.add_node(dst_name, group=node_type, **regulon.context)
            graph.add_edge(src_name, dst_name, weight=edge_strength, interaction=edge_type, **regulon.context)
    nx.readwrite.write_graphml(graph, fname)
