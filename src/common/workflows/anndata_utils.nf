// read the dataset id from an anndata file
process getDatasetId {
  container 'mumichae/scib-base:0.1'
  cache 'deep'

  input:
    file h5ad
  
  output:
    tuple stdout, file("$h5ad")

  script:
    """
    #!/usr/bin/env python
    import anndata
    ad = anndata.read_h5ad('$h5ad', backed=True)
    assert 'dataset_id' in ad.uns
    print(ad.uns['dataset_id'], end='')
    """
}

// read the method id from an anndata file
process getMethodId {
  container 'mumichae/scib-base:0.1'
  cache 'deep'
  
  input:
    file h5ad
  
  output:
    tuple stdout, file("$h5ad")

  script:
    """
    #!/usr/bin/env python
    import anndata
    ad = anndata.read_h5ad('$h5ad', backed=True)
    assert 'method_id' in ad.uns
    print(ad.uns['method_id'], end='')
    """
}