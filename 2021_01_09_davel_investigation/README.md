## Summary

This study compared the outcome of reading several molecules from SDF using OpenEye and RDKit ToolkitWrappers in the 0.8.X OpenFF toolkit. The OpenFF Molecule objects were compared after loading using several methods. Details are available in the notebook in this directory. 

One major outcome was the discovery that the [MiniDrugBank.sdf](https://github.com/openforcefield/openff-toolkit/blob/0.8.3/openforcefield/data/molecules/MiniDrugBank.sdf) file in the OpenFF toolkit has several ambiguously-defined or SDF/MOL specification-violating molecules. This is likely due to its conversion from mol2 to sdf in the past. One conclusion of this study is that a different molecule-loading test set should be found.

See also

https://github.com/openforcefield/openff-toolkit/issues/824
