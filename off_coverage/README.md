# off_coverage 1.0

NOTE: this project was previously hosted on [Sourcehut](https://hg.sr.ht/~dalke/off_coverage).

## Goal

The goal of this project is to identify a subset of compound records
which cover the same range of feature space as a larger data set.

## Basic concepts - feature space coverage

Each compound has an identifier. Each identifier is associated with
zero or more features. A feature is simply a label.

A feature label might mean "contains unspecified tetrahedral
chiralirty" or "does not contains unspecified chirarlity" or "causes
RDKit to complain about a valence error" or "OEChem and RDKit report a
different number of atoms" or "cannot be isomorphically matched to
itself".

Features may also come from code coverage, such as "lines {298, 400,
and 582} of toolkits.py were executed in any order", or code tracing,
such as "lines 298, 299, 300, and 299" of toolkits.py were executed in
exactly that order.

Given a set of compounds, each with features, can we select a subset
of the compounds such that all of the features are present in at least
compound? Can we further prioritize the result to minimize the number
of records, the number of features present in each one, the number of
atoms, etc.?

## The approach - Z3 

Yes. The
[set cover problem](https://en.wikipedia.org/wiki/Set_cover_problem)
is well-studied, with a number of available tools. This project uses
the
[Z3 Theorem Prover](https://en.wikipedia.org/wiki/Z3_Theorem_Prover).

The overall approach is to develop tools to generate features for a
given data set, convert that into something that Z3 can accept, and
let Z3 determine the (or at least a) coverage.

Z3 has a programming interface to many languages, including Python.
This project currently generates a file in "SMT" format, and includes
a decoder to parse the expected Z3 output into something easier to parse.


## Background - The driving problem


The Open Force Field Initiative uses
[MiniDrugBank](https://github.com/openforcefield/MiniDrugBank). It was
created to meet the needs for developing tools for chemical
perception. It contains selected structures from DrugBank Release
Version 5.0.1, and further processed and converted.

The 371 structures in MiniDrugBank.sdf are used during development to
check that new algorithms can handle a baseline set of chemical
structures. They also used in the unit tests to check for breaking
changes.

This data set may be both too large for some cases and too small for
others.

For example, the unit test take a while to run, in part because all
371 structures are tested in many test cases. Is there a smaller test
set which exercises the same functionality? This set might minimize
the number of records used, at the expense of more complicated
structures with many different feartures.


On the other hand, the baseline requirements have changed as OpenFF's
understanding of chemical nuances have evolved. MiniDrugBank doesn't
exercise some of the ways that other (non-Open Babel) toolkits
represent stereochemistry in 2D. There is a clear need to add new test
cases, but frequent processing all (appropriate) structures from
ChEMBL or PubChem is too time-intensive. Can these be processed once
to find good candidates?

Developers might prefer this test set to contain more compounds, each
with simpler features, to make it easier to diagnose problems.

## off_coverage.py 

The program `off_coverage.py` command-line tool implements several
subcommands to generate feature files using RDKit, OpenEye's OEChem,
and the OpenFF wrappers for those two toolkits, to merge feature
files, to convert a feature file to SMT format for Z3, and to parse
Z3's output to a simple list.


## The feature file format 

The feature file format is a text format with one line per record.

Each line contains an identifier followed by a tab followed by
space-separated features and/or feature weights.

(The initial tab delimiter is because some record identifiers may
contain a space character.)

A "feature" is a sequence of characters matching
`/[A-Za-z_][A-Za-z0-9_]+/`. These are the set of features that must be
present in the subset.

A "feature weight" is also a sequence of characters matching
`/[A-Za-z_][A-Za-z0-9_]+/`, followed by an `=`, followed by an
integer. These are used to generate weights during set coverage
minization.

## A first example

I'll demonstrate with the test file `small.smi` containing the
following three structures:

```
CN1C(=O)N(C)C(=O)C(N(C)C=N2)=C12 caffeine
C1C(C)=C(C=CC(C)=CC=CC(C)=CCO)C(C)(C)C1 vitamin a
O=C(C)Oc1ccccc1C(=O)O aspirin
```

I'll use the 'rdkit' command to have OpenFF use RDKit to the SMILES
file and generate an OpenFF topology using `from_file_obj()`.

```
% python off_coverage.py rdkit small.smi --atom-count
caffeine	rd_natoms=24 rd_parse_ok
vitamin a	rd_off_undef_stereo_bond rd_parse_err
aspirin	rd_natoms=21 rd_parse_ok
```

On the first output line, the `rd_parse_ok` means that OpenFF's
RDKitToolkitWrapper was able to parse the record into an RDKit
molecule and the `rd_natoms=21` means there were 24 atoms in total
(including hydrogens) in the OpenFF Molecule.

On the second output line, the `rd_off_undef_stereo_bond` means OpenFF
could not convert the molecule into its Molecule because there was an
undefined stereo bond.

I'll do it again, this time saving the results to a file:

```
% python off_coverage.py rdkit small.smi --atom-count -o small.feats
```

I'll then set up a Z3 solver using the features, find the optimal
cover, and show the results:

```
% python off_coverage.py minimize small.feats
Reading 'small.feats' ...
Creating Z3 model ...
Searching for a solution ...
Finished searching. Writing current model.
Solution: optimal Selected: 2/3
1	vitamin a
1	aspirin
```

or, without progress information and with a list of identifiers which
aren't included:

```
% python off_coverage.py minimize small.feats --quiet --all
Solution: optimal Selected: 2/3
0	caffeine
1	vitamin a
1	aspirin
```

## Changing the weight for each identifier

The default goal minimizes the number of identifiers used; each record
has a weight of 1. The weight value can be changed with `minimize
--weight`, as in the following:

```
% python off_coverage.py minimize small.feats --weight rd_natoms --quiet
 Unable to evaluate 'rd_natoms' for id 'vitamin a': name 'rd_natoms' is not defined
 Available names are: ['rd_off_undef_stereo_bond', 'rd_parse_err']
```

Okay, not *quite* that. The problem is `rd_natoms` isn't defined for
record "vitamin a". The value after `--weight` is a Python expression
evaluated for each entry in the feature file. Each feature and feature
weight is available as a Python variable with the same name. Features
have a value of True and feature weights have the corresponding
value (that is, the feature weight `rd_natoms=21` is available as the
Python variable name `rd_natoms` with integer value 21).

In addition, the special variable '_' is a dictionary mapping feature
or feature weight name to its value. This can be useful if your
feature name is a Python reserved word because you can write `_["in"]`
to get access to the feature named `in`.

This can also be used to test if a variable exists, using a dictionary
contains test (`key in dict`):

```
% python off_coverage.py minimize small.feats --weight 'rd_natoms if "rd_parse_ok" in _ else 0' --quiet
Solution: optimal Selected: 2/3
1	vitamin a
1	aspirin
```

or by using `dict.get()` to return a default value if the name doesn't exist:

```
% python off_coverage.py minimize small.feats --weight '_.get("rd_natoms", 0)' --quiet
Solution: optimal Selected: 2/3
1	vitamin a
1	aspirin
```

Both examples set up Z3 to weight each structure by `rd_natoms`, which
is a feature weight containing the number of atoms. The above will
minimize the total number of atoms, which won't change anything
because the default output coincidentally chose that solution.

I'll let's change things and minimize the negative weight:

```
% python off_coverage.py minimize small.feats --weight '-_.get("rd_natoms", 0)' -q
Solution: optimal Selected: 3/3
1	caffeine
1	vitamin a
1	aspirin
```

Oops! That doesn't work. This selects all of the molecules, because
that gives the largest negative number!

Try again, this time subtracting from a large number:

```
% python off_coverage.py minimize small.feats --weight '1000-_.get("rd_natoms", 0)' -q
Solution: optimal Selected: 2/3
1	caffeine
1	vitamin a
```

Yes, this is the correct solution.

Two other special variables are available, in addition to the `_`
mentioned earlier. The special variable `_features` is a set
containing only the feature names. For example, you can strongly
weight the result to prefer a small number of features per record by
using something like `len(_features)**2`. The special variable `ID`
stores the record id.

Note: the `--weight` expression must return an integer.

## RDKit feature generation

By default the "rdkit" command will:

* Convert the RDKit molecule into an OpenFF Molecule using
RDKitToolkitWrapper, which applies OpenFF's full normalization (add
hydrogens, perceive stereo, etc.)

* Skip any molecules with an atom containing element #0 (typically
wildcard atoms ("`*`") in SMILES or R-groups in an SDF).

* Parse the RDKit warning messages and include a feature for each
different type of message (the `--parse-warnings`) option

* Include the rd_natoms feature weight  (the `--atom-count` option)

There is experimental code to generate features based on circular
fingerprints (`--circular`). That does not seem useful as OpenFF
doesn't need test cases for all possible circular environments.

There is experimental code to accept a pattern definition file for
substructure searching (`--patterns`). This is more likely to be
useful in the future, but is probably broken and definitely not useful
now.

## OEChem feature generation

By default the "openeye" command will:

* Convert the OEGraphMol into an OpenFF Molecule using
OpenEyeToolkitWrapper, which applies OpenFF's full normalization (add
hydrogens, perceive stereo, etc.)

* Skip any molecules with an atom containing element #0 (typically
wildcard atoms ("`*`") in SMILES or R-groups in an SDF).

* Parse the OEChem warning messages and include a feature for each
different type of message (the `--parse-warnings`) option

* Include the oe_natoms feature weight  (the `--atom-count` option)

There is experimental code to generate features based on circular
fingerprints (`--circular`). That does not seem useful as OpenFF
doesn't need test cases for all possible circular environments.

There is experimental code to accept a pattern definition file for
substructure searching (`--patterns`). This is more likely to be
useful in the future, but is probably broken and definitely not useful
now.

## Common options

- If the input is an SD file then the default uses the record title as
the identifier. Use `--id-tag` to get the id from a tag in the SD record.

- By default only single component structures are processed. Use
`--multiple-components` to process multi-component structures.

- The `--novel` flag only outputs a new feature record if at least one
of the features hasn't been seen more than `N` times (by default, 1)
Use `--novel-count` to change the value of N. If `--novel-count` is
specied then `--novel` does not also need to be specified.

- Use `--quiet` to disable status information from the off_coverage.py
tool.

- Use `-o` or `--output` have the output go to the named file instead
of stdout.

- Use `--description` to save a mapping of feature name to a more
human-readable description to the named file.

- Use `--trace` to include code-coverage-based features.

- The `--ids` option is completely untested. The idea is to use or
exclude records if the id is specified in a given file.

There should probably be a way to filter molecules by atom count
and/or molecular weight and/or other features.

## Cross-compare feature generation

By default the `xcompare` subcommand will:

- Parse SMILES and SDF records and pass the record to the underlying
OpenFF toolkit wrappers. This gives an way to check if a toolkit
failed on a specific record, and record that failure.

- Compute atom and bond counts (`--counts`) using the OpenFF molecules
generated by each toolkit wrapper and set features if they match or
don't match. By default this will also add the `rd_natoms` and
`oe_natoms` weight features. This can be disabled with
`--no-atom-count`.

- Compare the two wrapper-generated OpenFF Molecules for isomorphism
and report success or failure (`--isomorphic`).

- Convert both OpenFF topology objects (if available) back to SMILES,
using both wrappers for each, and compare that the toolkit generate
canonically correct SMILES for the given wrapper (`--smiles`).

The one non-default feature is:

- Compere the InChI's generated from each wrapper-generated OpenFF
Molecule using the two toolkit wrappers (`--inchi`).

This doesn't seem to add more information than `--smiles` because most
of the issues came from converting the toplogy object to the given
chemistry toolkit molecule.

## xcompare using CHEMBL113 (caffeine), and `--description`

I'll use `xcompare` to analyze
`[CHEMBL113.sdf](https://www.ebi.ac.uk/chembl/compound_report_card/CHEMBL113/)`,
which is the ChEMBL record for caffeine. First, I'll download it from
ChEMBL. Annoyingly, their download doesn't include the "`$$$$`"
terminator for SDF records, so I'll have to add that to the downloaded
record:

```
% curl -O 'https://www.ebi.ac.uk/chembl/api/data/molecule/CHEMBL113.sdf'
% printf '\n$$$$\n' >> CHEMBL113.sdf
```

As another nuisance, ChEMBL doesn't store a usable id in the SDF
record's title line, the feature files need a usable id, and
`off_coverage` raises an exception if that happens:


```
% python off_coverage.py xcompare CHEMBL113.sdf
Traceback (most recent call last):
  File "off_coverage.py", line 3435, in <module>
    main()
  File "off_coverage.py", line 3432, in main
    parsed_args.command(parsed_args.subparser, parsed_args)
  File "off_coverage.py", line 2922, in xcompare_command
    raise ValueError(f"Missing id for record #{recno}")
ValueError: Missing id for record #1
```

But that's okay, the `off_coverage` commands support `--id-tag` to get
the record name from the named tag.

```
% python off_coverage.py xcompare CHEMBL113.sdf --id-tag chembl_id | fold -s
CHEMBL113	oe_natoms=24 oe_parse_ok rd_natoms=24 rd_parse_ok
xcmp_isomorphic_ok xcmp_natom_ok xcmp_nbond_ok xcmp_oe2oe_rd2oe_smi_ok
xcmp_oe2oe_rd2rd_smi_natoms_ok xcmp_oe2oe_smi_ok xcmp_oe2rd_rd2rd_smi_ok
xcmp_oe2rd_smi_ok xcmp_rd2oe_smi_ok xcmp_rd2rd_smi_ok
```

That's rather a lot of opaque feature names. What do they mean?

You can get a mapping from the feature name to its full description
using the `--description` option, which saves those unique
descriptions to the named file.

```
% python off_coverage.py xcompare CHEMBL113.sdf --id-tag chembl_id \
        --description caffeine.txt > /dev/null
% cat caffeine.txt
oe_parse_ok	OpenEyeToolkitWrapper could parse the record
rd_parse_ok	RDKitToolkitWrapper could parse the record
xcmp_natom_ok	OpenEye and RDKit wrappers have the same atom counts
xcmp_nbond_ok	OpenEye and RDKit wrappers have the same bond counts
xcmp_isomorphic_ok	OE and RD topologies are isomorphic
xcmp_oe2oe_smi_ok	Can convert OE topology to SMILES using OE wrapper
xcmp_rd2rd_smi_ok	Can convert RD topology to SMILES using RD wrapper
xcmp_oe2oe_rd2rd_smi_natoms_ok	RD and OE same-toolkit SMILES have the same atom counts
xcmp_oe2rd_smi_ok	Can convert OE topology to SMILES using RD
xcmp_rd2oe_smi_ok	Can convert RD topology to SMILES using RD
xcmp_oe2oe_rd2oe_smi_ok	OE generated SMILES from RD and OE are the same
xcmp_oe2rd_rd2rd_smi_ok	RD generated SMILES from RD and OE are the same
```

(Note that there's a tab between the feature name and its description.)

## Code coverage feature generation with `--trace`

All three of these feature generation subcommands support the
`--trace` command to include a feature based on which lines of code
are executed in the `openff` module or any of its submodules.

The above test conditions treat the code as a black box as it only
looks at return values, exceptions, and warning/error messages. It
doesn't know which code is actually being executed. There may be some
important code which is silently excuted by some molecules but not
others. How do we detect that?

Python's
[`sys.settrace`](https://docs.python.org/3/library/sys.html?highlight=settrace#sys.settrace)
function lets debuggers and coverage tools figure out which lines (and byte code) are being executed.

The `--trace` command uses that hook to record which lines from the
toolkits module is being executed for each calculation function, and
adds a unique execution signature to the features using the `trace_`
prefix, as in:

```
% python off_coverage.py xcompare small.smi -o small.feats --trace on
% awk '{print substr($0, 1, 500)}' small.feats | fold -b -w 70
caffeine	oe_natoms=24 oe_parse_ok rd_natoms=24 rd_parse_ok trace_00d69
73c4b864113f064f69bf34c85c7c098d003019ab009607997e7298de0b9 trace_0115
07fcb297a7a2d89f14bc2eb3dd3e733380b938dbea3d9ffa37f797c015db trace_012
14d40909181e3d1796e5dbdd55d569bda9ceb6fb6d9bf354969805154d065 trace_04
02e5367624eddc8be80d9958c92000ad55fe29b067318b3ca16f0295842165 trace_0
7a22e9d04590ecb2e5ad642bf7d6775fa0f6e44685e8023c2cdbd3a6a5eaa47 trace_
0888b65e011891915642e4d72654bf0584605e4b5d59147c3a3d88cb759ac41f trace
_090d48ea2
vitamin a	oe_off_undef_stereo_bond oe_parse_err rd_off_undef_stereo_bo
nd rd_parse_err trace_06581a188d3aa9679f7f20506e9c1f37415a05ba21e24e80
980a2a596e5e3df0 trace_08b38952442bfdd7c7a35fedc5a2dd9b79e31fb3c5169dd
44efc1496904ec8d9 trace_2ea6f5fe593ffdf636ad2a2c62f5d77673e25f87bd36fe
fda99c5b1f106a9f2e trace_4d3399336374d10b6bdd30383ea08549c059497b65ec9
ecafb4e31b9bc913fb3 trace_55aa2222bb49d36033cb6afcc9aef2efcf64e1ccde47
714e8d2fb1002e3b43eb trace_691e1289e9d2965305025e67871008ef0f1f29c1d9b
f313c0ce6c
aspirin	oe_natoms=21 oe_parse_ok rd_natoms=21 rd_parse_ok trace_01214d
40909181e3d1796e5dbdd55d569bda9ceb6fb6d9bf354969805154d065 trace_0402e
5367624eddc8be80d9958c92000ad55fe29b067318b3ca16f0295842165 trace_07a2
2e9d04590ecb2e5ad642bf7d6775fa0f6e44685e8023c2cdbd3a6a5eaa47 trace_088
8b65e011891915642e4d72654bf0584605e4b5d59147c3a3d88cb759ac41f trace_0b
f87c3627e7686331e275ad3477bbe16341eb580e7dda6d8e0ebe9cba5be648 trace_0
ddb2ffa35d246898ca0f4ff428776fe6c65f464a6c5cacbb1463e2c0268280e trace_
0ff775f6d7
```

There are a *lot* of signatures that I didn't show. But, what are
they?

There are a several ways to generate a signature.

Perhaps the simplest way to think about is to track which line numbers
are executed for each module (`--trace mod-cover`). The resulting set
of line numbers gives the line coverage for the module. These can be
sorted and combined with the module name to give a canonical
representation for that module's code coverage, then converted to a
unique signature (in this case, as a sha256 hash).

Another is to track the exact sequence of line numbers executed in the
module (`--trace mod-sequence`). This much stricter definition of
equivalence doesn't seem useful.

A third is to track the pairs of line numbers (`--trace mod-pairs`),
which captures some of the branching information not present in simple
line coverage. (This was influenced by the
[American Fuzzy Lop](https://lcamtuf.coredump.cx/afl/) fuzzer.)

Suppose, however, that a given function is called several times.
Module-level coverage merges them into one when they might have
distinct execution patterns.

Instead, another available option is to track each function call
separately. These are available as `func-cover`, `func-sequence` and
`func-pairs`.

The `--trace on` option is an alias for `--trace func-pairs`.

### trace descriptions

Each trace feature also generates a description to the `--description`
file, so we can out what, say, something like:

```
trace_0402e5367624eddc8be80d9958c92000ad55fe29b067318b3ca16f0295842165 
```

means.

```
% grep trace_0402e53676 small.txt  | fold
trace_0402e5367624eddc8be80d9958c92000ad55fe29b067318b3ca16f0295842165	<mod: op
enff.toolkit.topology.molecule func: _add_bond<pairs: [(0,3932),(3932,3933),(39
33,3934),(3934,3946),(3946,3950),(3950,3951),(3951,3952),(3952,3953),(3953,3954)
,(3954,3955),(3955,3956),(3956,3958),(3958,3959),(3959,3961),(3961, -7923)]>>
```

This says it's a trace of the `_add_bond` function in the module
`openff.toolkit.topology.molecule`, using the "pairs" trace method.

The unique pairs for this function are given as line numbers, with a
couple of exceptions. The `0` in `(0,3932)` indicates this is the
first line executed in the function. The `-7923` in the last tuple
indicates the function ended with a return from line `3961`.

Why is it negative? There are two ways for a function to end; with a
`return` or with a `raise`. Both of these have a line number
(`lineno`). With a `return`, the value `-2*lineno-1` is used. With
a `raise`, the value `-2*lineno-2`is used.

With tedious examination it's possible to trace some of the execution
order. I think this really needs a GUI tool for visualization.

The other trace methods use different representations. I hope it's
easy to figure out how it works. You'll have to look at the source
code for full details.

### Trace coverage sizes

The different trace methods generate different sets of unique trace
features. Here is a summary of the selection size for a xcompare of
MiniDrugBank.sdf (see below):

- `mod-sequence` selected: 370/371 
- `mod-pairs` selected: 74/371 
- `mod-cover` selected: 60/371
- `func-sequence` selected: 354/371
- `func-pairs` selected: 35/371
- `func-cover` selected: 29/371

In general, the module-level trace methods require greater code path
similarity than the function-level methods; and `sequence` requires
greater code path similarity than `pairs`, which requires greater code
path similarity than `cover`.


## Putting it all together

I used the default feature generation `xcompare` to generate features
for MiniDrugBank:

```
% python off_coverage.py xcompare --trace on MiniDrugBank.sdf -o MiniDrugBank.feats
```

That generates output like the following two lines (severely trimmed for readability):

```
% head -2 MiniDrugBank.feats | awk '{print substr($0, 1, 900)}' | fold -w 75 -s
DrugBank_5354	oe_natoms=44 oe_parse_ok rd_natoms=44 rd_parse_ok
trace_00d6973c4b864113f064f69bf34c85c7c098d003019ab009607997e7298de0b9
trace_011507fcb297a7a2d89f14bc2eb3dd3e733380b938dbea3d9ffa37f797c015db
trace_01214d40909181e3d1796e5dbdd55d569bda9ceb6fb6d9bf354969805154d065
trace_0402e5367624eddc8be80d9958c92000ad55fe29b067318b3ca16f0295842165
trace_0888b65e011891915642e4d72654bf0584605e4b5d59147c3a3d88cb759ac41f
trace_090d48ea2e1d45ee01639c83431e11dd7125ed4afe510664112c6ed101421dd4
trace_0bd8730eb2c3f36e0174dbea6febe48c15f9614844f037c298629b92cac113cf
trace_0d9c32972dc842798a13a9d0c69acfc082a46dfab9dbc4eeb246074dc067ee2c
trace_0ddb2ffa35d246898ca0f4ff428776fe6c65f464a6c5cacbb1463e2c0268280e
trace_0ff775f6d7cbcb51675394b12592c743aed60517ca4f224f7fab0c520ac62bb5
trace_11c00f044b4643a8d80afd72cd61fd5d38ede45b9e3ce9ed80ec961b736f20f8
trace_1548bedae311d53ef28c2d128fa29fa7fc6ee7cd2bdef8587
DrugBank_2791	oe_natoms=19 oe_parse_ok rd_natoms=19 rd_parse_ok
trace_01214d40909181e3d1796e5dbdd55d569bda9ceb6fb6d9bf354969805154d065
trace_0402e5367624eddc8be80d9958c92000ad55fe29b067318b3ca16f0295842165
trace_0888b65e011891915642e4d72654bf0584605e4b5d59147c3a3d88cb759ac41f
trace_0ddb2ffa35d246898ca0f4ff428776fe6c65f464a6c5cacbb1463e2c0268280e
trace_0ff775f6d7cbcb51675394b12592c743aed60517ca4f224f7fab0c520ac62bb5
trace_1548bedae311d53ef28c2d128fa29fa7fc6ee7cd2bdef858756a775e31344710
trace_171eced97b9e1d829dd54162dbc7115f8ff0f122cd8cba73c2e97c36c3c7d429
trace_196db8571afc1e27581147b1547db4f6ec9223a13c2d699e12cd990df7649257
trace_1af8befd8ee92bbe29643bcdf4a2ef4febae1ae13402fadc794c24cf09694350
trace_1be3393c06f6db9a4a08aec57e1bf28e4720b762c316d3f2878ebf26bbc91e16
trace_1f4309b6c6e16104688d0960bcdfcb668e3f77d833e0a51a8126652a9ed69504
trace_22f8f57ec2cfd2531e0933fe1f2f7f9cd985da01f0cacd704
 ```

I used Z3 to minimize the result, which it did in about four seconds:

```
% python off_coverage.py minimize MiniDrugBank.feats -q
Solution: optimal Selected: 35/371
1	DrugBank_5354
1	DrugBank_5387
1	DrugBank_246
1	DrugBank_423
1	DrugBank_3479
1	DrugBank_3547
1	DrugBank_3581
1	DrugBank_3693
1	DrugBank_6182
1	DrugBank_6295
1	DrugBank_4161
1	DrugBank_1538
1	DrugBank_6531
1	DrugBank_4249
1	DrugBank_1700
1	DrugBank_4323
1	DrugBank_4346
1	DrugBank_6775
1	DrugBank_4468
1	DrugBank_1849
1	DrugBank_4593
1	DrugBank_4662
1	DrugBank_7108
1	DrugBank_2095
1	DrugBank_2148
1	DrugBank_2178
1	DrugBank_2186
1	DrugBank_2210
1	DrugBank_4959
1	DrugBank_2429
1	DrugBank_2538
1	DrugBank_2563
1	DrugBank_2642
1	DrugBank_2684
1	DrugBank_2728
```

In hand-wavy theory this means the test set can go from 371 structures
to only 35!


## first line of `minimize` output, and --timeout

If the Z3 `minimize` step finishes then the first line of the output
will look something like:

```
Solution: optimal Selected: 72/371
```

The second word, "optimal", indicates Z3 found an optimal solution.
The "72" says that 72 records are needed, out of 371 input record.

It can take Z3 a long time to find a minimal set. You can interrupt it
with ^C, or use the `--timeout` option to set a maximum search time,
in milliseconds.

In either case, off_coverage will still attempt to report a set of
features to use.

If Z3 has a solution, but not one that it can prove is the minimal
solution, then `off_coverage` will print the word "approximate"
instead of "optimal".

Sometimes Z3 doesn't even have valid solution. In that case,
"off_coverage" will fill in the missing parts using a simple greedy
algorithm: for each missing feature it find the first record
containing that feature, and includes it.

If that happens then the second word will be "greedy" instead of
"optimal" or "approximate". Here's an example using a subset of the
xcompare features from ChEBI, where I first give it 100 milliseconds
and then give it 30 seconds:

```
% python off_coverage.py minimize chebi.feats --timeout 100 | head -1
Reading 'chebi.feats' ...
Creating Z3 model ...
Searching for a solution ...
Finished searching. Writing current model.
Solution: greedy Selected: 724/28793
% time python off_coverage.py minimize chebi.feats --timeout 30000 | head -1
Reading 'chebi.feats' ...
Creating Z3 model ...
Searching for a solution ...
Finished searching. Writing current model.
Solution: optimal Selected: 629/28793
```

Several subcommands support a `--novel` option, which implements the
"greedy" algorithm including the `merge` subcommand. I'll use that to
show the 100ms timeout is too short:

```
% python off_coverage.py merge chebi.feats chebi.feats --novel -q | wc -l
     724
```

(The `merge` command requires two feature files, so I specified
chebi.feats twice.)

On the other hand, a 5 minute timeout is enough to find that 629
records is the minimal size.

## merge command

The `merge` command merges two feature datasets into one. There are
several different ways to merge the two data set. See its `--help` for
details.

In this example, suppose I want to replace MiniDrugBank records with
ChEBI records, if possible, but keep essential MiniDrugBank records.

First, I'll merge merge the two datasets ("MiniDrugBank.feats" and
"chebi.feats"), and only use those features from chebi.feats which are
also features in MiniDrugBank.feats:


```
% python off_coverage.py merge MiniDrugBank.feats chebi.feats -o merge.feats --no-new-features
```

I'll then minimize using a weight of 10 for ids starting with
"DrugBank" and 1 for all other records.

```
% python off_coverage.py minimize merge.feats --weight '10 if ID.startswith("DrugBank") else 1' -q
```

This found that 10 of the MiniDrugBank records in the 72 molecule
subset could be replaced by ChEBI records. (But I'll need to examine
things more closely as `--no-new-features` is barely tested.)

## subtract command

This tool was used to
[reduce the number of test cases in OpenFF](https://github.com/openforcefield/openff-toolkit/pull/999)
based on code coverage analysis.

Sometimes this could be done by generating trace features and finding
a minimal solution. Other cases required a differential analysis,
where I wanted the coverage for Y, but computing Y also required
computing X. My approach was generate the coverage for X, generate the
coverage for X+Y, then `subtract` the former features from the latter,
leaving only the unique features for Y.

I'll use that to see which functions the `--inchi` tests adds to the
`xcompare` analysis.

```
% python off_coverage.py xcompare small.smi --trace on -o small_default.feats
% python off_coverage.py xcompare small.smi --trace on -o small_default_inchi.feats--default --inchi
% python off_coverage.py subtract small_default_inchi.feats small_default.feats -o small_inchi.feats
```

The resulting `small_inchi.feats` file is significantly smaller than
the other two feature files:

```
% wc small_default.feats
       3     287   18420 small_default.feats
% wc small_default_inchi.feats
       3     297   18824 small_default_inchi.feats
% wc small_inchi.feats
       3      18     484 small_inchi.feats
```

It's so small that I'll display the (folded) content here:

```
% fold -s small_inchi.feats
caffeine	oe_natoms=24 rd_natoms=24
trace_098e8177c862aa994d8062aabb1f4c7204b6e6f6cd79a21c9285f5e52cbb9ff1
trace_e0c125114890b6370fcf3bb75a6d9e3de010d2768b50a43a763a87309bef9196
xcmp_oe2oe_inchi_ok xcmp_oe_rd_inchi_ok xcmp_rd2rd_inchi_ok
vitamin a
aspirin	oe_natoms=21 rd_natoms=21
trace_098e8177c862aa994d8062aabb1f4c7204b6e6f6cd79a21c9285f5e52cbb9ff1
trace_e0c125114890b6370fcf3bb75a6d9e3de010d2768b50a43a763a87309bef9196
xcmp_oe2oe_inchi_ok xcmp_oe_rd_inchi_ok xcmp_rd2rd_inchi_ok
```

There are only two new function calls! I'll re-do the `--default
--inchi` analysis to also include the code description, and see what
those functions are:

```
% python off_coverage.py xcompare small.smi --trace on -o small_default_inchi.feats \
     --default --inchi --description small_default_inchi.description
% egrep '098e8177c862|e0c125114890' small_default_inchi.description
trace_e0c125114890b6370fcf3bb75a6d9e3de010d2768b50a43a763a87309bef9196	<mod: openff.toolkit.utils.rdkit_wrapper func: to_inchi <pairs: [(0,1826),(1826,1828),(1828,1829),(1829,1832),(1832,1833),(1833,-3667)]>>
trace_098e8177c862aa994d8062aabb1f4c7204b6e6f6cd79a21c9285f5e52cbb9ff1	<mod: openff.toolkit.utils.openeye_wrapper func: to_inchi <pairs: [(0,1385),(1385,1387),(1387,1389),(1389,1395),(1395,1397),(1397,-2795)]>>
```

These make identical function calls, and inspection shows the
`--inchi` tests add nothing beyond detecting failures in the
underlying toolkit, so we can reduce the entire test suite to a single
case:

```
% python off_coverage.py minimize small_inchi.feats -q
Solution: optimal Selected: 1/3
1	caffeine
```

## Programmatic interface

If you put `off_coverage.py` in an importable directory, or install it
via the included `setup.py`, then you can use some of the tracing code
in your own software.

Be aware that it was not developed as a tool for general use. If it
works outside of the OpenFF context it was developed for - cool!


What code paths does Python's
[urlparse](https://docs.python.org/3/library/urllib.parse.html#urllib.parse.urlparse)
use to parse different URLs? The following uses off_coverage for a
coverage analysis of a list of URLs. It sets up a CallCoverage
instance, which can create a context manager to be report coverage for
a named test case.

```
import off_coverage

from urllib.parse import urlparse

urls = [
    "https://openforcefield.org/",
    "https://github.com/openforcefield/openff-toolkit/tree/master/examples#index-of-provided-examples",
    "http://dalkescientific.com/",
    "http://localhost:8080",
    "https://chemfp.com",
    "ftp://ftp.ncbi.nlm.nih.gov/pubchem/",
    "https://www.ebi.ac.uk/chembl/",
    "https://www.google.com/search?q=chemfp",
    "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/",
]

# Record coverage analysis to "urlparse.feats".
# There will be one record per line.
# The default uses trace_type="func-pairs".
coverage = off_coverage.open_call_coverage(
    "urlparse.feats",
    modules = [
        "urllib.",  # limit coverage to modules under urllib
        ]
        )

for url in urls:
    # Start tracing.
    # Use the url as the record name.
    with coverage.record(url) as state:
        # Add a 'n' feature weight using the URL length
        state.add_feature(f"n={len(url)}")
        # Parse the url.
        urlparse(url)
```

I'll run it the program generate the coverage information:

```
% python url_check.py
```

This generated the output file `urlparse.feats` and the description
output file `urlparse.feats.description`. I'll use the feature file to
generate a minimal set:

```
% python off_coverage.py minimize urlparse.feats -q
Solution: optimal Selected: 4/9
1	https://github.com/openforcefield/openff-toolkit/tree/master/examples#index-of-provided-examples
1	http://localhost:8080
1	ftp://ftp.ncbi.nlm.nih.gov/pubchem/
1	https://www.google.com/search?q=chemfp
```

To help figure out what's going on, here's a program to show which
URLs trigger the different code paths:

```
from collections import defaultdict
d = defaultdict(list)
for line in open("urlparse.feats"):
    id, _, rest = line.rstrip("\n").partition("\t")
    for feature in rest.split():
        if feature.startswith("trace_"):
            d[feature].append(id)
for trace, ids in sorted(d.items()):
    print(trace)
    for id in ids:
        print("  ", id)
```

This generates:

```
trace_00506cf40fe57c67f606f344c6ddccdbe14708c32a323fabe05adb9bb101e4ea
   https://github.com/openforcefield/openff-toolkit/tree/master/examples#index-of-provided-examples
trace_0a3380a6bacc6f08bcf5fc06eb75c84847760a9a8058a5d823065f807a23967b
   https://openforcefield.org/
   https://github.com/openforcefield/openff-toolkit/tree/master/examples#index-of-provided-examples
   http://dalkescientific.com/
   http://localhost:8080
   https://chemfp.com
   ftp://ftp.ncbi.nlm.nih.gov/pubchem/
   https://www.ebi.ac.uk/chembl/
   https://www.google.com/search?q=chemfp
   ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/
trace_2f7a951665e32c43625f0c07b81bbed5bfed4e77f5ed03d0e0004220263193ce
   https://openforcefield.org/
   https://github.com/openforcefield/openff-toolkit/tree/master/examples#index-of-provided-examples
   https://chemfp.com
   ftp://ftp.ncbi.nlm.nih.gov/pubchem/
   https://www.ebi.ac.uk/chembl/
   https://www.google.com/search?q=chemfp
   ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/
trace_42802f38c3da956be66675d79c10c2e66d1c7bd560c9b6e46bde8f826672140a
   http://dalkescientific.com/
   http://localhost:8080
trace_43c7927193fe14c0dae2ab8cd84e6bc18e33c36d21e1e169bf9b27c36a73a0ec
https://www.google.com/search?q=chemfp
```

What is "trace_00506cf...", which is only visited by the GitHub URL?

```
% grep trace_00506cf urlparse.feats.description | fold
trace_00506cf40fe57c67f606f344c6ddccdbe14708c32a323fabe05adb9bb101e4ea	<mod: ur
llib.parse func: urlsplit <pairs: [(0,418),(418,419),(419,420),(420,421),(421,42
2),(422,424),(424,426),(426,427),(427,428),(428,429),(429,444),(444,445),(444,45
0),(445,444),(450,451),(451,453),(453,455),(455,456),(456,457),(457,458),(458,46
0), (460,461),(461,462),(462,464),(464,465),(465,466),(466,467),(467,-935)]>>
```

More analysis shows the unique branch points are:

```
(460, 461), (461, 462)
```

This corresponds to the following lines from `urllib.parse`:

```
460:   if allow_fragments and '#' in url:
461:        url, fragment = url.split('#', 1)
```

Doing this manually was a pain. This is where a visualization tool
would be helpful.

Anyway, there are multiple solutions. lf I prefer one with smaller URL
lengths, I'll ask to minize also by `n`, which was a feature weight I
added in `url_check.py`:

```
% python off_coverage.py minimize urlparse.feats -q --weight n
Solution: optimal Selected: 4/9
1	https://github.com/openforcefield/openff-toolkit/tree/master/examples#index-of-provided-examples
1	http://localhost:8080
1	https://chemfp.com
1	https://www.google.com/search?q=chemfp
```

Or, if I really wanted ftp URLS, I can weight that one post hoc by
looking for "ftp" in the special "ID" variable:

```
% python off_coverage.py minimize urlparse.feats -q --weight '("ftp" not in ID)*100'
Solution: optimal Selected: 4/9
1	https://github.com/openforcefield/openff-toolkit/tree/master/examples#index-of-provided-examples
1	http://localhost:8080
1	ftp://ftp.ncbi.nlm.nih.gov/pubchem/
1	https://www.google.com/search?q=chemfp
```



## Legal

This program was written by Andrew Dalke <dalke@dalkescientific.com>,
under contract for the Open Force Field Initiative.

### MIT LICENSED

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Changes

### Version 1.0 - 9 July 2021

- change the structure input routines to always use OpenFF to parse
the structure records

- added `--description` option to get a mapping from feature names to
a more readable description.

- new subcommand to `subtract` one set of features from another

- added programmatic access to the call coverage utility

- improved error handling for ChEBI test cases

- disallow inputs with atomic number of 0 (eg, `*` in SMILES, or
R-groups in SDF records).

### Version 0.9 - 8 June 2021

- fixed `--circular` bug reported by Rajarshi Guha

- changed the function trace implementation so each call of the same function is traced independently

- added option to trace pairs of line numbers, which should better approximate branch coverage

- added `--no-new-features` to the `merge` command

- switched to Z3's Python API instead of writing SMT output, having
the user manually run z3, and decoding the result. The new `minimize`
command replaces the old `z3_create` and `z3_decode` commands.

- The minimize `--weight` expression can now access the record id as
`ID`.

### Version 0.8 - 4 June 2021

Initial release

