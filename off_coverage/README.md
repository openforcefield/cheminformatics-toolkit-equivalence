# off_coverage 0.9

NOTE: this is the first release at [Open Force Field's GitHub repo](https://github.com/openforcefield/cheminformatics-toolkit-equivalence) insted of a [personal Sourcehut repo](https://hg.sr.ht/~dalke/off_coverage).


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

I'll use the 'rdkit' command to generate atom counts and include
features which came from converting the RDMol into an OpenFF topology
using `from_object()`:

```
% python off_coverage.py rdkit small.smi --atom-count --openff 
caffeine	rd_natoms=24 rd_openff_good rd_parse
vitamin a	rd_natoms=48 rd_off_undef_stereo_bond rd_parse
aspirin	rd_natoms=21 rd_openff_good rd_parse
```

On the first output line, the `rd_parse` means that RDKit was able to
parse the record into an RDKit molecule, the `rd_openff_good` means
that OpenFF was able to convert the molecule into a topology object,
and the `rd_natoms=21` means there were 21 atoms in total (including
hydrogens).

On the second output line, the `rd_off_undef_stereo_bond` means OpenFF
could not convert the molecule into topology object because there was
an undefined stereo bond.

I'll do it again, this time saving the results to a file:

```
% python off_coverage.py rdkit small.smi --openff --atom-count -o small.feats
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
1	vitamin a
1	aspirin
```

This example sets up Z3 to weight each structure by `rd_natoms`, which
is a feature weight containing the number of atoms. The above will
minimize the total number of atoms, which won't change anything
because the default output coincidentally chose that solution.

I'll let's change things and minimize the negative weight:

```
% python off_coverage.py minimize small.feats --weight -rd_natoms
usage: off_coverage.py z3_create [-h] [--weight EXPR] [--output OUTPUT]
                                 filename
off_coverage.py minimize: error: argument --weight: expected one argument
```

Ooops! That didn't work. The argument parser (Python's argparse)
thinks `--rd_natoms` is a command-line option, not a parameter to
weight.

But that's okay - the function actually takes a Python expression, so I'll subtract from 0:

```
% python off_coverage.py minimize small.feats --weight 0-rd_natoms --quiet
Solution: optimal Selected: 3/3
1	caffeine
1	vitamin a
1	aspirin
```

Oops! That still doesn't work. This selects all of the molecules,
because that gives the largest negative number!

Try again, this time subtracting from a large number:

```
% python off_coverage.py minimize small.feats --weight 10000-rd_natoms -q
Solution: optimal Selected: 2/3
1	caffeine
1	vitamin a
```

Yes, this is the correct solution.

The weight expression can use any of the defined feature weights. Each
feature name is also available (as True). This is useful to have
default for missing cases, for example: `1000 if rd_parse_fail else
rd_natoms` will use the value of 1000 if the `rd_parse_fail` feature
is present, presumably indicating that RDKit was unable to parse the
record.

Alternatively, the special variable `_` refers to the dictionary of
features and feature weights, so you could do `rd_natoms if
"rd_natoms" in _ else 1000`.

The special variable `_features` is a set containing only the feature
names. For example, you can strongly weight the result to prefer a
small number of features per record by using something like
`len(_features)**2`. The special variable `ID` stores the record id.

Note: the weight expression must return an integer.

## RDKit feature generation

By default the "rdkit" command will:

* Parse the RDKit warning messages and include a feature for each
different type of message (the `--parse-warnings`) option

* Include the rd_natoms feature weight  (the `--atom-count` option)

* Convert the RDKit molecule into an OpenFF topology (the `--openff`
option)

There is experimental code to generate features based on circular
fingerprints (`--circular`). That does not seem useful as OpenFF
doesn't need test cases for all possible circular environments.

There is experimental code to accept a pattern definition file for
substructure searching (`--patterns`). This is more likely to be
useful in the future, but is probably broken and definitely not useful
now.

Note: The rdkit standardization currently only does:

```
def rd_standardize_mol(mol):
    return Chem.AddHs(mol)
```

Future work should make it be more like the OpenFF standardizer.

## OEChem feature generation

By default the "openeye" command will:

* Parse the OEChem warning messages and include a feature for each
different type of message (the `--parse-warnings`) option

* Include the oe_natoms feature weight  (the `--atom-count` option)

* Convert the OEGraphMol into an OpenFF topology (the `--openff`
option)

There is experimental code to generate features based on circular
fingerprints (`--circular`). That does not seem useful as OpenFF
doesn't need test cases for all possible circular environments.

There is experimental code to accept a pattern definition file for
substructure searching (`--patterns`). This is more likely to be
useful in the future, but is probably broken and definitely not useful
now.

Note: The OEChem standardization currently only does:

```
def oe_standarize_mol(mol):
    if not oechem.OEAddExplicitHydrogens(mol):
        if oe_msg_handler is not None:
            oe_msg_handler.messages.append(
                "WARNING: Unable to add explicit hydrogens"
                )
    oechem.OEAssignAromaticFlags(mol, oechem.OEAroModel_MDL)
    oechem.OEPerceiveChiral(mol)
    return mol
```

Future work should make it be more like the OpenFF standardizer.

## Common options

- If the input is an SD file then the default uses the record title as
the identifier. Use `--id-tag` to get the id from a tag in the SD record.

- By default only single component structures are processed. Use
`--multiple-components` to process multi-component structures.

- The `--novel` flag only outputs a new feature record if at least one
of the features hasn't been seen more than `N` times (by default, 1)
Use `--novel-count` to change the value of N.

- Use `--quiet` to disable status information from the off_coverage.py
tool.

- Use `-o` or `--output` have the output go to the named file instead
of stdout.


- The `--ids` option is completely untested. The idea is to use or
exclude records if the id is specified in a given file.

There should probably be a way to filter molecules by atom count
and/or molecular weight and/or other features.

## Cross-compare feature generation

By default the `xcompare` subcommand will:

- Parse SMILES and SDF records and pass the record to the underlying
chemistry toolkits. This gives an way to check if a toolkit failed on
a specific record, and record that failure.

- Compute atom and bond counts (`--counts`) using both toolkits and
set features if they match or don't match. By default this will also
add the `rd_natoms` and `oe_natoms` weight features. This can be
disabled with `--no-atom-count`.

- Generate an OpenFF topology for both RDKit and OpenEye toolkit
wrappers, and record success or failure. If successful, compare the
two for isomorphism and report success or failure (`--isomorphic`).

- Convert both OpenFF topology objects (if available) back to SMILES,
using both wrappers for each, and compare that the toolkit generate
canonically correct SMILES for the given wrapper (`--smiles`).

- There is an incomplete `--inchi` option which does the same with
InChI. This didn't seem to add more information than `--smiles`
because most of the issues came from converting the toplogy object to
the given chemistry toolkit molecule.

## Code coverage feature generation with `--trace`

All three of these feature generation subcommand support the `--trace`
command to include a feature based on which lines of code are executed
in the `openff.toolkit.utils.toolkits` module.

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
% fold -b -w 70 small.feats
caffeine	oe_natoms=24 oe_openff_good oe_parse parse rd_natoms=24 rd_op
enff_good rd_parse trace_ccf31b724b03e766be0d42cc975fda4679a3cdd79c2b6
2e768f1889208913253 trace_de76635cb81cc901b30cc8d69f4c1faa3d9e9a0d6104
90caddcbeda8f5a8bf3b xcmp_isomorphic_err xcmp_natom_ok xcmp_nbond_ok x
cmp_smi_oe2_ok xcmp_smi_oe_ok xcmp_smi_oe_rd_ok xcmp_smi_ok xcmp_smi_r
d2_ok xcmp_smi_rd_oe_ok xcmp_smi_rd_ok
vitamin a	oe_natoms=48 oe_off_undef_stereo_bond oe_parse parse rd_nato
ms=48 rd_off_undef_stereo_bond rd_parse trace_1d70229208662e766d61f9a9
cbcf34074f11751b4b1efb19e93a2584bc9a3b55 xcmp_natom_ok xcmp_nbond_ok
aspirin	oe_natoms=21 oe_openff_good oe_parse parse rd_natoms=21 rd_ope
nff_good rd_parse trace_7644c766cd90ffd5561c41d558839a3636f70b91189f5d
b129028be298c39fcb trace_de76635cb81cc901b30cc8d69f4c1faa3d9e9a0d61049
0caddcbeda8f5a8bf3b xcmp_isomorphic_ok xcmp_natom_ok xcmp_nbond_ok xcm
p_smi_oe2_ok xcmp_smi_oe_ok xcmp_smi_oe_rd_ok xcmp_smi_ok xcmp_smi_rd2
_ok xcmp_smi_rd_oe_ok xcmp_smi_rd_ok
```

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

These can result in different sets of unique trace features. Here is a
summary of the selection size for a xcompare of MiniDrugBank.sdf (see
below):

- `mod-cover` selected: 72/371 
- `mod-sequence` selected: 342/371
- `mod-pairs` selected: 72/371
- `func-cover` selected: 31/371 
- `func-sequence` selected: 336/371
- `func-pairs` selected: 27/371

The option `--trace on` is an alias for `--trace mod-cover`. The
default is to not include tracing, which is the same as the option
`--trace off`.

I do not know if `func-pairs` is a better option than `mod-pairs`. I
do not know why `func-pairs` selects fewer cases then `func-cover`
since I think that should never happen.

Well hey, something for future investigation.

## Putting it all together

I used the default feature generation `xcompare` to generate features
for MiniDrugBank:

```
% python off_coverage.py xcompare --trace on MiniDrugBank.sdf -o MiniDrugBank.feats
```

That generates output like the following two lines (folded for readability):

```
% head -2 MiniDrugBank.feats | fold -w 75 -s
DrugBank_5354	oe_natoms=44 oe_openff_good oe_parse parse rd_natoms=44
rd_openff_good rd_parse
trace_181f38a7bb91f03268f985b9c235a4f405563ebcff8ed2eace58162ada4c0c7f
trace_8320c8e8a89100dd4f471ed9d311b60d04aa94704de0d3d81dc2d7f9ae3870b7
xcmp_isomorphic_ok xcmp_natom_ok xcmp_nbond_ok xcmp_smi_oe2_ok
xcmp_smi_oe_ok xcmp_smi_oe_rd_ok xcmp_smi_ok xcmp_smi_rd2_ok
xcmp_smi_rd_oe_ok xcmp_smi_rd_ok
DrugBank_2791	oe_natoms=19 oe_openff_good oe_parse parse rd_natoms=19
rd_openff_good rd_parse
trace_08be8cfd08795036334181c48bbd00b04e2002fa4b223cec3385abc50266a863
trace_cad8d98ca19220fb7dbb98fc19229d24d527c2407a0885572f092f5751165a54
xcmp_isomorphic_ok xcmp_natom_ok xcmp_nbond_ok xcmp_smi_oe2_ok
xcmp_smi_oe_ok xcmp_smi_oe_rd_ok xcmp_smi_ok xcmp_smi_rd2_ok
xcmp_smi_rd_oe_ok xcmp_smi_rd_ok
 ```

I used Z3 to minimize the result, which it did in under a second:

```
% python off_coverage.py minimize MiniDrugBank.feats -q
Solution: optimal Selected: 72/371
1	DrugBank_5354
1	DrugBank_2800
1	DrugBank_5418
1	DrugBank_104
1	DrugBank_5516
1	DrugBank_5523
1	DrugBank_2967
1	DrugBank_246
1	DrugBank_2991
1	DrugBank_3028
1	DrugBank_3046
1	DrugBank_3087
1	DrugBank_390
1	DrugBank_5804
1	DrugBank_3346
1	DrugBank_3358
1	DrugBank_3406
1	DrugBank_5900
1	DrugBank_5902
1	DrugBank_3479
1	DrugBank_3503
1	DrugBank_3547
1	DrugBank_3565
1	DrugBank_6032
1	DrugBank_914
1	DrugBank_977
1	DrugBank_6182
1	DrugBank_3817
1	DrugBank_3954
1	DrugBank_6355
1	DrugBank_4074
1	DrugBank_1448
1	DrugBank_1449
1	DrugBank_4138
1	DrugBank_4161
1	DrugBank_1538
1	DrugBank_6531
1	DrugBank_1564
1	DrugBank_6533
1	DrugBank_4215
1	DrugBank_4217
1	DrugBank_1598
1	DrugBank_1608
1	DrugBank_4249
1	DrugBank_1637
1	DrugBank_1661
1	DrugBank_6647
1	DrugBank_4323
1	DrugBank_1721
1	DrugBank_1722
1	DrugBank_1742
1	DrugBank_6722
1	DrugBank_4468
1	DrugBank_4515
1	DrugBank_6865
1	DrugBank_4580
1	DrugBank_1971
1	DrugBank_4662
1	DrugBank_7108
1	DrugBank_7124
1	DrugBank_2140
1	DrugBank_2148
1	DrugBank_2186
1	DrugBank_2237
1	DrugBank_4959
1	DrugBank_2429
1	DrugBank_5154
1	DrugBank_2563
1	DrugBank_2570
1	DrugBank_2584
1	DrugBank_2585
1	DrugBank_2684
```

In hand-wavy theory this means the test set can go from 371 structures
to only 72!

Or down to 27 if `func-pairs` gives valid coverage.

Of course, this is only the first pass, and it doesn't include any of
the topology considerations. Still, it seems a promising path.

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

