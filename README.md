# PedigreeBase

## Quick start

This package provides a simple set of functions to deal with a pedigree structure, inbreeding coefficients, and relationship matrices in quantitative genetics.

## Example

```julia
# pedigree list using integer code

# define a pedigree list
pedlist = [0 0 1 1 3 1 5; 0 0 0 2 4 4 6]

# or read a file to load the pedigree list
pedlist = read_ped("pedigree.txt",integer=true)

# reordering pedigree so that parents precede their progeny
perm,invp = find_ped_order(pedlist)
permute_ped!(invp,pedlist)

# -------
# characters in a pedigree file
pedlist,idtable = read_ped("pedigree.characters.txt")

# reordering
perm,invp = find_ped_order(pedlist)
permute_ped!(invp,pedlist,idtable)
# -------

# inbreeding
f = get_inb(pedlist)

# A-matrix
A = get_nrm(pedlist)

# A-inverse
Ainv = get_nrminv(pedlist, f)
```
