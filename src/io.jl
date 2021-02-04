"""
    pedlist,idtable = read_ped(filename [, integer, order, hasupg])
    pedlist,idtable = read_ped(T, filename [, integer, order, hasupg])

Reads a pedigree file and store the pedigree list to an integer matrix `pedlist`
(2 x *n*; where *n* is the number of animals found in the pedigree). The argument
`T` defines the data type of `pedlist` and the default is Int. This function
supports two types of ID: character and integer.

By default, this function reads the original as characters and converts it to an
integer. The table is in a dictionary `idtable`, and it comes as the 2nd return-value.
If the ID is an integer, you can put the option `integer=true`, where the function
simply reads it as an integer without any conversion. No table is created so you will
have a single return value, `pedlist`.

With integer coding, the function treats all sires and dams not in the animal field
as real animals and creates new entries for them. It is not appropriate when your
pedigree file is "renumbered" by RENUMF90 with unknown parent groups (UPGs) where
UPG code is larger than the ID for real animals. The option `hasupg=true` does
nothing for UPG in sire/dam fields so that such UPG-code stays in `pedlist` as is.
The pedigree file has a space-separated text format with at least three fields: animal,
sire, and dam ID. By default, the first three fields are assigned to these IDs.

The optional argument `order` changes the position of fields. The default is
`order=[1,2,3]` and the 1st field for animal, the 2nd for sire, and the 3rd for dam.

### Examples
```jldoctest
julia> using PedigreeTools;

# for character data
julia> pedlist,idtable = read_ped("pedigree.txt");

# for character data with Int32 storage
julia> pedlist,idtable = read_ped(Int32,"pedigree.txt");

# for integer coding
julia> pedlist = read_ped("pedigree.txt", integer=true);
```
"""
function read_ped(Ti::DataType, filename::String; integer=false, order::Vector{To}=[1,2,3], hasupg=false) where {To<:Integer}
   if count(order.<1)>0
      error("`order` has 0 or negative value(s)")
   end
   if !(Ti<:Integer)
      error("non integer datatype specified")
   end
   if integer
      return read_ped_integer(Ti, filename, order=order; hasupg=hasupg)
   else
      return read_ped_general(Ti, filename, order=order)
   end
end

function read_ped(filename::String, Ti::DataType=Int; integer=false, order::Vector{To}=[1,2,3], hasupg=false) where {To<:Integer}
   return read_ped(Ti, filename, integer=integer, order=order, hasupg=hasupg)
end

# a general pedigree file; including characters
# change the name from ReadPedigreeGeneral to read_ped_general
function read_ped_general(Ti::DataType, filename::String; order::Vector{To}=[1,2,3]) where To<:Integer
   idtable = Dict{String,Ti}()
   idanim = Dict{String,Ti}()
   #pedtuple = Array{Tuple{Ti,Ti,Ti},1}()
   pedtuple = Array{NTuple{3,Ti},1}()
   line = ""
   pedline = zeros(Ti,3)
   lastid = zero(Ti)

   # checks
   if !isfile(filename)
      error("file not found")
   end

   # read the file
   idtable["0"] = 0
   open(filename,"r") do ios
      while(!eof(ios))
         line = strip(readline(ios))
         tokens = split(line)
         # animal
         id = tokens[order[1]]
         if haskey(idanim,id)
            @warn "ID: '"*id*"' appears in this file again. The last definision will be used."
            idanim[id] = idanim[id] + 1
         else
            idanim[id] = 1
         end
         # read the tokens
         for i=1:3
            id = tokens[order[i]]
            # no existing key
            if !haskey(idtable,id)
               lastid = lastid + 1
               idtable[id] = lastid
            end
            # set integer code
            pedline[i] = idtable[id]
         end
         # make a tuple
         push!(pedtuple,tuple(pedline...))
      end
   end

   # convert tuple to list
   pedlist = convert_pedigree_tuple_to_matrix(pedtuple,lastid)

   return pedlist,idtable
end

# an integer pedigree file; including characters
# change the name from ReadPedigreeInteger to read_ped_integer
function read_ped_integer(Ti::DataType, filename::String; order::Vector{To}=[1,2,3], hasupg=false) where To<:Integer
   #pedtuple = Array{Tuple{Ti,Ti,Ti},1}()
   pedtuple = Array{NTuple{3,Ti},1}()
   idanim = Dict{Ti,Ti}()
   lastid = zero(Ti)
   i = zero(Ti)
   k = zero(Ti)
   pedline = zeros(Ti,3)

   # checks
   if !isfile(filename)
      error("file not found")
   end

   # switch how to the maximum id
   find_maxid1(list::Vector{Ti},lastid) where Ti<:Integer = max(lastid,list[1])
   find_maxid2(list::Vector{Ti},lastid) where Ti<:Integer = max(lastid,maximum(list))
   if hasupg
      # maximum id only from the "animal" column
      find_maxid = find_maxid1
   else
      # maximum id only from all the columns
      find_maxid = find_maxid2
   end

   # read the file
   lastid = 0
   open(filename,"r") do ios
      while(!eof(ios))
         line = strip(readline(ios))
         tokens = split(line)
         # animal
         id = parse(Int,tokens[order[1]])
         if haskey(idanim,id)
            @warn "ID: '"*string(id)*"' appears in this file again. The last definision will be used."
            idanim[id] = idanim[id] + 1
         else
            idanim[id] = 1
         end
         # read the tokens
         for j=1:3
            try
               pedline[j] = parse(Int,tokens[order[j]])
            catch
               error("conversion error (non-integer value found): " * tokens[order[j]])
            end
         end
         # last id
         lastid = find_maxid(pedline,lastid)
         # make a tuple
         push!(pedtuple,tuple(pedline...))
      end
   end

   # convert tuple to list
   pedlist = convert_pedigree_tuple_to_matrix(pedtuple,lastid)

   return pedlist
end

# pedigree tuple to pedigree list
function convert_pedigree_tuple_to_matrix(pedtuple::Array{NTuple{3,Ti},1},lastid::Ti) where {Ti<:Integer}
   pedlist = zeros(Ti,2,lastid)
   for p in pedtuple
      if p[1]<1 || p[1]>lastid
         error("Invalid ID found")
      end
      pedlist[1:2,collect(p[1])] = collect(p[2:3])
   end
   return pedlist
end

# check pedigree
# change the name from CheckPedigree to check_ped
"""
    bool = check_ped(pedlist [, parentsfirst, warning])

Checks if the pedigree list has no apparent error. The pedigree list `pedlist` is
an integer matrix generated with `ReadPedigree` from a file. This function checks
the following errors.

* Sires (or dams) appear as dams (or sires).
* The list has a loop i.e., an animal's ancestor is also progeny of the animal.
* Parents precede their progeny i.e., IDs are chronologically sorted (optional; checked with `parentsfirst=true`).

It returns `true` with no error but `false` with a warning message where the error is (omit with `warining=false`).

### Examples

```jldoctest
julia> using PedigreeTools;

julia> pedlist,idtable = read_ped("pedigree.txt");

julia> check_ped(pedlist)
true

# checks if parents precede their progeny
julia> check_ped(pedlist, parentsfirst=true)
true
```
"""
function check_ped(pedlist::Matrix{Ti}; parentsfirst::Bool=false, warning::Bool=true) where Ti<:Integer
   n::Ti = size(pedlist,2)
   mark = zeros(Bool,n)

   # check sire/dam inconsistency
   mark = zeros(Bool,n)
   for i=1:n
      k = pedlist[1,i]
      if k>=1 && k<=n
         mark[k] = true
      end
   end
   for i=1:n
      k = pedlist[2,i]
      if k>=1 && k<=n
         if mark[k]
            if warning
               @warn "sire/dam found in the alternative sex: " * string(k)
            end
            return false
         end
      end
   end

   # check pedigree loop
   for i=1:n
      loop::Bool = hasloop(pedlist,i,i,0,false)
      if loop
         if warning
            @warn "pedigree loops found: " * string(i)
         end
         return false
      end
   end

   # check if sorted chronologically
   if parentsfirst
      for i=1:n
         unsorted::Bool = csorted(pedlist,i,i,false)
         if unsorted
            if warning
               @warn "pedigree not chronologically sorted: " * string(i)
            end
            return false
         end
      end
   end

   return true
end

# recursive tracing back of pedigree
function hasloop(pedlist::Matrix{Ti},ini::Ti,i::Ti,gen::Int,loop::Bool) where Ti<:Integer
   # immediate exit if it has a loop
   if loop==true
      return true
   end
   # existing animal
   if i>=1 && i<=size(pedlist,2)
      if gen>0 && i==ini
         # the same animal found as ancestor = pedigree loop
         return true
      else
         # tracing back
         return hasloop(pedlist,ini,pedlist[1,i],gen+1,loop) |
                hasloop(pedlist,ini,pedlist[2,i],gen+1,loop)
      end
   else
      # unknown base animal; do nothing
      return loop
   end
end

# check if chronologically sorted
function csorted(pedlist::Matrix{Ti},ini::Ti,i::Ti,unsorted::Bool) where Ti<:Integer
   # immediate exit if unsorted
   if unsorted==true
      return true
   end
   # existing animal
   if i>=1 && i<=size(pedlist,2)
      if i > ini
         # ancestor with larger ID than the staring animal = unsorted
         return true
      else
         # tracing back
         return csorted(pedlist,ini,pedlist[1,i],unsorted) |
                csorted(pedlist,ini,pedlist[2,i],unsorted)
      end
   else
      # unknown base animal; do nothing
      return unsorted
   end
end

"""
    newpedlist,newidtable = set_upg(pedlist, idtable, upgtable)

Assigns unknown parent groups (UPG) to unknown parents following the rule
in `upgtable` for a pedigree list `pedlist` and the original ID table `idtable`
See the following example;
5 UPG and the animals with IDs "A", "B", and "C" are grouped into the 1st UPG and "D" and "E" are grouped into the 2nd UPG.

```
upgtable = Dict("A"=>1, "B"=>1, "C"=>1, "D"=>2, "E"=>2)
```

The "animals" defined as UPG in `upgtable` are eliminated from the pedigree list.
If there are parent animals, they are replaced with integer UPG-code larger than the number of real animals. The table `idtable` is also updated accordingly.
"""
function set_upg(pedlist::Matrix{Ti},idtable::Dict{String,Ti},upgtable::Dict{String,Tj}; sortkeys=false) where {Ti<:Integer,Tj}
   #newpedtuple = Array{Tuple{Int,Int,Int},1}()
   newpedtuple = Array{NTuple{3,Ti},1}()
   invidtable = Dict{Ti,String}()
   newidtable = Dict{String,Ti}()
   upgdef = Dict{Tj,Ti}()
   pedline = zeros(Ti,3)
   lastcode = zero(Ti)
   # sort
   if sortkeys
      sortfunc(x) = sort(collect(x))
   else
      sortfunc = x->x
   end
   # the number of UPGs
   nupg=0
   for upgid in sortfunc(keys(upgtable))
      if upgid!="0" && haskey(idtable,upgid)
         nupg = nupg + 1
      end
   end
   # the number of real animals
   n = size(pedlist,2)
   m = n - nupg
   if m<1
      error("No real animals")
   end
   lastupg = m
   # make a new table with upg
   invidtable[0] = "0"
   newidtable["0"] = 0
   for id in sortfunc(keys(idtable))
      if id!="0"
         code = idtable[id]
         invidtable[code] = id
         if haskey(upgtable,id)
            # this id is upg: assign a large number (>m)
            if !haskey(newidtable,id)
               # grup ID
               group = upgtable[id]
               if haskey(upgdef,group)
                  # exsiting definition
                  newidtable[id] = upgdef[group]
               else
                  # new definition
                  lastupg = lastupg + 1
                  newidtable[id] = lastupg
                  upgdef[group] = lastupg
               end
            end
         else
            # this id is real
            # no existing key
            if !haskey(newidtable,id)
               lastcode = lastcode + 1
               newidtable[id] = lastcode
            end
         end
      end
   end
   # make a new pedigree list
   for id in sortfunc(keys(idtable))
      if id!="0"
         # animal code with id
         code = idtable[id]
         # store only if it is a real animal
         if !haskey(upgtable,id)
            # sire id and dam id
            sid = invidtable[pedlist[1,code]]
            did = invidtable[pedlist[2,code]]
            # new code
            pedline[1] = newidtable[id]
            pedline[2] = newidtable[sid]
            pedline[3] = newidtable[did]
         end
         # make a tuple
         push!(newpedtuple,tuple(pedline...))
      end
   end

   # convert tuple to list
   newpedlist = convert_pedigree_tuple_to_matrix(newpedtuple,lastcode)

   return newpedlist,newidtable
end
