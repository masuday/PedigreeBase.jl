"""
    perm,invp = find_ped_order(pedlist [,start])

Traverses the pedigree and finds a permutation so that the ancestors precede their progeny.
The resulting vector `perm` associates the current ID with the new ID i.e., `perm[newid]=oldid`.
The vector `invp` has the inverse permutation i.e. `invp[oldid]=newid`.
By default, the function takes each animal as the starting animal.
When supplying an optional vector `start` with starting IDs, the tracing starts only with these progeny.
In this case, `perm` and `invp` will include `0` for non-tracked animals.
"""
function find_ped_order(pedlist::Matrix{Ti},start::Vector{Ti}=zeros(Ti,0)) where Ti<:Integer
   n::Ti = size(pedlist,2)
   perm = zeros(Ti,n)
   invp = zeros(Ti,n)
   newid = zero(Ti)
   if length(start)>0
      for k=1:length(start)
         i = start[k]
         if i<1 || i>n
            error("starting ID out of range")
         end
         newid = find_ped_order_within_family!(pedlist,i,newid,0,invp,perm)
      end
   else
      for i=1:n
         newid = find_ped_order_within_family!(pedlist,i,newid,0,invp,perm)
      end
   end
   return perm,invp
end

# assign new (permutaed) id using recursion within family
# perm[newid]  = origid
# invp[origid] = newid
function find_ped_order_within_family!(pedlist::Matrix{Ti},aid::Ti,lastaid::Ti,gen::Ti,invp::Vector{Ti},perm::Vector{Ti}) where Ti<:Integer
   sid = pedlist[1,aid]
   did = pedlist[2,aid]
   n::Ti = size(pedlist,2)
   newaid::Ti = lastaid

   # already registered in invp
   if invp[aid] > 0
      return newaid
   end
   # too many generations
   if gen>64
      error("too many generations back")
   end
   # not registered in invp
   if (sid<1||sid>n) && (did<1||did>n)
      # both unknown: increment new id

   elseif sid<1||sid>n
      # sire unknown: track the dam first; then register this animal
      newaid = find_ped_order_within_family!(pedlist,did,newaid,gen+1,invp,perm)
   elseif did<1||did>n
      # dam unknown: track the sire first; then register this animal
      newaid = find_ped_order_within_family!(pedlist,sid,newaid,gen+1,invp,perm)
   else
      # both known: just track both parents
      newaid = find_ped_order_within_family!(pedlist,sid,newaid,gen+1,invp,perm)
      newaid = find_ped_order_within_family!(pedlist,did,newaid,gen+1,invp,perm)
   end
   # not registered yet
   if invp[aid] < 1
      newaid = newaid + 1
      invp[aid] = newaid
      perm[newaid] = aid
   end
   return newaid
end

# a function return the renumbered ID for a parent
function get_new_parent_code(oldpid::Ti,n::Ti,pedlist::Matrix{Ti},invp::Vector{Ti}) where Ti<:Integer
   if oldpid<1
      # missing parent
      return 0
   elseif oldpid>n
      # unknown parent group
      return oldpid
   else
      # renumbered parent ID
      return invp[oldpid]
   end
end

# apply the permutation to pedigree
# perm[newid]  = origid
# invp[origid] = newid
"""
    permute_ped!(invp, pedlist [,idtable])

This function "renumbers" the animals by a permutation vector.
It updates the code in `pedlist` according to the permutation list in `invp` as `invp[oldid]=newid`, which is generated with some functions like `find_ped_order`.
The function also replaces IDs in `idtable` (optional).
A code for unknown-parent group, which is usually greater than the maximum ID of real animals, remains in the new pedigree.
"""
function permute_ped!(invp::Vector{Ti},pedlist::Matrix{Ti},idtable::Dict{String,Ti}=Dict{String,Ti}()) where Ti<:Integer
   # incomplete permutation vector
   n::Ti = length(invp)
   m::Ti = count((invp.>=1) .& (invp.<=n))
   if length(invp)!=size(pedlist,2)
      error("dimension mismatch between pedlist and invp")
   end
   if m!=n
      error("incomplete permutation vector(s); use ExtractPedigree() to extract the subset pedigree")
   end

   # swap the entry
   tmplist = zeros(Ti,2,n)
   for oldaid=1:n
      # animal
      newaid = invp[oldaid]
      # sire
      oldsid = pedlist[1,oldaid]
      newsid = get_new_parent_code(oldsid,n,pedlist,invp)
      # dam
      olddid = pedlist[2,oldaid]
      newdid = get_new_parent_code(olddid,n,pedlist,invp)
      # update
      tmplist[1:2,newaid] = [newsid; newdid]
   end
   pedlist[1:2,1:n] = tmplist[1:2,1:n]

   # change the ID table (optional)
   for id in keys(idtable)
      oldaid = idtable[id]
      if oldaid < 1
         newaid = 0
      elseif oldaid > n
         newaid = oldaid
      else
         newaid = invp[oldaid]
      end
      idtable[id] = newaid
   end
end

# extract the pedigree based on the permutation vector
# perm[newid]  = origid
# invp[origid] = newid
"""
    newlist,newtable = extract_ped(invp, pedlist ,idtable)
    newlist          = extract_ped(invp, pedlist)

Extract the code in `pedlist` according to the permutation list in `invp` as `invp[oldid]=newid`, which is generated with some functions like `find_ped_order`.
With optional `idtable`, it will also return with the new table.
The function re-assigns the code to animals in `newlist`.
The behavior is similar to `permute_ped`, but this function supports reducing the amount of pedigree.
"""
function extract_ped(invp::Vector{Ti},pedlist::Matrix{Ti},idtable::Dict{String,Ti}=Dict{String,Ti}()) where Ti<:Integer
   # animals in the full pedigree
   n::Ti = length(invp)
   # animals in the subset
   m::Ti = count((invp.>=1) .& (invp.<=n))
   if length(invp)!=size(pedlist,2)
      error("dimension mismatch between pedlist and invp")
   end

   # create a new subset of pedigree
   newlist = zeros(Ti,2,m)
   for oldaid=1:n
      # animal
      newaid = invp[oldaid]
      if newaid>=1 && newaid<=m
         # sire
         oldsid = pedlist[1,oldaid]
         newsid = get_new_parent_code(oldsid,n,pedlist,invp)
         # dam
         olddid = pedlist[2,oldaid]
         newdid = get_new_parent_code(olddid,n,pedlist,invp)
         # adjustment needed for unknown parent groups
         # convert oldupg>n to newupg>m by newid-n+m
         newsid = newsid>n ? newsid-n+m : newsid
         newdid = newdid>n ? newdid-n+m : newdid
         # update
         newlist[1:2,newaid] = [newsid, newdid]
      end
   end

   # change the ID table
   if length(idtable)>0
      newtable::Dict{String,Ti} = Dict{String,Ti}()
      for id in keys(idtable)
         # conversion with unknown parent groups
         oldaid = idtable[id]
         if oldaid < 1
            newaid = 0
         elseif oldaid > n
            newaid = oldaid-n+m
         else
            newaid = invp[oldaid]
         end
         newtable[id] = newaid
      end
      return newlist,newtable
   else
      return newlist
   end
end

"""
    invs = invsubidx(subidx [,n])

Returns the "inverse permutation" of sun-indices (`subidx[newid]=oldid`) resulting in `invs[oldid]=newid`.
The `subid` is similar to the permutation vector `perm`, but `subid` doesn't have to include all the index; `length(subid)` can be smaller than the maximum ID.
Additional argument `n` specifies the maximum ID in the pedigree.
"""
function invsubidx(subidx::Vector{T},n::Integer) where{T}
   invs = zeros(T,n)
   for i=1:length(subidx)
      if subidx[i]>0
         invs[subidx[i]] = i
      end
   end
   return invs
end

"""
    r = get_upg_range(pedlist)

    Finds the range of IDs for unknown parent groups.
"""
function get_upg_range(pedlist)
   return (size(pedlist,2)+1):maximum(pedlist)
end
