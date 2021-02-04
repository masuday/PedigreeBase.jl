# return perm[newid]=oldid sorted by idtable
# perm[newid]  = origid
# invp[origid] = newid
function sort_by_dict(idtable::Dict{String,Int},n::Int)
   perm = zeros(Int,n)
   invp = zeros(Int,n)
   newid = 0
   for key in sort(collect(keys(idtable)))
      origid = idtable[key]
      if key!="0" && origid<=n
         newid = newid + 1
         perm[newid] = origid
         invp[origid] = newid
      end
   end
   return perm,invp
end

function get_qb(upgtable)
   nunk = upgtable.count
   nupg = length(unique(values(upgtable)))
   Qb = zeros(nunk,nupg)
   i = 0
   for key in sort(collect(keys(upgtable)))
      i = i + 1
      group = upgtable[key]
      Qb[i,group] = 1.0
   end
   return Qb
end
