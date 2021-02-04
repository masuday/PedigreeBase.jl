"""
    f = get_inb(pedlist [,check, sorthere, method])
    f = get_inb(Tv, pedlist [,check, sorthere, method])

Computes inbreeding coefficients for animals in a pedigree list `pedlist`.
The list is an integer matrix (2 x *n*, where *n* is the number of animals).
The results are in the vector `f` with the type `Tv` (Float64 by default).
Unknown parent groups (UPGs) will be replaced with a missing code (`0`).

The pedigree should be coded as parents precede their progeny.
All descendent mush have a larger code than their ancestors.
By default, the function checks if the pedigree is "permutaed", and stops if it is wrong.
You can omit the check by `check=false`.
If you want to permute the pedigree chronologically within this function, put `sorthere=true`.

The default method is Meuwissen and Luo (1992).
This is the only method to be supported in this version.

### Examples
```jldoctest
julia> using PedigreeTools;

julia> pedlist,idtable = read_ped("pedigree.txt");

# assuming animals are chronologically sorted.
julia> f = get_inb(pedlist);

# sort animals here.
julia> f = get_inb(pedlist, sorthere=true);
```
"""
function get_inb(Tv::DataType, pedlist::Matrix{Ti};
                 check=true, sorthere=false, method::String="MeuwissenAndLuo") where Ti<:Integer
   # make a copy of pedigree
   ped = normalize_pedigree(pedlist)

   # pedigree checks
   if check
      if sorthere
         # being sorted so no checks for `parentsfirst`
         if !check_ped(ped)
            error("error in pedigree")
         end
      else
         # assuming to be already sorted
         if !check_ped(ped,parentsfirst=true)
            error("error in pedigree, or not sorted")
         end
      end
   end
   # computations
   if method=="MeuwissenAndLuo"
      # This algorithm rewrites the pedigree list so
      # we use the copy of it to keep the original unchanged.
      if sorthere
         perm,invp = find_ped_order(ped)
         permute_ped!(invp,ped)
         f = kernel_meuwissen_and_luo!(Tv,ped)
         f .= f[invp]
         return f
      else
         return kernel_meuwissen_and_luo!(Tv,ped)
      end
      #return get_inb_ml(Tv,pedlist)
   else
      throw(ArgumentError("unsupported method for inbreeding computations"))
   end
end

# default = Float64
function get_inb(pedlist::Matrix{Ti}, Tv::DataType=Float64;
                 check=false, sorthere=false, method::String="MeuwissenAndLuo") where Ti<:Integer
   return get_inb(Tv, pedlist, check=check, sorthere=sorthere, method=method)
end

# normalize pedigree
# replace a negative/large ID with 0
function normalize_pedigree(pedlist::Matrix{Ti}) where {Ti<:Integer}
   n::Ti = size(pedlist,2)
   ped = zeros(Ti,2,n)
   sid::Ti = 0
   did::Ti = 0
   # normalize pedigree; replace negative/large id with 0
   for i=1:n
      sid = pedlist[1,i]
      did = pedlist[2,i]
      if sid<1 || sid>n
         sid = 0
      end
      if did<1 || did>n
         did = 0
      end
      ped[1,i] = sid
      ped[2,i] = did
   end
   return ped
end

# kernel of Meuwissen and Luo (1992)
# the argument "ped" being rewritten
function kernel_meuwissen_and_luo!(Tv::DataType, ped::Matrix{Ti}) where Ti<:Integer
   n::Ti = size(ped,2)
   #f = OffsetArray{Tv}(0:n)
   f = OffsetArray{Tv}(undef,0:n)
   i::Ti = 0
   j::Ti = 0
   k::Ti = 0
   s0::Ti = 0
   d0::Ti = 0
   ks::Ti = 0
   kd::Ti = 0
   fi::Tv = 0.0
   r::Tv = 0.0
   point = zeros(Ti,n)
   T = zeros(Tv,n)
   B = zeros(Tv,n)

   # start
   f[0] = -1
   for i=1:n
      s0 = ped[1,i]
      d0 = ped[2,i]
      ped[1,i] = max(s0,d0)
      ped[2,i] = min(s0,d0)
      B[i] = 0.5 - 0.25*(f[s0]+f[d0])

      if s0==0 || d0==0
         f[i] = 0.0
      else
         fi = -1.0
         T[i] = 1.0
         j = i
         while(j!=0)
            k = j
            r = 0.5 * T[k]
            ks = ped[1,k]
            kd = ped[2,k]
            if ks != 0
               while(point[k]>ks)
                  k = point[k]
               end
               T[ks] = T[ks] + r
               if ks != point[k]
                  point[ks] = point[k]
                  point[k] = ks
               end
               if kd != 0
                  while(point[k]>kd)
                     k = point[k]
                  end
                  T[kd] = T[kd] + r
                  if kd != point[k]
                     point[kd] = point[k]
                     point[k] = kd
                  end
               end
            end
            fi = fi + T[j]*T[j]*B[j]
            T[j] = 0.0
            k = j
            j = point[j]
            point[k] = 0
         end
         f[i] = fi
      end
   end
   return f[1:n]
end
