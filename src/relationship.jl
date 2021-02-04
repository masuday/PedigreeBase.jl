"""
    A = get_nrmsub(pedlist [, check=false, sorthere=true])
    A = get_nrmsub(T, pedlist [, check=false, sorthere=true])

Calculates the numerator (or additive) relationship matrix, so-called *A-matrix* using the recursive ("tabular") method from the pedigree list `pedlist`.
The pedigree must be sorted as the ancestors precede their progeny unless putting `sorthere=true`.
The resulting matrix is dense. You can specify the type of `A` as `T` (Float64 by default).
By default, this function checks if the pedigree is sorted chronologically or not; use `check=false` to turn it off.
With the option `sorthere=true`, the function reorders the pedigree internally, and returns the correct A-matrix with the original order.
However, it creates a temporary A-matrix with the reordered index, and it needs twice more memory.
"""
function get_nrm(Tv::DataType, pedlist::Matrix{Ti}; check=true, sorthere=false) where Ti<:Integer
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
            error("error in pedigree, or being not sorted")
         end
      end
   end
   # tabular method
   if sorthere
      # copy the pedigree list
      ped = normalize_pedigree(pedlist)
      perm,invp = find_ped_order(ped)
      permute_ped!(invp,ped)
      B = tabular_method_dense(Tv, ped)
      return B[invp,invp]
   else
      return tabular_method_dense(Tv, pedlist)
   end
end
function get_nrm(pedlist::Matrix{Ti}, Tv=Float64::DataType; check=false, sorthere=false) where Ti<:Integer
   return get_nrm(Tv, pedlist, check=check, sorthere=sorthere)
end

# lower triangular only
function tabular_method_dense(Tv::DataType, pedlist::Matrix{Ti}) where Ti<:Integer
   n::Ti = size(pedlist,2)
   sid::Ti = 0
   did::Ti = 0
   A = zeros(Tv,n,n)
   for j=1:n
      # IDs
      sid = pedlist[1,j]
      did = pedlist[2,j]
      # diagonal
      A[j,j] = 1 + aij(A,n,sid,did)/2
      # off diagonal (i>j i.e. sid<j & did<j)
      for i=j+1:n
         sid = pedlist[1,i]
         did = pedlist[2,i]
         A[i,j] = (aij(A,n,j,sid)+aij(A,n,j,did))/2
         A[j,i] = A[i,j]
      end
   end
   return A
end

# refers to A(i,j)
function aij(A::Matrix{Tv}, n::Ti, i::Ti, j::Ti) where {Tv<:AbstractFloat,Ti<:Integer}
   if i<1 || i>n
      val = zero(Tv)
   elseif j<1 || j>n
      val = zero(Tv)
   else
      val = A[i,j]
   end
   return val
end

"""
    Ainv = get_nrminv(pedlist, f; rank=0)
    Ainv = get_nrminv(pedlist; f, rank=0)

Computes the inverse of numerator (additive) relationship matrix using the "Henderson's" method from the pedigree list `pedlist` and inbreeding coefficients `f`.
Without `f`, this function assumes zero inbreeding for all individuals.
See Henderson (1976) and Quaas (1976) for details.
The output is SparseMatrixCSC with the same type of elements with `f` while SparseMatrixDict (i.e., a hash matrix) is used as a temporary matrix.
If you want to fix the size of A-inverse, specify `rank`.
The function creates the inverse with either a bigger size: `rank` or the maximum ID in the pedigree.
"""
function get_nrminv(pedlist::Matrix{Ti}, f::Vector{Tv}; rank=0) where {Ti<:Integer,Tv<:AbstractFloat}
   # there are UPGs if n<m
   n::Ti = size(pedlist,2)
   m::Ti = max(maximum(pedlist),n,rank)
   AInv = SparseMatrixDict{Tv,Int}(m,m)
   w::Array{Tv,1} = [1.0, -0.5, -0.5]
   triplet = zeros(Ti,3)
   for k=1:n
      sid = pedlist[1,k]
      did = pedlist[2,k]
      triplet .= [k, sid, did]
      if sid<1 || sid>n
         fs = 0.0
         ms = 1.0
      else
         fs = f[sid]
         ms = 0.0
      end
      if did<1 || did>n
         fd = 0.0
         md = 1.0
      else
         fd = f[did]
         md = 0.0
      end
      b = 4.0/((1.0-fs)*(1.0+ms)+(1.0-fd)*(1.0+md))
      for j=1:3
         for i=1:3
            aj = triplet[j]
            ai = triplet[i]
            #if (ai>=1 && ai<=n) && (aj>=1 && aj<=n)
            if ai>=1 && aj>=1
               val = b*w[i]*w[j]
               AInv[ai,aj] = AInv[ai,aj] + val
            end
         end
      end
   end
   return sparse(AInv)
end
function get_nrminv(pedlist::Matrix{Ti}; f::Vector{Tv}=Float64[], rank=0) where {Ti<:Integer,Tv<:AbstractFloat}
   if length(f)==0
      inb = zeros(Float64,maximum(pedlist))
   else
      inb = f
   end
   return get_nrminv(pedlist,inb, rank=rank)
end

# change UPG code to 0
function normalize_id_unchanged(pid,maxid)
   return pid
end
function normalize_id_removing_upg(pid,maxid)
   return ifelse(pid>maxid,0,pid)
end

"""
    Q = get_qm(pedlist)
    Q = get_qm(T, pedlist)

Computes the Q-matrix relating unknown parent groups (UPG) to individuals, as suggested by Quaas (1988).
The pedigree list `pedlist` should include a UPG code that is greater than the maximum ID of real animals.
The type of `Q` is `T` (Float64 by default).
"""
function get_qm(Tv::DataType, pedlist::Matrix{Ti}) where Ti<:Integer
   # real animals
   n::Ti = size(pedlist,2)
   # UPGs
   nupg::Ti = maximum(pedlist) - n
   # no UPGs
   if nupg < 1
      error("No groups in the pedigree list")
   end
   # start computations
   Q = zeros(Tv,n,nupg)
   Qp = zeros(Tv,nupg)
   for i=1:n
      Qp .= zeros(Tv,nupg)
      calculate_upg_q!(i,pedlist,n,1.0,Qp)
      Q[i,:] = Qp
   end
   return Q
end
function get_qm(pedlist::Matrix{Ti},Tv::DataType=Float64) where Ti<:Integer
   return get_qm(Tv,pedlist)
end

function calculate_upg_q!(aid::Ti, pedlist::Matrix{Ti}, maxanim::Ti, g::Tv, Qp::Vector{Tv}) where {Ti<:Integer,Tv<:AbstractFloat}
   if aid > maxanim
      group = aid - maxanim
      Qp[group] = Qp[group] + g
   elseif aid > 0
      calculate_upg_q!(pedlist[1,aid],pedlist,maxanim,g/2,Qp)
      calculate_upg_q!(pedlist[2,aid],pedlist,maxanim,g/2,Qp)
   end
end

"""
    P = get_pm(pedlist)
    P = get_pm(T, pedlist)

Computes a P-matrix that links parents with their progeny, as shown in Quaas (1988).
The type of `P` is `T` (Float64 by default).
"""
function get_pm(Tv::DataType, pedlist::Matrix{Ti}; verbose::Bool=true) where Ti<:Integer
   # real animals
   n = size(pedlist,2)
   # UPGs
   nupg = maximum(pedlist) - n
   # zero index
   if minimum(pedlist)<1
      if verbose
         @warn "pedigree with 0 index; not all animals assined to UPG"
      end
      P = zeros(Tv,n,n)
   else
      P = zeros(Tv,n,n+nupg)
   end
   # start computations
   for i=1:n
      sid = pedlist[1,i]
      did = pedlist[2,i]
      if sid>0
         P[i,sid] = 0.5
      end
      if did>0
         P[i,did] = 0.5
      end
   end
   return P
end
function get_pm(pedlist::Matrix{Ti}, Tv::DataType=Float64; verbose::Bool=true) where Ti<:Integer
   return get_pm(Tv,pedlist,verbose=verbose)
end
