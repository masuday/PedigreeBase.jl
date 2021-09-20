module PedigreeBase

using SparseArrays
using OffsetArrays
using SparseMatrixDicts

export read_ped, check_ped, set_upg
export find_ped_order, permute_ped!, extract_ped, get_upg_range, subset_ped, mgs_ped
export get_inb, get_inbupg_code
export get_nrm, get_nrminv, get_mgsnrminv, get_qm, get_pm

include("io.jl")
include("permutation.jl")
include("inbreeding.jl")
include("relationship.jl")

end
