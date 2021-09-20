using PedigreeBase
using LinearAlgebra
using Test

include("functions.jl")

@testset "Reading pedigree" begin
   # character file
   pedlist,idtable = read_ped("ped.henderson.char")
   @test pedlist == [0 0 1 1 3 1 5; 0 0 0 2 4 4 6]
   @test idtable == Dict("A"=>1,"B"=>2,"C"=>3,"D"=>4,"E"=>5,"F"=>6,"G"=>7,"0"=>0)
   try
      pedlist = read_ped("ped.henderson.char",integer=true)
      @test false
   catch
      # must be failed
      @test true
   end
   # integer file
   pedlist = read_ped("ped.henderson.int",integer=true)
   @test pedlist == [0 0 1 1 3 1 5; 0 0 0 2 4 4 6]
   # shuffled file
   pedlist,idtable = read_ped("ped.henderson.char.2",order=vec([2 4 5]))
   @test pedlist==[2 0 0 5 2 4 2; 3 0 0 1 0 7 1]
   @test idtable == Dict("A"=>2,"B"=>3,"C"=>5,"D"=>1,"E"=>4,"F"=>7,"G"=>6,"0"=>0)
   pedlist = read_ped("ped.henderson.int.2",order=vec([2 3 1]),integer=true)
   @test pedlist==[0 0 1 1 3 1 5; 0 0 0 2 4 4 6]
end

@testset "Wrong file format" begin
   # wrong file name
   @test_throws ErrorException read_ped("ThereIsNotSuchAFile.")
   # invalid order
   @test_throws BoundsError read_ped("ped.henderson.char.2",order=vec([2 8 5]))
   @test_throws ErrorException read_ped("ped.henderson.char.2",order=vec([2 4 0]))
   # missing columns
   @test_throws BoundsError read_ped("ped.henderson.char.error1")
   @test_throws BoundsError read_ped("ped.henderson.int.error1",integer=true)
   # non-integer data
   @test_throws ErrorException read_ped("ped.henderson.int.error2",integer=true)
   # duplication
   println("The following warning is expected.")
   @test_logs (:warn,"ID: 'A' appears in this file again. The last definision will be used.") read_ped("ped.henderson.char.error3")
   @test_logs (:warn,"ID: '1' appears in this file again. The last definision will be used.") read_ped("ped.henderson.int.error3",integer=true)
end

@testset "Logical error in pedigree" begin
   # okay
   pedlist,idtable = read_ped("ped.henderson.char")
   @test check_ped(pedlist)
   pedlist = read_ped("ped.henderson.int",integer=true)
   @test check_ped(pedlist)
   # mixed sire/dam
   pedlist,idtable = read_ped("ped.henderson.char.error4")
   #@test !check_ped(pedlist)
   @test_logs (:warn,"sire/dam found in the alternative sex: 1") check_ped(pedlist)
   pedlist = read_ped("ped.henderson.int.error4",integer=true)
   @test_logs (:warn,"sire/dam found in the alternative sex: 1") check_ped(pedlist)
   # loops
   pedlist,idtable = read_ped("ped.henderson.char.error5")
   @test_logs (:warn,"pedigree loops found: 3") check_ped(pedlist)
   pedlist = read_ped("ped.henderson.int.error5",integer=true)
   @test_logs (:warn,"pedigree loops found: 3") check_ped(pedlist)
   # progeny precede ancestor
   pedlist,idtable = read_ped("ped.henderson.char.3")
   @test_logs (:warn,"pedigree not chronologically sorted: 1") check_ped(pedlist, parentsfirst=true)
   pedlist = read_ped("ped.henderson.int.3",integer=true)
   @test_logs (:warn,"pedigree not chronologically sorted: 1") check_ped(pedlist, parentsfirst=true)
end

@testset "Permutation of pedigree" begin
   # character pedigree
   pedlist,idtable = read_ped("ped.henderson.char.3")
   perm,invp = find_ped_order(pedlist)
   permute_ped!(invp,pedlist,idtable)
   @test check_ped(pedlist, parentsfirst=true)
   # integer pedigree
   pedlist = read_ped("ped.henderson.int.3", integer=true)
   perm,invp = find_ped_order(pedlist)
   permute_ped!(invp,pedlist)
   @test check_ped(pedlist, parentsfirst=true)
   # make a subset
   pedlist,idtable = read_ped("ped.henderson.char.3")
   perm,invp = find_ped_order(pedlist,vec([7 3]))
   newpedlist,newidtable = extract_ped(invp,pedlist,idtable)
   @test check_ped(newpedlist, parentsfirst=true)
   # the same but different starting animals
   pedlist,idtable = read_ped("ped.henderson.char.3")
   perm,invp = find_ped_order(pedlist,vec([5 2 7]))
   newpedlist,newidtable = extract_ped(invp,pedlist,idtable)
   @test check_ped(newpedlist, parentsfirst=true)
end

@testset "Pedigree with UPG" begin
   pedlist,idtable = read_ped("ped.henderson.char.4")
   upgtable = Dict("x"=>1,"y"=>2)
   pedlist,idtable = set_upg(pedlist,idtable,upgtable)
   perm,invp = find_ped_order(pedlist)
   permute_ped!(invp,pedlist,idtable)
   @test check_ped(pedlist, parentsfirst=true)

   pedlist,idtable = read_ped("ped.henderson.char.4")
   upgtable = Dict("x"=>1,"y"=>2)
   pedlist,idtable = set_upg(pedlist,idtable,upgtable)
   perm,invp = find_ped_order(pedlist,vec([6]))
   newpedlist,newidtable = extract_ped(invp,pedlist,idtable)
   @test check_ped(newpedlist, parentsfirst=true)
end

@testset "Inbreeding coefficients" begin
   # reference
   pedlist1,idtable1 = read_ped("ped.henderson.char")
   f1 = get_inb(pedlist1, check=true)
   # shuffled data
   pedlist2,idtable2 = read_ped("ped.henderson.char.2",order=vec([2 4 5]))
   perm,invp = find_ped_order(pedlist2)
   permute_ped!(invp,pedlist2,idtable2)
   f2 = get_inb(pedlist2, check=true)
   @test f1 ≈ f2
   # reordered data
   pedlist3,idtable3 = read_ped("ped.henderson.char.3")
   perm,invp = find_ped_order(pedlist3)
   permute_ped!(invp,pedlist3,idtable3)
   f3 = get_inb(pedlist3, check=true)
   @test f1 ≈ f3

   # integer case
   pedlist=read_ped("ped.henderson.int",integer=true)
   s = pedlist[1,:]
   d = pedlist[2,:]
   f1a = get_inb(s,d)
   f1b = get_inb(pedlist)
   @test f1a ≈ f1b
end

@testset "Inbreeding coefficients with sorting" begin
   # reference
   pedlist1,idtable1 = read_ped("ped.henderson.char")
   f1 = get_inb(pedlist1, check=true)
   # try again Case 3 keeping the original order
   pedlist3,idtable3 = read_ped("ped.henderson.char.3")
   perm,invp = find_ped_order(pedlist3)
   permped3 = copy(pedlist3)
   permute_ped!(invp,permped3)
   f3 = get_inb(permped3, check=true)
   f3 .= f3[invp]
   nerr=0
   for k in sort(collect(keys(idtable1)))
      v=idtable1[k]
      if v>0
         #println((k,v,f1[v]))
         if f1[v]!=f3[idtable3[k]]
            nerr = nerr + 1
         end
      end
   end
   @test nerr==0

   # built-in sorting
   f3s = get_inb(pedlist3, check=true, sorthere=true)
   @test f3s ≈ f3
end

@testset "Numerator relationship matrix" begin
   # reference
   pedlist1,idtable1 = read_ped("ped.henderson.char")
   A1 = get_nrm(pedlist1)
   f1 = get_inb(pedlist1, check=true)
   n = size(A1,1)
   p1,ip1 = sort_by_dict(idtable1,n)
   A1 = A1[p1,p1]
   @test f1==diag(A1) .- 1.0
   # shuffled data
   pedlist2,idtable2 = read_ped("ped.henderson.char.2",order=vec([2 4 5]))
   perm,invp = find_ped_order(pedlist2)
   permute_ped!(invp,pedlist2,idtable2)
   A2 = get_nrm(pedlist2)
   p2,ip2 = sort_by_dict(idtable2,n)
   A2 = A2[p2,p2]
   @test A1==A2
   # shuffled with UPGs
   pedlist3,idtable3 = read_ped("ped.henderson.char.4")
   upgtable3 = Dict("x"=>1,"y"=>2)
   pedlist3,idtable3 = set_upg(pedlist3,idtable3,upgtable3)
   perm,invp = find_ped_order(pedlist3)
   permute_ped!(invp,pedlist3,idtable3)
   A3 = get_nrm(pedlist3)
   p3,ip3 = sort_by_dict(idtable3,n)
   A3 = A3[p3,p3]
   @test A1==A3
end

@testset "A-inverse" begin
   # reference
   pedlist1,idtable1 = read_ped("ped.henderson.char")
   A1 = get_nrm(pedlist1)
   f1 = get_inb(pedlist1, check=true)
   n = size(A1,1)
   p1,ip1 = sort_by_dict(idtable1,n)
   A1 = A1[p1,p1]
   Ainv1 = Matrix( get_nrminv(pedlist1, f1) )
   Ainv1 = Ainv1[p1,p1]
   @test Ainv1 ≈ inv(A1)
   # shuffled data
   pedlist2,idtable2 = read_ped("ped.henderson.char.2",order=vec([2 4 5]))
   perm,invp = find_ped_order(pedlist2)
   permute_ped!(invp,pedlist2,idtable2)
   A2 = get_nrm(pedlist2)
   f2 = get_inb(pedlist2, check=true)
   p2,ip2 = sort_by_dict(idtable2,n)
   A2 = A2[p2,p2]
   Ainv2 = Matrix( get_nrminv(pedlist2, f=f2) )
   Ainv2 = Ainv2[p2,p2]
   @test Ainv2 ≈ inv(A2)
   # shuffled with UPGs
   pedlist3,idtable3 = read_ped("ped.henderson.char.4")
   upgtable3 = Dict("x"=>1,"y"=>2)
   pedlist3,idtable3 = set_upg(pedlist3,idtable3,upgtable3)
   perm,invp = find_ped_order(pedlist3)
   permute_ped!(invp,pedlist3,idtable3)
   A3 = get_nrm(pedlist3)
   f3 = get_inb(pedlist3, check=true)
   p3,ip3 = sort_by_dict(idtable3,n)
   A3 = A3[p3,p3]
   Ainv3 = Matrix( get_nrminv(pedlist3, f=f3) )
   Ainv3 = Ainv3[p3,p3]
   @test Ainv3 ≈ inv(A3)
end

@testset "A-inverse with UPG" begin
   pedlist=read_ped("ped.quaas.int",integer=true,hasupg=true)
   f = get_inb(pedlist)
   Ai = Matrix(get_nrminv(pedlist,f=f))
   M = zeros(6,6)
   M[1,1]=4/3; M[1,2]=-2/3; M[1,5]=-1/2; M[1,6]=-1/6;
   M[2,2]=11/6; M[2,3]=1/2; M[2,4]=-1; M[2,6]=-2/3;
   M[3,3]=3/2; M[3,4]=-1; M[3,5]=-1/2; M[3,6]=-1/2;
   M[4,4]=2; M[5,5]=1/2; M[5,6]=1/2; M[6,6]=5/6;
   #M = M + M' - diagm(diag(M))
   M = M + M' - Matrix(Diagonal(diag(M)))
   @test M ≈ Ai
   upgs = get_upg_range(pedlist)
   @test upgs == 5:6
end

@testset "Q matrix" begin
   # regular process
   pedlist,idtable = read_ped("ped.quaas.char")
   upgtable = Dict("a"=>1,"b"=>2,"c"=>2,"d"=>1,"e"=>2)
   pedlist,idtable = set_upg(pedlist,idtable,upgtable)
   perm,invp = find_ped_order(pedlist)
   newpedlist,newidtable = extract_ped(invp,pedlist,idtable)
   P1 = get_pm(newpedlist)
   Q1 = get_qm(newpedlist)
   cp = [newidtable["a"]-4, newidtable["b"]-4]
   rp,rip = sort_by_dict(newidtable,size(pedlist,2))

   # sorted by the original ID
   pedlist,idtable = read_ped("ped.quaas.char")
   @test_logs (:warn,"pedigree with 0 index; not all animals assined to UPG") get_pm(pedlist)
   P0 = get_pm(pedlist, verbose=false)
   p,ip = sort_by_dict(idtable,size(pedlist,2))
   Pb = P0[p,p][1:4,5:9]
   P = P0[p,p][1:4,1:4]
   Qb = get_qb(upgtable)
   Q0 = inv(Matrix(1.0I,4,4)-P)*Pb*Qb
   @test inv(Matrix(1.0I,4,4)-P) ≈ Matrix(1.0I,4,4) + P + P*P
   @test Q0 ≈ Q1[rp,cp]
end

@testset "Subset pedigree" begin
   #           1  2  3  4  5  6  7  8  9 10 11  12 13 14 15
   pedlist = [ 0  0  0  0  0  0  2  1  2  7  7  11 11  9 11;
               0  0  0  0  0  0  5  4  3  6  4   8 10 13 10]
   idlist = [14]
   perm,invp,subpedlist = subset_ped(pedlist,idlist)
   f = get_inb(pedlist)
   subf = get_inb(subpedlist)
   @test all( map(x->f[perm[x]],findall(perm.>0)) .== subf )
   @test all( map(x->subf[invp[x]],findall(invp.>0)) .== subf )
end

@testset "MGS pedigree" begin
   # Henderson (1975)
   # for all pedigree including females
   #           1  2 X_3 Y_4 Z_5 3_6 4_7 5_8 6_9
   pedlist = [ 0  0  1   2   2   2   1   0   7
               0  0  0   0   0   3   4   5   0]
   mgslist = mgs_ped(pedlist)
   ref_mgs = [ 0  0  1   2   2   2   1   0   7
               0  0  0   0   0   1   2   2   0]
   @test all(mgslist .== ref_mgs)

   # for all pedigree including males only
   #           1  2 X_3 Y_4 Z_5 3_6 4_7 5_8 6_9
   pedlist = [ 0  0  1   2   2   2   1   0   7
               0  0  0   0   0   3   4   5   0]
   males =   [true, true, false, false, false, true, true, true, true]
   mgslist = mgs_ped(pedlist, males)
   ref_mgs = [ 0  0  0   0   0   2   1   0   7
               0  0  0   0   0   1   2   2   0]
   @test all(mgslist .== ref_mgs)

   Ainv = get_mgsnrminv(mgslist, males)
   idx = [1,2,6,7,8,9]
   ref_Ainv = [
      21.8182 5.4545  -5.4545  -10.9091    0 0
      5.4545  22.8182 -10.9091 -5.4545  -4.0 0
     -5.4545 -10.9091  21.8182  0        0   0
     -10.9091 -5.4545  0        26.8182  0  -10.0
     0        -4.0     0        0       16.0 0
     0         0       0        -10.0   0    20.0
   ]
   @test isapprox(ref_Ainv,Ainv[idx,idx]*15,atol=1e-3)
end
