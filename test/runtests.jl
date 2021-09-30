module TestingOIFITS

using Test, Printf
using OIFITS, FITSIO
using OIFITS: extname, get_format, is_same, fix_name

const dir = @__DIR__

const files = ("contest-2004-obj1.oifits" ,"contest-2004-obj2.oifits",
               "contest-2008-binary.oifits", "contest-2008-obj1-H.oifits",
               "contest-2008-obj1-J.oifits", "contest-2008-obj1-K.oifits",
               "contest-2008-obj2-H.oifits", "contest-2008-obj2-J.oifits",
               "contest-2008-obj2-K.oifits", "contest-2008-obj1-J.oifits")

const counter = Ref(0)

const quiet = true

const tempfile = let tup = mktemp()
    close(tup[2])
    tup[1]
end

println("tempfile = \"$tempfile\"")

function tryload(dir, file)
    global counter
    try
        data = OIDataSet(joinpath(dir, file))
        counter[] += 1
        quiet || @info "file \"$file\" successfully loaded"
        return true
    catch
        @warn "failed to load \"$file\""
        return false
    end
end

@testset "Strings" begin
    @test  is_same("", "")
    @test  is_same("", " ")
    @test  is_same("\t", "")
    @test  is_same(" ZCMa", " zcma")
    @test  is_same(" ZCMa", " zcma ")
    @test !is_same("ZCMa", " zcma")
    @test  is_same("SAO-206462", "sao-206462")
    @test  is_same("SAO-206462", "sao-206462\t")
    @test  is_same("SAO 206462", "sao 206462")
    @test  is_same("SAO 206462", "sao 206462\t")
    @test !is_same("SAO 206462", "sao\t206462")

    @test fix_name("") === ""
    @test fix_name(" ") === ""
    @test fix_name("a\t ") === "A"
    @test fix_name("\tbcd\t ") === "\tBCD"
    @test fix_name("\téëï\t ") === "\tÉËÏ"
end

get_field_type(T::Type{<:OIDataBlock}, sym::Symbol) = fieldtype(T, sym)
get_field_type(::Type{<:OI_TARGET}, sym::Symbol) =
    fieldtype((sym === :revn ? OI_TARGET : OITargetEntry), sym)

@testset "OI type definitions and formats" begin
    for T in (OI_ARRAY,
              OI_CORR,
              OI_FLUX,
              OI_INSPOL,
              OI_T3,
              OI_TARGET,
              OI_VIS,
              OI_VIS2,
              OI_WAVELENGTH)

        @testset "$(extname(T))" begin
            for rev in 3:-1:1
                spec = get_format(T, rev; throw_errors=false)
                if spec === nothing
                    continue
                end
                for s in spec
                    @test s.type ∈ (:A, :C, :D, :E, :I, :J, :L)
                    S = get_field_type(T, s.symb)
                    rank = ndims(s)
                    if rank == 0
                        # Keyword.
                        if s.type === :A
                            @test S === String
                        elseif s.type === :C
                            @test S <: Complex{<:AbstractFloat}
                        elseif s.type === :E || s.type === :D
                            @test S <: AbstractFloat
                        elseif s.type === :I || s.type === :J
                            @test S <: Integer
                        elseif s.type === :L
                            @test S <: Bool
                        end
                    else
                        # Column.
                        if s.type === :A
                            @test rank == 1
                            @test S === (T === OI_TARGET ? String : Vector{String})
                        else
                            @test rank == (T === OI_TARGET ? 1 : ndims(S))
                            if s.type === :C
                                @test eltype(S) <: Complex{<:AbstractFloat}
                            elseif s.type ∈ (:D, :E)
                                @test eltype(S) <: AbstractFloat
                            elseif s.type ∈ (:I, :J)
                                @test eltype(S) <: Int
                            elseif s.type === :L
                                @test eltype(S) <: Bool
                            end
                        end
                    end
                end
            end
        end
    end
end

@testset "Read/write OI-FITS files" begin
    # Each file infdividually.
    for file in files
        A = OIDataSet(joinpath(dir, file))
        @test isa(A, OIDataSet)
        isfile(tempfile) && rm(tempfile; force=true)
        write(tempfile, A)
        B = read(OIDataSet, tempfile)
        @test isa(B, OIDataSet)
        @test_throws ErrorException write(tempfile, B)
    end
    A = OIDataSet(map(x -> joinpath(dir, x), files)...)
    @test isa(A, OIDataSet)
    @test A.instr[1].name ===  A.instr[1].insname
    @test A.array[1].name ===  A.array[1].arrname
    name = A.instr[end].name
    @test isa(A.instr[name], OI_WAVELENGTH)
    @test A.instr[name] === A.instr[end]
    @test A.instr[uppercase(name*" ")] === A.instr[end]
    @test A.instr[lowercase(name*" ")] === A.instr[end]
    @test_throws KeyError A.instr[" "*name]
    @test get(A.instr, " "*name, undef) === undef
    name = A.array[end].name
    @test isa(A.array[name], OI_ARRAY)
    @test A.array[name] === A.array[end]
    @test A.array[uppercase(name*" ")] === A.array[end]
    @test A.array[lowercase(name*" ")] === A.array[end]
    @test_throws KeyError A.array[" "*name]
    @test get(A.array, " "*name, undef) === undef
end

@testset "OI_TARGET methods" begin
    data = OIDataSet(joinpath(dir,files[1]))
    A = data.target
    @test isa(A, OI_TARGET)
    @test eltype(A) === OITargetEntry
    @test length(A) == length(A.list)
    @test ndims(A) == 1
    @test IndexStyle(A) === IndexLinear()
    @test firstindex(A) == 1
    @test lastindex(A) == length(A)
    @test size(A) == (length(A),)
    @test size(A,1) == length(A)
    @test size(A,2) == 1
    @test axes(A) == (1:length(A),)
    @test axes(A,1) == 1:length(A)
    @test eachindex(A) === Base.OneTo(length(A))
    @test keys(Int, A) === Base.OneTo(length(A))
    for (key, val) in zip(keys(String, A), values(A))
        @test A[key] === val
    end
    i = 0
    for tgt in A
        i += 1
        @test isa(tgt, OITargetEntry)
        @test A[i] === tgt
        @test haskey(A, i) == true
        for key in (tgt.target,
                    uppercase(tgt.target),
                    lowercase(tgt.target),
                    tgt.target*"  \t")
            @test haskey(A, key)
            @test !haskey(A, key*"+")
            @test A[key] === tgt
            @test get(A, key, nothing) === tgt
            @test get(A, key*"+", nothing) === nothing
        end
    end
end

isfile(tempfile) && rm(tempfile; force=true)

end # module

nothing
