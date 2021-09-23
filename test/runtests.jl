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
        data = OIData(joinpath(dir, file))
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
get_field_type(::Type{<:OITarget}, sym::Symbol) =
    fieldtype((sym === :revn ? OITarget : OITargetEntry), sym)

@testset "OI type definitions and formats" begin
    for T in (OIArray,
              OICorr,
              OIFlux,
              OIInsPol,
              OIT3,
              OITarget,
              OIVis,
              OIVis2,
              OIWavelength)

        @testset "$(extname(T))" begin
            for rev in 3:-1:1
                spec = get_format(T, rev; throwerrors=false)
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
                            @test S === (T === OITarget ? String : Vector{String})
                        else
                            @test rank == (T === OITarget ? 1 : ndims(S))
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
    for file in files
        A = OIData(joinpath(dir, file))
        @test isa(A, OIData)
        isfile(tempfile) && rm(tempfile; force=true)
        write(tempfile, A)
        B = read(OIData, tempfile)
        @test isa(B, OIData)
        @test_throws ErrorException write(tempfile, B)
    end
end

@testset "OI_TARGET methods" begin
    A = OIData(joinpath(dir,files[1]))
    @test isa(A.target, OITarget)
    @test length(A.target) == length(A.target.rows)
    i = 0
    for tgt in A.target
        i += 1
        @test isa(tgt, OITargetEntry)
        @test A.target[i] === tgt
        for key in (tgt.target,
                    uppercase(tgt.target),
                    lowercase(tgt.target),
                    tgt.target*"  \t")
            @test A.target[key] === tgt
        end
    end
end

isfile(tempfile) && rm(tempfile; force=true)

end # module

nothing
