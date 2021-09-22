module TestingOIFITS

using Test, Printf
using OIFITS, FITSIO
using OIFITS: extname, get_format, is_same, fix_name

dir = @__DIR__

files = ("contest-2004-obj1.oifits" ,"contest-2004-obj2.oifits",
         "contest-2008-binary.oifits", "contest-2008-obj1-H.oifits",
         "contest-2008-obj1-J.oifits", "contest-2008-obj1-K.oifits",
         "contest-2008-obj2-H.oifits", "contest-2008-obj2-J.oifits",
         "contest-2008-obj2-K.oifits", "contest-2008-obj1-J.oifits")

counter = 0

quiet = true

function tryload(dir, file)
    global counter
    try
        data = OIData(joinpath(dir, file))
        counter += 1
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

@testset "OI type definitions and formats" begin
    F = Float64
    for T in (OIArray{F},
              OICorr{F},
              OIFlux{F},
              OIInsPol{F},
              OIT3{F},
              OITarget{F},
              OIVis{F},
              OIVis2{F},
              OIWavelength{F})

        @testset "$(extname(T))" begin
            for rev in 3:-1:1
                spec = get_format(T, rev; throwerrors=false)
                if spec === nothing
                    continue
                end
                for s in spec
                    @test s.type ∈ (:A, :C, :D, :E, :I, :J, :L)
                    S = fieldtype(T, s.symb)
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
                            @test S === Vector{String}
                        else
                            @test rank == ndims(S)
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

#=
@testset "Low-level Interface" begin
    # Test string comparisons/conversion "à la FITS".
    @test OIFITS.to_upper('a') == 'A'
    @test OIFITS.to_upper('b') == 'B'
    @test OIFITS.to_upper('z') == 'Z'
    @test OIFITS.to_upper('A') == 'A'
    @test OIFITS.to_upper('B') == 'B'
    @test OIFITS.to_upper('Z') == 'Z'
    @test OIFITS.to_upper('é') == 'é'
    @test OIFITS.is_same("Beta  ", "BETA") == true
    @test OIFITS.is_same(" Beta  ", "BETA") == false
    @test OIFITS.is_same("beta", " Beta  ") == false
    @test OIFITS.is_same("beta  ", "BeTa  ") == true
    @test OIFITS.fix_name("beTa  ") == "BETA"

    fits =  FITS(joinpath(dir, first(files)), "r")

    hdr = read_header(fits[1])
    @test OIFITS.get_hdu_type(hdr) === :image_hdu

    hdr["STRVAL"] = "a string"
    hdr["INTVAL"] = 3
    hdr["FLTVAL"] = 1.6

    @test_throws KeyError OIFITS.get_key_index(hdr, "DUMMY")
    @test OIFITS.get_key_index(hdr, "DUMMY", π) === π
    @test OIFITS.get_key_index(hdr, "SIMPLE", π) === 1

    @test_throws ErrorException OIFITS.get_integer(hdr, "STRVAL")
    @test isa(OIFITS.get_integer(hdr, "SIMPLE", :none), Int)
    @test isa(OIFITS.get_integer(hdr, "NAXIS"), Int)
    @test isa(OIFITS.get_integer(hdr, "NAXIS", :none), Int)

    @test_throws ErrorException OIFITS.get_logical(hdr, "NAXIS", true)
    @test OIFITS.get_logical(hdr, "SIMPLE") === true
    @test isa(OIFITS.get_logical(hdr, "SIMPLE", :none), Bool)

    @test_throws ErrorException OIFITS.get_float(hdr, "STRVAL", 1.0)
    @test isa(OIFITS.get_float(hdr, "NAXIS"), Float64)
    @test isa(OIFITS.get_float(hdr, "NAXIS", :none), Float64)

    @test_throws ErrorException OIFITS.get_logical(hdr, "DUMMY")
    @test OIFITS.get_logical(hdr, "DUMMY", nothing) === nothing
    @test typeof(OIFITS.get_logical(hdr, "SIMPLE")) === Bool
    @test typeof(OIFITS.get_logical(hdr, "SIMPLE", nothing)) === Bool

    @test_throws ErrorException OIFITS.get_comment(hdr, "DUMMY")
    @test OIFITS.get_comment(hdr, "DUMMY", nothing) === nothing
    @test typeof(OIFITS.get_comment(hdr, "SIMPLE")) === String
    @test typeof(OIFITS.get_comment(hdr, "SIMPLE", nothing)) === String

    hdr = read_header(fits[2])
    @test OIFITS.get_hdu_type(hdr) in (:image_hdu, :binary_table, :ascii_table)

    @test_throws ErrorException OIFITS.get_value(hdr, "DUMMY")
    @test OIFITS.get_value(hdr, "DUMMY", nothing) === nothing
    @test typeof(OIFITS.get_value(hdr, "XTENSION")) === String
    @test typeof(OIFITS.get_value(hdr, "XTENSION", nothing)) === String
    @test typeof(OIFITS.get_value(hdr, "BITPIX")) === Int
    @test typeof(OIFITS.get_value(hdr, "BITPIX", nothing)) === Int

    @test_throws ErrorException OIFITS.get_integer(hdr, "DUMMY")
    @test OIFITS.get_integer(hdr, "DUMMY", nothing) === nothing
    @test typeof(OIFITS.get_integer(hdr, "NAXIS")) === Int
    @test typeof(OIFITS.get_integer(hdr, "NAXIS", nothing)) === Int

    @test_throws ErrorException OIFITS.get_float(hdr, "DUMMY")
    @test OIFITS.get_float(hdr, "DUMMY", nothing) === nothing
    #@test typeof(OIFITS.get_float(hdr, "BSCALE")) === Float64
    @test typeof(OIFITS.get_float(hdr, "BSCALE", 1.0)) === Float64

    @test_throws ErrorException OIFITS.get_string(hdr, "DUMMY")
    @test_throws ErrorException OIFITS.get_string(hdr, "NAXIS")
    @test OIFITS.get_string(hdr, "DUMMY", nothing) === nothing
    @test typeof(OIFITS.get_string(hdr, "XTENSION")) === String
    @test typeof(OIFITS.get_string(hdr, "XTENSION", nothing)) === String

    # Test OIFITS.to_logical
    @test OIFITS.to_logical(true) === true
    @test OIFITS.to_logical(0) === false
    for tmp in ([true], (true,false))
        @test OIFITS.to_logical(tmp) === tmp
    end
    @test OIFITS.to_logical([0]) == [false]
    @test OIFITS.to_logical((1,0)) == (true,false)

    # Test OIFITS.to_integer
    @test OIFITS.to_integer(1) === 1
    @test OIFITS.to_integer(0x01) === 1
    for tmp in ([7], (11,))
        @test OIFITS.to_integer(tmp) === tmp
    end
    @test OIFITS.to_integer([0x05,0x0003]) == [5,3]
    @test OIFITS.to_integer((0x05,0x0003)) === (5,3)

    # Test OIFITS.to_float
    @test OIFITS.to_float(1.5) === 1.5
    @test OIFITS.to_float(1) === 1.0
    for tmp in ([7.1], (11.3,4.7))
        @test OIFITS.to_float(tmp) === tmp
    end
    @test OIFITS.to_float([0x05,0x0003]) == [5.0,3.0]
    @test OIFITS.to_float((0x05,0x0003)) === (5.0,3.0)

    # Test OIFITS.to_complex
    @test OIFITS.to_complex(1.5 + 3.1im) === 1.5 + 3.1im
    @test OIFITS.to_complex(1 - 2im) === 1.0 - 2.0im
    for tmp in ([1.3 - 7.1im], (11.3 + 0.0im, 4.7im))
        @test OIFITS.to_complex(tmp) === tmp
    end
    @test OIFITS.to_complex([0x05,0x0003]) == [5.0+0.0im,3.0+0.0im]
    @test OIFITS.to_complex((0x05,0x0003)) === (5.0+0.0im,3.0+0.0im)

    # Test OIFITS.to_string
    for tmp in ("hello", ["so", "wonderful"], ("world", ))
        @test OIFITS.to_string(tmp) === tmp
    end
    mesg = "hello world!"
    @test OIFITS.to_string(:hello) == "hello"
    @test OIFITS.to_string(SubString(mesg, 7:11)) == "world"
    @test OIFITS.to_string([SubString(mesg, 7:11),
                            SubString(mesg, 1:5)]) == ["world", "hello"]
    @test OIFITS.to_string([:hello,:world]) == ["hello","world"]
    @test OIFITS.to_string((:hello,SubString(mesg, 7:11))) == ("hello","world")

    # Test OIFITS.is_logical
    @test OIFITS.is_logical(false) == true
    @test OIFITS.is_logical(1) == false
    @test OIFITS.is_logical([false]) == true
    @test OIFITS.is_logical((false,)) == true

    # Test OIFITS.is_integer
    @test OIFITS.is_integer(0x33) == true
    @test OIFITS.is_integer(nothing) == false
    @test OIFITS.is_integer([0x12]) == true
    @test OIFITS.is_integer((1,0xff)) == true

    # Test OIFITS.is_float
    @test OIFITS.is_float(0x33) == true
    @test OIFITS.is_float(1im) == false
    @test OIFITS.is_float([0.1]) == true
    @test OIFITS.is_float((1,0xff)) == true

    # Test OIFITS.is_complex
    @test OIFITS.is_complex(1 + 0im) == true
    @test OIFITS.is_complex(1im) == true
    @test OIFITS.is_complex(1.2im) == true
    @test OIFITS.is_complex("hello") == false
    @test OIFITS.is_complex([0,1im]) == true
    @test OIFITS.is_complex((1,0+1im)) == false
    @test OIFITS.is_complex((3.0im,0+1im)) == true

    # Test OIFITS.is_string
    @test OIFITS.is_string(SubString(mesg, 7:11)) == true
    @test OIFITS.is_string(:hello) == false
    @test OIFITS.is_string(["hello"]) == true
    @test OIFITS.is_string(("a",)) == true

    # Check OIFITS.get_hdu_type
    @test OIFITS.get_hdu_type("BINTABLE") === :binary_table
    @test OIFITS.get_hdu_type("IMAGE") === :image_hdu
    @test OIFITS.get_hdu_type("TABLE") === :ascii_table
    @test OIFITS.get_hdu_type(TableHDU) === :binary_table
    @test OIFITS.get_hdu_type(ImageHDU) === :image_hdu
    @test OIFITS.get_hdu_type(ASCIITableHDU) === :ascii_table
    @test OIFITS.get_hdu_type(3) === :unknown
    @test OIFITS.get_hdu_type(fits[1]) === :image_hdu

    # Check OIFITS.to_fieldname
    @test OIFITS.to_fieldname("TIME_INT") === :time_int
    @test OIFITS.to_fieldname(:eff_wave) === :eff_wave

    # Check OIFITS.message
    @test OIFITS.warn("This is a test.") === nothing
    @test OIFITS.inform("This is another test.") === nothing
    @test OIFITS.message("This is yet another test.") === nothing

    # Check OIFITS.type_to_bitpix
    @test OIFITS.type_to_bitpix(UInt8)   ==   8
    @test OIFITS.type_to_bitpix(Int16)   ==  16
    @test OIFITS.type_to_bitpix(Int32)   ==  32
    @test OIFITS.type_to_bitpix(Int64)   ==  64
    @test OIFITS.type_to_bitpix(Float32) == -32
    @test OIFITS.type_to_bitpix(Float64) == -64

    # Check OIFITS.bitpix_to_type
    @test OIFITS.bitpix_to_type(  8) === UInt8
    @test OIFITS.bitpix_to_type( 16) === Int16
    @test OIFITS.bitpix_to_type( 32) === Int32
    @test OIFITS.bitpix_to_type( 64) === Int64
    @test OIFITS.bitpix_to_type(-32) === Float32
    @test OIFITS.bitpix_to_type(-64) === Float64
    @test_throws ErrorException OIFITS.bitpix_to_type(0)

    # Check OIFITS.eqcoltype_to_type and OIFITS.type_to_eqcoltype
    types = Vector{DataType}(undef, 200)
    for id in eachindex(types)
        types[id] = OIFITS.eqcoltype_to_type(id)
    end
    for T in (UInt8, Int8, Bool, String, Cushort, Cshort, Cuint, Cint,
              Culong, Clong, Cfloat, Int64, Cdouble,
              Complex{Cfloat}, Complex{Cdouble})
        @test types[OIFITS.type_to_eqcoltype(T)] === T
    end

    for i in 2:length(fits)
        if OIFITS.get_hdu_type(fits[i]) == :binary_table
            try
                extname, comment = read_key(fits[i], "EXTNAME")
                if extname == "OI_TARGET"
                    dict1 = OIFITS.read_table(fits[i])
                    dict2 = OIFITS.read_table(fits[i],
                                              key -> key ∈ ("RAEP0", "DECEP0"))
                    @test haskey(dict1, "TARGET") == true
                    @test haskey(dict2, "TARGET") == false
                    @test haskey(dict2, "RAEP0") == true
                    @test haskey(dict2, "RAEP0.units") == true
                    @test dict2["RAEP0.units"] == "deg"
                    break
                end
            catch ex
                nothing
            end
        end
    end
end
=#

@testset "Read OI FITS files" begin
    for file in files
        @test tryload(dir, file) == true
    end
end

#=
@testset "High-level Interface" begin
    master = OIFITS.load(joinpath(dir, first(files)))
    @test eltype(master) <: OIDataBlock
    @test size(master) == (length(master),)
    @test axes(master) === (Base.OneTo(length(master)),)
    @test eachindex(master) === Base.OneTo(length(master))
    for (key,db) in master.array
        @test isa(db, OIArray)
        @test db.extname == "OI_ARRAY"
    end
    for (key,db) in master.instr
        @test isa(db, OIWavelength)
        @test db.extname == "OI_WAVELENGTH"
        @test typeof(db.eff_wave) <: Vector
    end
    for (key,db) in master.correl
        @test isa(db, OICorrelation)
        @test db.extname == "OI_CORREL"
    end
    let db = master.target
        @test isa(db, OITarget)
        @test db.extname == "OI_TARGET"
    end

    # Select first target by name and by id.
    @test isa(OIFITS.select_target(master, first(master.target.target)),
              OIMaster)
    @test isa(OIFITS.select_target(master, first(master.target.target_id)),
              OIMaster)
    @test OIFITS.select_target(master, "") === nothing

    # Select first target by wavelength range.
    @test isa(OIFITS.select_wavelength(master, 0, Inf), OIMaster)
    subset = OIFITS.select_wavelength(master, w -> false)
    @test isa(subset, OIMaster)
    for db in subset
        @test isa(subset, Union{OIWavelength,OIVis,OIVis2,OIT3,
                                OIFlux,OIPolarization}) == false
    end
    for db in master
        if isa(db, Union{OIWavelength,OIVis,OIVis2,OIT3,
                         OIFlux,OIPolarization})
            @test OIFITS.select_wavelength(db, w -> false) === nothing
        end
    end

    # Test selection by extension names.
    sel = ("OI_VIS", "OI_VIS2", "OI_T3")
    for db in OIFITS.select(master, sel...)
        @test db.extname in sel
        # Check that `extname` is read-only.
        extname = db.extname
        @test isa(extname, String)
        @test_throws ErrorException db.extname = ""
        @test db.extname === extname
    end

    # Accessors (we check that the same object is returned).
    @test OIFITS.get_target(master) === master.target
    @test isa(OIFITS.get_targets(master), Vector{<:AbstractString})

    @test OIFITS.get_instrument(master) === master.instr
    @test isa(OIFITS.get_instruments(master), Vector{<:AbstractString})
    let db = OIFITS.get_instrument(master, first(OIFITS.get_instruments(master)))
        @test isa(db, OIWavelength)
        @test OIFITS.get_instrument(db) === db
    end

    @test OIFITS.get_array(master) === master.array
    @test isa(OIFITS.get_arrays(master), Vector{<:AbstractString})
    let db = OIFITS.get_array(master, first(OIFITS.get_arrays(master)))
        @test isa(db, OIArray)
        @test OIFITS.get_array(db) === db
    end

    for db in master
        @test OIFITS.get_extname(db) == db.extname
        @test OIFITS.Builder.last_revision(db) ≥ db.revn
        if isa(db, Union{OIVis,OIVis2,OIT3,OIFlux,OIPolarization})
            @test OIFITS.get_eff_wave(db) === db.instr.eff_wave
            @test OIFITS.get_eff_band(db) === db.instr.eff_band
            @test OIFITS.get_instrument(db) === db.instr
            @test OIFITS.get_array(db) === db.array
        end
        if isa(db, OITarget)
            @test OIFITS.get_target(db) === db.target
        elseif isa(db, Union{OIVis,OIVis2,OIT3,OIFlux,OIPolarization})
            @test OIFITS.get_target(db) === db.owner.target
        else
            @test OIFITS.get_target(db) === nothing
        end
        if isa(db, OIWavelength)
            @test OIFITS.get_instrument(db) === db
        elseif isa(db, Union{OIVis,OIVis2,OIT3,OIFlux,OIPolarization})
            @test OIFITS.get_instrument(db) === db.instr
            @test isa(OIFITS.get_instrument(db), OIWavelength)
        else
            @test OIFITS.get_instrument(db) === nothing
        end
        if isa(db, OIArray)
            @test OIFITS.get_array(db) === db
        elseif isa(db, Union{OIVis,OIVis2,OIT3,OIFlux,OIPolarization})
            @test OIFITS.get_array(db) === db.array
        else
            @test OIFITS.get_array(db) === nothing
        end
    end

    nwaves = 5
    instr1 = OIWavelength(eff_wave=rand(nwaves), eff_band=rand(nwaves),
                          insname="Instr1")
    @test isa(instr1, OIWavelength)
    @test isa(instr1, OIWavelength{Float64})
    @test OIFITS.float_type(instr1) === Float64
    @test instr1.revn == OIFITS.Builder.last_revision(instr1)
    @test OIFITS.is_attached(instr1) == false
    @test iswritable(instr1, :extname) == false
    @test iswritable(instr1, :insname) == true
    @test isreadonly(instr1, :extname) == true
    @test isreadonly(instr1, :insname) == false
    @test instr1.extname == "OI_WAVELENGTH"
    instr2 = OIWavelength{Float32}(eff_wave=rand(nwaves), eff_band=rand(nwaves),
                                   insname="Instr2")
    @test isa(instr2, OIWavelength)
    @test isa(instr2, OIWavelength{Float32})
    @test OIFITS.float_type(instr2) === Float32
    @test instr2.revn == OIFITS.Builder.last_revision(instr2)
    @test OIFITS.is_attached(instr2) == false

    owner = OIMaster(instr1, instr2)
    @test OIFITS.float_type(owner) === Float64
    @test "INSTR1" in keys(owner.instr)
    @test "INSTR2" in keys(owner.instr)
    @test OIFITS.is_attached(instr1) == true
    @test OIFITS.is_attached(instr1, owner) == true
    @test OIFITS.is_attached(instr2) == false

end
=#

end # module

nothing
