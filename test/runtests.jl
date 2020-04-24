module TestOIFITS

using Test, Printf
using OIFITS, FITSIO

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
        db = OIFITS.load(joinpath(dir, file))
        counter += 1
        quiet || @info "file \"", file, "\" successfully loaded"
        return true
    catch
        @warn "failed to load \"", file, "\""
        return false
    end
end

@testset "Load OI-FITS" begin
    for file in files
        @test tryload(dir, file) == true
    end
end

@testset "Low-level Interface" begin
    fits =  FITS(joinpath(dir, first(files)), "r")
    hdr = read_header(fits[1])

    @test_throws ErrorException OIFITS.get_logical(hdr, "DUMMY")
    @test OIFITS.get_logical(hdr, "DUMMY", nothing) === nothing
    @test typeof(OIFITS.get_logical(hdr, "SIMPLE")) === Bool
    @test typeof(OIFITS.get_logical(hdr, "SIMPLE", nothing)) === Bool

    @test_throws ErrorException OIFITS.get_comment(hdr, "DUMMY")
    @test OIFITS.get_comment(hdr, "DUMMY", nothing) === nothing
    @test typeof(OIFITS.get_comment(hdr, "SIMPLE")) === String
    @test typeof(OIFITS.get_comment(hdr, "SIMPLE", nothing)) === String

    hdr = read_header(fits[2])

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
    @test OIFITS.get_string(hdr, "DUMMY", nothing) === nothing
    @test typeof(OIFITS.get_string(hdr, "XTENSION")) === String
    @test typeof(OIFITS.get_string(hdr, "XTENSION", nothing)) === String

end

@testset "High-level Interface" begin
    master = OIFITS.load(joinpath(dir, first(files)))
    @test eltype(master) <: OIDataBlock
    @test size(master) == (length(master),)
    @test axes(master) === (Base.OneTo(length(master)),)
    @test eachindex(master) === Base.OneTo(length(master))
    for (key,val) in master.array
        @test val.extname == "OI_ARRAY"
    end
    for (key,val) in master.instr
        @test val.extname == "OI_WAVELENGTH"
        @test typeof(val.eff_wave) <: Vector
    end
    for (key,val) in master.correl
        @test val.extname == "OI_CORREL"
    end
    #@test master.instr[first(master.instr)].extname == "OI_WAVELENGTH"
    @test master.target.extname == "OI_TARGET"

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
                                OISpectrum,OIPolarization}) == false
    end
    for db in master
        if isa(db, Union{OIWavelength,OIVis,OIVis2,OIT3,
                         OISpectrum,OIPolarization})
            @test OIFITS.select_wavelength(db, w -> false) === nothing
        end
    end

    # Test selection by extension names.
    sel = ("OI_VIS", "OI_VIS2", "OI_T3")
    for db in OIFITS.select(master, sel...)
        @test db.extname in sel
    end

end

end # module

nothing
