module TestOIFITS

using OIFITS, Test, Printf

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

@testset "Interface" begin
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
end

end # module

nothing
