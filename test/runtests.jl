module TestOIFITS

using OIFITS
using Compat
using Compat.Test
using Compat: @debug, @error, @info, @warn

dir = dirname(@__FILE__)

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

@testset "OIFITS.load" begin
    for file in files
        @test tryload(dir, file) == true
    end
end

end

nothing
