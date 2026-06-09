@testset "Project.toml" begin
  using TOML

  project = TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))

  @testset "required fields present" begin
    @test haskey(project, "package") || haskey(project, "name")
    @test haskey(project, "deps")
    @test haskey(project, "compat")
  end

  @testset "AtomicAndPhysicalConstants is a dependency" begin
    @test haskey(project["deps"], "AtomicAndPhysicalConstants")
  end

  @testset "core dependencies present" begin
    deps = project["deps"]
    @test haskey(deps, "Optim")
    @test haskey(deps, "QuadGK")
    @test haskey(deps, "LambertW")
  end

  @testset "julia compat declared" begin
    @test haskey(project["compat"], "julia")
  end

  @testset "no hardcoded physical constants in source" begin
    src_dir = joinpath(@__DIR__, "..", "src")
    patterns = [r"8\.617", r"8\.314", r"6\.022", r"6\.02214", r"1\.602"]
    for (root, dirs, files) in walkdir(src_dir)
      for fname in files
        endswith(fname, ".jl") || continue
        content = read(joinpath(root, fname), String)
        for pat in patterns
          @test !occursin(pat, content) broken=false
        end
      end
    end
  end

end
