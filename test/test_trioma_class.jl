@testset "update_attribute! and inspect" begin

  # TriomaClass (from TriomaCore) is the canonical inspectable/updatable
  # container used to exercise inspect and update_attribute!.

  # --------------------------------------------------------------------------
  @testset "inspect - all attributes" begin
    obj = TriomaClass(a=100, b=200)
    out = IOBuffer()
    inspect(obj; io=out)
    s = String(take!(out))
    @test occursin("a: 100", s)
    @test occursin("b: 200", s)
  end

  @testset "inspect - specific attribute" begin
    obj = TriomaClass(a=100, b=200)
    out = IOBuffer()
    inspect(obj; variable_name="a", io=out)
    s = String(take!(out))
    @test  occursin("a: 100", s)
    @test !occursin("b: 200", s)
  end

  @testset "inspect - case insensitive" begin
    obj = TriomaClass(a=100)
    out = IOBuffer()
    inspect(obj; variable_name="A", io=out)
    s = String(take!(out))
    @test occursin("a: 100", s)
  end

  @testset "inspect - nested struct" begin
    child  = TriomaClass(value=999)
    parent = TriomaClass(value=50, child=child)
    out = IOBuffer()
    inspect(parent; io=out)
    s = String(take!(out))
    @test occursin("value: 50",  s)
    @test occursin("child",      s)
    @test occursin("value: 999", s)
  end

  # --------------------------------------------------------------------------
  @testset "update_attribute! - direct field" begin
    obj = TriomaClass(a=100)
    update_attribute!(obj, "a", 999)
    @test obj.a == 999
  end

  @testset "update_attribute! - nested field" begin
    child  = TriomaClass(a=100)
    parent = TriomaClass(child=child)
    update_attribute!(parent, "a", 999)
    @test parent.child.a == 999
  end

  @testset "update_attribute! - nonexistent raises ArgumentError" begin
    obj = TriomaClass(a=100)
    @test_throws ArgumentError update_attribute!(obj, "nonexistent", 999)
    try
      update_attribute!(obj, "nonexistent", 999)
    catch e
      @test occursin("nonexistent", e.msg)
    end
  end

  @testset "update_attribute! - n_pipes propagates to nested objects" begin
    child1 = TriomaClass(n_pipes=5)
    child2 = TriomaClass(n_pipes=5)
    parent = TriomaClass(n_pipes=5, child1=child1, child2=child2)
    update_attribute!(parent, "n_pipes", 10)
    @test parent.n_pipes == 10
    @test parent.child1.n_pipes == 10
    @test parent.child2.n_pipes == 10
  end

  @testset "update_attribute! - multiple nesting levels" begin
    leaf   = TriomaClass(value=100)
    middle = TriomaClass(leaf=leaf)
    root   = TriomaClass(middle=middle)
    update_attribute!(root, "value", 999)
    @test root.middle.leaf.value == 999
  end

end
