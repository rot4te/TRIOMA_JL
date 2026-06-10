@testset "update_attribute! and inspect" begin

  # --------------------------------------------------------------------------
  # Helpers: concrete mutable struct for testing inspect and update_attribute!
  # --------------------------------------------------------------------------
  mutable struct TestObj
    a   ::Union{Int, Nothing}
    b   ::Union{Int, Nothing}
    n_pipes::Union{Int, Nothing}
    child ::Union{TestObj, Nothing}
    child1::Union{TestObj, Nothing}
    child2::Union{TestObj, Nothing}
    value ::Union{Int, Nothing}
    leaf  ::Union{TestObj, Nothing}
    middle::Union{TestObj, Nothing}
  end
  TestObj(; a=nothing, b=nothing, n_pipes=nothing,
            child=nothing, child1=nothing, child2=nothing,
            value=nothing, leaf=nothing, middle=nothing) =
    TestObj(a, b, n_pipes, child, child1, child2, value, leaf, middle)

  # --------------------------------------------------------------------------
  @testset "inspect - all attributes" begin
    obj = TestObj(a=100, b=200)
    out = IOBuffer()
    inspect(obj; io=out)
    s = String(take!(out))
    @test occursin("a: 100", s)
    @test occursin("b: 200", s)
  end

  @testset "inspect - specific attribute" begin
    obj = TestObj(a=100, b=200)
    out = IOBuffer()
    inspect(obj; variable_name="a", io=out)
    s = String(take!(out))
    @test  occursin("a: 100", s)
    @test !occursin("b: 200", s)
  end

  @testset "inspect - case insensitive" begin
    obj = TestObj(a=100)
    out = IOBuffer()
    inspect(obj; variable_name="A", io=out)
    s = String(take!(out))
    @test occursin("a: 100", s)
  end

  @testset "inspect - nested struct" begin
    child  = TestObj(value=999)
    parent = TestObj(value=50, child=child)
    out = IOBuffer()
    inspect(parent; io=out)
    s = String(take!(out))
    @test occursin("value: 50",  s)
    @test occursin("child",      s)
    @test occursin("value: 999", s)
  end

  # --------------------------------------------------------------------------
  @testset "update_attribute! - direct field" begin
    obj = TestObj(a=100)
    update_attribute!(obj, "a", 999)
    @test obj.a == 999
  end

  @testset "update_attribute! - nested field" begin
    child  = TestObj(a=100)
    parent = TestObj(child=child)
    update_attribute!(parent, "a", 999)
    @test parent.child.a == 999
  end

  @testset "update_attribute! - nonexistent raises ArgumentError" begin
    obj = TestObj(a=100)
    @test_throws ArgumentError update_attribute!(obj, "nonexistent", 999)
    try
      update_attribute!(obj, "nonexistent", 999)
    catch e
      @test occursin("nonexistent", e.msg)
    end
  end

  @testset "update_attribute! - n_pipes propagates to nested objects" begin
    child1 = TestObj(n_pipes=5)
    child2 = TestObj(n_pipes=5)
    parent = TestObj(n_pipes=5, child1=child1, child2=child2)
    update_attribute!(parent, "n_pipes", 10)
    @test parent.n_pipes == 10
    @test parent.child1.n_pipes == 10
    @test parent.child2.n_pipes == 10
  end

  @testset "update_attribute! - multiple nesting levels" begin
    leaf   = TestObj(value=100)
    middle = TestObj(leaf=leaf)
    root   = TestObj(middle=middle)
    update_attribute!(root, "value", 999)
    @test root.middle.leaf.value == 999
  end

end
