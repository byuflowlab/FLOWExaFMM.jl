using Test
import FLOWExaFMM

# Hello world test
println("Running Hello, World test...")
@test FLOWExaFMM.greet()=="hello, world"
