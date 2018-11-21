decode(x::Float32) = (b=bitstring(x); (b[1], b[2:9], b[10:32]))
decode(x::Float64) = (b=bitstring(x); (b[1], b[2:12], b[13:64]))

function ex1()
    a = Float32(1/3)
    b = Float64(1/3)
    c = Float64(a)
    
    println("\tliczba dziesiętnie\tznak\tcecha\t\tmantysa")
    println("a\t", a, "\t\t", decode(a)[1],"\t",decode(a)[2],"\t",decode(a)[3])
    println("b\t", b, "\t", decode(b)[1],"\t",decode(b)[2],"\t",decode(b)[3])
    println("c\t", c, "\t", decode(c)[1],"\t",decode(c)[2],"\t",decode(c)[3])
end

ex1()
