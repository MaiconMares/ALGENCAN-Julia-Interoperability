mutable struct MyDataPtr
  counters::Vector{Int64}

  MyDataPtr(counters) = new(zeros(5))
end

function main()
  using Libdl
  lib = Libdl.dlopen("./libfortran_struct.so")

  pdata = MyDataPtr([])
  pdata_ptr = Ref(pdata)

  ccall(
    (:__fortran_struct_MOD_inc_counter, "libfortran_struct.so"), Cvoid,
    (Ptr{MyDataPtr},), pdata_ptr)

  println(pdata.counters)

  Libdl.dlclose(lib)
end