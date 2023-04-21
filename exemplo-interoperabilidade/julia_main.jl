function sum_callback(x::Float64, y::Float64)::Float64
  println("Callback called!")
  return x + y
end

function main()
  lib = Libdl.dlopen("./libfortran_client.so")
  # compute_media_ = Libdl.dlsym(lib, :compute_media_)

  a::Float64 = 10.0
  b::Float64 = 20.0

  sum_callback_ptr = @cfunction(sum_callback, Float64, (Ref{Float64}, Ref{Float64}))

  res = ccall((:__fortran_client_MOD_compute_media, "libfortran_client.so"), Float64,
              (Ptr{Cvoid}, Ref{Float64}, Ref{Float64}), sum_callback_ptr, a, b)

  # Or just res = ccall(compute_media_, Float64, (Ref{Float64}, Ref{Float64}), a, b)

  println(res)
  Libdl.dlclose(lib)

end