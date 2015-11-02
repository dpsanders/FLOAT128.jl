
function timeit(fn::Function, arg::Real, n::Int)
   mn = 9999999.9; mx=md=0.0;
   for i in 1:n
      a = time_ns()
      b = fn(arg)
      c = time_ns()-a
      if (c < mn) mn = c end;
      if (c > mx) mx = c end;
      md = md+c;
   end
   md1 = div(md,n)
   Int(mn),Int(md1),Int(mx)
end

